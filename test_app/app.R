################### LIBRARIES ############
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(mrgsolve)
library(dplyr)
library(purrr)
library(rhandsontable)

############### STATIC DATA ############

# Load PBPK models for different species
RatPBPK.code <- readRDS("../models/ratPBPK.RDS")
ratpbpk <- mcode("Ratpbpk", RatPBPK.code)

MicePBPK.code <- readRDS("../models/micePBPK.RDS")
micepbpk <- mcode("Micepbpk", MicePBPK.code)

humanPBPK.code <- readRDS("../models/humanPBPK.RDS")
humanpbpk <- mcode("Humanpbpk", humanPBPK.code)

monkeyPBPK.code <- readRDS("../models/monkeyPBPK.RDS")
monkeypbpk <- mcode("Monkeypbpk", monkeyPBPK.code)

# Load best fit parameters
rat_params <- readRDS("../models/rat_mcmc.rds")
rat_best <- rat_params$bestpar
mouse_params <- readRDS("../models/mouse_mcmc.rds")
mouse_best <- mouse_params$bestpar
human_params <- readRDS("../models/human_mcmc.rds")
human_best <- human_params$bestpar
monkey_params <- readRDS("../models/monkey_mcmc.rds")
monkey_best <- monkey_params$bestpar

# general TK params for non PBPK-models
tk_params <- readRDS("../Additional files/Datasets/general TK/tk_params.rds")

# # Define species models and parameters
# species_models <- list(
#   Rat = list(model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
#   Mouse = list(model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
#   Human = list(model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
#   Monkey = list(model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
# )

# Define species parameters lookup
species_models <- list(
  Rat = list(best_params = rat_best, default_bw = 0.3),
  Mouse = list(best_params = mouse_best, default_bw = 0.025),
  Human = list(best_params = human_best, default_bw = 82.3),
  Monkey = list(best_params = monkey_best, default_bw = 3.5)
)


# Parameters for PK models
default_params <- list(
  HLe_invivo = 6,       # Half-life in hours
  Vd = 50,              # Volume of distribution (L/kg)
  kabs = 0.5            # Absorption rate constant (1/h)
)

# Initialize editable table data with unrestricted input fields
initial_table_data <- data.frame(
  PFAS = factor(c("PFHxS", "PFHxA", "PFOA", "PFOS"),
                levels = c("PFHxS", "PFHxA", "PFOS",
                           "PFOA", "PFNA", "PFHpA")),
  Species = factor(c("Rat", "Mouse", "Rat", "Mouse"),
                   levels = c("Rat", "Mouse", "Human", "Monkey")),
  Sex = factor(c("Female", "Male")),
  Route = factor(c("Intravenous", "Intragastric", "Oral", "Oral") ,
                 levels = c("Intravenous", "Intragastric", "Intraperitoneal", "Oral", "Dermal")),
  Model_Type = factor(c("PBPK", "single-compartment", "two-compartment", "PBPK")),
  Body_Weight_kg = c(0.3, 0.025, 0.3, 0.025),
  Dose_mg_per_kg = c(100, 100, 100, 100),
  Interval_Hours = as.integer(c(12, 12, 12, 12)),
  Exposure_Duration_Days = as.integer(c(3, 3, 3, 3)),
  Time_Serum_Collected_hr = as.integer(c(24, 48, 24, 48)),
  Serum_Concentration_mg_L = as.numeric(c(NA, NA, NA, NA)),
  stringsAsFactors = FALSE
)

# Capture the original column types and levels
original_types <- sapply(initial_table_data, class)
original_levels <- lapply(initial_table_data, levels)

# Reactive value to store experiment table data
experiment_data <- reactiveVal(initial_table_data)

############# PK Model Functions ############
simplified_conc <- function(D_per_dose, Vd, ke, tau, n, t, exposure_duration_hr) {
  if (t <= exposure_duration_hr) {
    doses_given <- floor(t / tau)
    C <- (D_per_dose / Vd) * ((1 - exp(-doses_given * ke * tau)) / (1 - exp(-ke * tau))) * exp(-ke * (t %% tau))
  } else {
    time_since_last_dose <- t - exposure_duration_hr
    C <- (D_per_dose / Vd) * exp(-ke * time_since_last_dose)
  }
  return(C)
}

full_conc <- function(D_per_dose, Vd, ke, ka, tau, n, t, exposure_duration_hr) {
  if (t <= exposure_duration_hr) {
    doses_given <- floor(t / tau)
    C <- (D_per_dose * ka) / (Vd * (ka - ke)) * (
      (exp(-ke * (t %% tau)) / (1 - exp(-ke * tau))) - 
        (exp(-ka * (t %% tau)) / (1 - exp(-ka * tau)))
    )
  } else {
    time_since_last_dose <- t - exposure_duration_hr
    C <- (D_per_dose * ka) / (Vd * (ka - ke)) * (
      (exp(-ke * time_since_last_dose)) - 
        (exp(-ka * time_since_last_dose))
    )
  }
  return(C)
}

########################## USER INTERFACE ################

ui <- dashboardPage(
  dashboardHeader(title = "PFAS PK Model Simulation"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Experiment Inputs", tabName = "inputs", icon = icon("flask")),
      menuItem("Model Parameters", tabName = "parameters", icon = icon("sliders-h")),
      menuItem("Results", tabName = "results", icon = icon("chart-line"))
      
    )
  ),
  
  dashboardBody(
    tabItems(
      ############ Input Tab Item #####
      tabItem(tabName = "inputs",
              fluidRow(
                box(title = "Input Experimental Conditions", status = "primary", solidHeader = TRUE, width = 12,
                    p("Enter values for each species, PFAS chemical, model type, body weight, dose, interval, duration, and serum collection details."),
                    rHandsontableOutput("experiment_table")#,
                    #actionButton("add_row", "Add Row")
                ),
                actionButton("run", "Run Simulation",
                             style = "color: #fff; background-color: #077336; border-color: #2e6da4; font-size: 20px; padding: 10px 15px;")
              )
      ),
      ############ Parameters Tab Item #####
      tabItem(tabName = "parameters",
              box(status = "primary", width = 12, collapsible = T,
              fluidRow(
                tabBox(width = 12,
                       tabPanel(title = "Simple TK Models",
                                rHandsontableOutput("params_table")  # Dynamically generate UI for each species' parameters
                                ),
                       tabPanel(title = "PBPK Models",
                                uiOutput("species_params")  # Dynamically generate UI for each species' parameters
                       ),
                       )
                )
              )
              ),
      ############ Results Tab Item #####
      tabItem(tabName = "results",
              fluidRow(
                box(title = "Concentration Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("concentrationPlot")),
                box(title = "Modeled vs. Measured Concentration", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("scatterPlot"))
              ),
              fluidRow(
                h3("Modeled and Measured Concentrations"),
                DTOutput("concentration_at_time")
              )
      )
      )
    )
  )

################# SERVER ############

server <- function(input, output, session) {

######## Experiment data table #####  
  # Render editable datatable
  output$experiment_table <- renderRHandsontable({
    rhandsontable(experiment_data()) %>% 
      hot_cols(strict = T,
               allowInvalid = F)
  })
  
  # Observe changes in the rhandsontable and update reactive value
  observeEvent(input$experiment_table, {
    # Get the updated data from the input
    updated_data <- hot_to_r(input$experiment_table)
    
    # Reapply original data types and levels
    for (col in names(updated_data)) {
      if (original_types[col] == "factor") {
        # Convert back to factor with original levels
        updated_data[[col]] <- factor(updated_data[[col]], levels = original_levels[[col]], ordered = is.ordered(initial_table_data[[col]]))
      } else if (original_types[col] == "integer") {
        # Convert back to integer
        updated_data[[col]] <- as.integer(updated_data[[col]])
      } else if (original_types[col] == "numeric") {
        # Convert back to numeric
        updated_data[[col]] <- as.numeric(updated_data[[col]])
      } else if (original_types[col] == "character") {
        # Convert back to character
        updated_data[[col]] <- as.character(updated_data[[col]])
      }
    }
    # Update the reactive value with the new data
    experiment_data(updated_data)
  })


############ Run Simulation #####
  # Run simulation based on experiment data
  simulation_results <- eventReactive(input$run, {
    results <- lapply(1:nrow(experiment_data()), function(i) {
      row <- experiment_data()[i,]
      species <- row$Species
      model_type <- row$Model_Type
      bw <- row$Body_Weight_kg
      dose <- row$Dose_mg_per_kg
      interval <- row$Interval_Hours
      exposure_duration <- row$Exposure_Duration_Days * 24
      
      if (model_type == "PBPK") {
        model <- species_models[[species]]$model
        ex <- ev(ID = 1, amt = dose * bw, ii = interval, addl = floor(exposure_duration / interval) - 1, cmt = "AST")
        tgrid <- tgrid(0, exposure_duration + 96, 1)
        output <- model %>%
          update(atol = 1E-8, maxsteps = 10000) %>%
          mrgsim_d(data = ex, tgrid = tgrid)
        data.frame(Species = species, PFAS = row$PFAS, Time = output$time / 24, Concentration = output$Plasma)
      } else {
        ke <- log(2) / default_params$HLe_invivo
        ka <- if (model_type == "two-compartment") default_params$kabs else NULL
        n_doses <- floor(exposure_duration / interval)
        times <- seq(0, exposure_duration + 96, by = 1)
        conc <- if (!is.null(ka)) {
          map_dbl(times, ~ full_conc(dose, default_params$Vd, ke, ka, interval, n_doses, .x, exposure_duration))
        } else {
          map_dbl(times, ~ simplified_conc(dose, default_params$Vd, ke, interval, n_doses, .x, exposure_duration))
        }
        data.frame(Species = species, PFAS = row$PFAS, Time = times / 24, Concentration = conc)
      }
    })
    do.call(rbind, results)
  })
  
  # Render concentration-time plot with dynamic xlim
  output$concentrationPlot <- renderPlotly({
    sim_results <- simulation_results()
    
    p <- ggplot(sim_results, aes(x = Time, y = Concentration, color = interaction(Species, PFAS))) +
      geom_line(size = 1) +
      labs(title = "Plasma Concentration Over Time",
           x = "Time (days)", 
           y = "Concentration (ug/ml)") +
      xlim(0, max(sim_results$Time, na.rm = TRUE) + 4) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Render scatter plot for modeled vs. measured concentrations
  output$scatterPlot <- renderPlotly({
    sim_results <- simulation_results()
    exp_data <- experiment_data() %>% filter(!is.na(Time_Serum_Collected_hr))
    
    measured_conc <- exp_data %>%
      mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
        subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
      })) %>%
      filter(!is.na(modeled_concentration))
    
    p <- ggplot(measured_conc, aes(x = modeled_concentration, y = Serum_Concentration_mg_L, color = Species)) +
      geom_point(size = 2) +
      labs(title = "Modeled vs. Measured Concentration",
           x = "Modeled Concentration (mg/L)",
           y = "Measured Concentration (mg/L)") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Display concentration at specified time points in a table
  output$concentration_at_time <- renderDT({
    sim_results <- simulation_results()
    exp_data <- experiment_data() %>% filter(!is.na(Time_Serum_Collected_hr))
    
    data <- exp_data %>%
      mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
        subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
      })) %>%
      select(PFAS, Species, Dose_mg_per_kg, Time_Serum_Collected_hr, Serum_Concentration_mg_L, modeled_concentration, Model_Type)
    
      dt <- datatable(data,
        rownames = F,
        extensions = 'Buttons', #enable buttons extension
        filter = "top",
        options = list(pageLength = 25, autoWidth = TRUE,  width = '100%', scrollX = TRUE,
                       dom = 'Blrtip', 
                       buttons = list(
                         # insert buttons with copy and print
                         # colvis includes the button to select and view only certain columns in the output table
                         # from https://rstudio.github.io/DT/extensions.html 
                         I('colvis'), 'copy', 
                         # code for the first dropdown download button. this will download only the current page only (depends on the number of rows selected in the lengthMenu)
                         # using modifier = list(page = "current")
                         # only the columns visible will be downloaded using the columns:":visible" option from:
                         list(extend = 'collection', buttons = list(list(extend = "csv", filename = "page",exportOptions = list(
                           columns = ":visible",modifier = list(page = "current"))),
                           list(extend = 'excel', filename = "page", title = NULL, 
                                exportOptions = list(columns = ":visible",modifier = list(page = "current")))),
                           text = 'Download current page'),
                         # code for the  second dropdown download button
                         # this will download the entire dataset using modifier = list(page = "all")
                         list(extend = 'collection',
                              buttons = list(list(extend = "csv", filename = "data",exportOptions = list(
                                columns = ":visible",modifier = list(page = "all"))),
                                list(extend = 'excel', filename = "data", title = NULL, 
                                     exportOptions = list(columns = ":visible",modifier = list(page = "all")))),
                              text = 'Download all data')),
                       # add the option to display more rows as a length menu
                       lengthMenu = list(c(10, 30, 50, -1),
                                         c('10', '30', '50', 'All'))),class = "display"
        ) %>% 
      formatStyle(
        target = 'row',
        backgroundColor = styleEqual(
          c("PBPK", "single-compartment", "two-compartment"), 
          c("#3A9AB2", "#6FB2C1", "#91BAB6")
        ),
        columns = "Model_Type"  # Specifies to base coloring on the data_available column
      ) 
    
    # Identify numeric columns
    numeric_cols <- names(data)[sapply(data, is.numeric)]
    
    # Apply formatting for significant digits to numeric columns
    for (col in numeric_cols) {
      dt <- dt %>% formatSignif(columns = col, digits = 3)
    }
    
     dt   
  })
  
  
############ Parameters Data Table #####
  # initalize input params
  initial_params_data <- tk_params
  
  # Reactive value to store params table data
  params_data <- reactiveVal(initial_params_data)
  
  # params editable table
  output$params_table <- renderRHandsontable({
    # read in experimental data to know what to serve up
    exp_data <- experiment_data() %>% 
      filter(Model_Type != "PBPK") %>%  # just get the simple TK data for this table
      distinct(PFAS, Species, Sex, Route, Model_Type) %>% 
      mutate(selected = "selected")
    
    # import params data from PFHpA repo
    params_data <- exp_data  %>% 
      left_join(params_data(),
                by = c("PFAS", "Species", "Sex", "Route")) %>% 
      filter(selected == "selected") %>% 
      select(-selected) %>% 
      arrange(desc(PFAS), Species, Sex)
    
    
    # JavaScript renderer function to set background color
    color_renderer <- "
  function (instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.TextRenderer.apply(this, arguments);
    if (cellProperties.readOnly) {
      td.style.background = '#D3D3D3';  // Light gray for read-only columns
    } else {
      td.style.background = '#FFFF99';  // Light yellow for editable columns
    }
  }
"
    rhandsontable(params_data) %>% 
      hot_cols(strict = T,
               allowInvalid = F) %>% 
      hot_col("PFAS", readOnly = T) %>% 
      hot_col("Sex", readOnly = T) %>% 
      hot_col("Species", readOnly = T) %>% 
      hot_col("Strain", readOnly = T) %>% 
      hot_col("Route", readOnly = T) %>% 
      hot_col("Model_Type", readOnly = T) %>% 
      hot_cols(renderer = color_renderer)  # Apply the custom renderer to all columns
  })
  
###### PBPK Params Table ####
  # Reactive values for model parameters and body weights
  reactive_params <- reactiveValues()
  reactive_bw <- reactiveValues(
    Rat = 0.3,
    Mouse = 0.025,
    Human = 82.3,
    Monkey = 3.5
  )
  
  # Observe experiment table and update UI based on conditions
  observe({
    table_data <- experiment_data()
    matching_rows <- subset(table_data, 
                            PFAS == "PFOS" & 
                              Sex == "Male" & 
                              Route == "Oral" & 
                              Model_Type == "PBPK")
    
    for (sp in unique(matching_rows$Species)) {
      if (!is.null(species_models[[sp]]) && is.null(reactive_params[[sp]])) {
        params <- species_models[[sp]]$best_params
        reactive_params[[sp]] <- data.frame(
          Parameter = names(params),
          Value = signif(sapply(params, exp), 4)
        )
        reactive_bw[[sp]] <- species_models[[sp]]$default_bw
      }
    }
  })
  
  # Render UI for body weight inputs
  output$body_weight_inputs <- renderUI({
    #req(input$experiment_table)
    table_data <- experiment_data()
    matching_rows <- subset(table_data, 
                            PFAS == "PFOS" & 
                              Sex == "Male" & 
                              Route == "Oral" & 
                              Model_Type == "PBPK")
    if (nrow(matching_rows) == 0) return(NULL)
    
    lapply(unique(matching_rows$Species), function(sp) {
      numericInput(
        inputId = paste0("bw_", sp), 
        label = paste(sp, "Body Weight (kg)"), 
        value = reactive_bw[[sp]]
      )
    })
  })
  
  # Render UI for species parameters
  output$species_params <- renderUI({
  #  req(input$experiment_table)
    table_data <- experiment_data()
    matching_rows <- subset(table_data, 
                            PFAS == "PFOS" & 
                              Sex == "Male" & 
                              Route == "Oral" & 
                              Model_Type == "PBPK")
    if (nrow(matching_rows) == 0) return(NULL)
    
    lapply(unique(matching_rows$Species), function(sp) {
      box(
        title = paste(sp, "Model Parameters"), status = "primary", solidHeader = TRUE, width = 12,
        DTOutput(paste0("params_", sp))
      )
    })
  })
  
  # Render and update species-specific parameter tables
  observe({
    for (sp in names(species_models)) {
      local({
        species <- sp
        output[[paste0("params_", species)]] <- renderDT({
        #  req(reactive_params[[species]])
          datatable(
            reactive_params[[species]], 
            editable = TRUE, 
            options = list(pageLength = 10, autoWidth = TRUE),
            rownames = FALSE
          )
        })
      })
    }
  })
  
  # Update reactive parameters when the DataTable is edited
  observe({
    for (sp in names(species_models)) {
      local({
        species <- sp
        observeEvent(input[[paste0("params_", species, "_cell_edit")]], {
          info <- input[[paste0("params_", species, "_cell_edit")]]
          params <- reactive_params[[species]]
          params[info$row, info$col] <- as.numeric(info$value)
          reactive_params[[species]] <- params
        })
      })
    }
  })

  
} # close server

shinyApp(ui = ui, server = server)

