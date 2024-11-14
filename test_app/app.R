library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(mrgsolve)
library(dplyr)
library(purrr)
library(rhandsontable)

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

# Define species models and parameters
species_models <- list(
  Rat = list(model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
  Mouse = list(model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
  Human = list(model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
  Monkey = list(model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
)


# Parameters for PK models
default_params <- list(
  HLe_invivo = 6,       # Half-life in hours
  Vd = 50,              # Volume of distribution (L/kg)
  kabs = 0.5            # Absorption rate constant (1/h)
)

# PK Models
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


ui <- dashboardPage(
  dashboardHeader(title = "PFAS PK Model Simulation"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Experiment Inputs", tabName = "inputs", icon = icon("flask")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Model Parameters", tabName = "parameters", icon = icon("sliders-h"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "inputs",
              fluidRow(
                box(title = "Interactive Experiment Table", status = "primary", solidHeader = TRUE, width = 12,
                    p("Enter values for each species, PFAS chemical, model type, body weight, dose, interval, duration, and serum collection details."),
                    rHandsontableOutput("experiment_table"),
                    actionButton("add_row", "Add Row")
                ),
                actionButton("run", "Run Simulation")
              )
      ),
      
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
      ),
      
      tabItem(tabName = "parameters",
              fluidRow(
                uiOutput("species_params")  # Dynamically generate UI for each species' parameters
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Initialize editable table data with unrestricted input fields
  initial_table_data <- data.frame(
    PFAS = factor(c("PFHxS", "PFHxA", "PFOA", "PFOS"),
                  levels = c("PFHxS", "PFHxA", "PFOS",
                             "PFOA", "PFNA", "PFHpA")),
    Species = factor(c("Rat", "Mouse", "Rat", "Mouse"),
                     levels = c("Rat", "Mouse", "Human", "Monkey")),
    Model_Type = factor(c("PBPK", "single-compartment", "two-compartment", "PBPK")),
    Body_Weight_kg = c(0.3, 0.025, 0.3, 0.025),
    Dose_mg_per_kg = c(100, 100, 100, 100),
    Interval_Hours = as.integer(c(12, 12, 12, 12)),
    Exposure_Duration_Days = as.integer(c(3, 3, 3, 3)),
    Time_Serum_Collected_hr = as.integer(c(24, 48, 24, 48)),
    Serum_Concentration_mg_L = as.numeric(c(NA, NA, NA, NA)),
    stringsAsFactors = FALSE
  )
  
  # Reactive value to store experiment table data
  experiment_data <- reactiveVal(initial_table_data)
  
  # Render editable datatable
  output$experiment_table <- renderRHandsontable({
    rhandsontable(experiment_data()) %>% 
      hot_cols(strict = T,
               allowInvalid = F)
  })
  
  # Update experiment data based on user edits
  observeEvent(input$experiment_table_cell_edit, {
    info <- input$experiment_table_cell_edit
    data <- experiment_data()
    data[info$row, info$col + 1] <- as.numeric(info$value)
    experiment_data(data)
  })
  
  # Add a new row to the datatable when the "Add Row" button is clicked
  observeEvent(input$add_row, {
    new_row <- data.frame(
      PFAS = "", Species = "", Model_Type = "", Body_Weight_kg = NA, Dose_mg_per_kg = NA,
      Interval_Hours = NA, Exposure_Duration_Days = NA, Time_Serum_Collected_hr = NA,
      Serum_Concentration_mg_L = NA, stringsAsFactors = FALSE
    )
    experiment_data(rbind(experiment_data(), new_row))
  })
  
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
    
    exp_data %>%
      mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
        subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
      })) %>%
      select(PFAS, Species, Dose_mg_per_kg, Time_Serum_Collected_hr, Serum_Concentration_mg_L, modeled_concentration, Model_Type) %>% 
      datatable()
  })
  
  # Corrected renderUI for species parameters
  output$species_params <- renderUI({
    print("Rendering species parameters...")
    print(input$species)  # Debugging output
    
    # Manually check if input$species is valid
    if (is.null(input$species) || length(input$species) == 0) {
      return(NULL)
    }
    
    # If input$species is valid, create the UI elements
    lapply(input$species, function(sp) {
      print(paste("Creating parameters for species:", sp))  # Debugging output
      box(
        title = paste(sp, "Model Parameters"), status = "primary", solidHeader = TRUE, width = 12,
        DTOutput(paste0("params_", sp))
      )
    })
  })
  
  observe({
    for (sp in input$species) {
      local({
        species <- sp
        output[[paste0("params_", species)]] <- renderDT({
          datatable(reactive_params[[species]], editable = TRUE, options = list(pageLength = 10, autoWidth = TRUE))
        })
      })
    }
  })
  
  # Update the reactive parameters when the DataTable is edited
  observe({
    for (sp in input$species) {
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

