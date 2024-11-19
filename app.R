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
library(readr)
library(readxl)
library(tidyverse)
library(cols4all)
library(shinyjs)


############### STATIC DATA ############

# Load PBPK models for different species
RatPBPK.code <- readRDS("models/ratPBPK.RDS")
ratpbpk <- mcode("Ratpbpk", RatPBPK.code)

MicePBPK.code <- readRDS("models/micePBPK.RDS")
micepbpk <- mcode("Micepbpk", MicePBPK.code)

humanPBPK.code <- readRDS("models/humanPBPK.RDS")
humanpbpk <- mcode("Humanpbpk", humanPBPK.code)

monkeyPBPK.code <- readRDS("models/monkeyPBPK.RDS")
monkeypbpk <- mcode("Monkeypbpk", monkeyPBPK.code)

# Load best fit parameters
rat_params <- readRDS("models/rat_mcmc.rds")
rat_best <- rat_params$bestpar
mouse_params <- readRDS("models/mouse_mcmc.rds")
mouse_best <- mouse_params$bestpar
human_params <- readRDS("models/human_mcmc.rds")
human_best <- human_params$bestpar
monkey_params <- readRDS("models/monkey_mcmc.rds")
monkey_best <- monkey_params$bestpar

# general TK params for non PBPK-models
tk_params <- readRDS("Additional files/Datasets/general TK/tk_params.rds")

# Define species models and parameters
species_models <- list(
  Rat = list(model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
  Mouse = list(model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
  Human = list(model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
  Monkey = list(model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
)

# # Define species parameters lookup
# species_models <- list(
#   Rat = list(best_params = rat_best, default_bw = 0.3),
#   Mouse = list(best_params = mouse_best, default_bw = 0.025),
#   Human = list(best_params = human_best, default_bw = 82.3),
#   Monkey = list(best_params = monkey_best, default_bw = 3.5)
# )


# Initialize editable table data with unrestricted input fields
initial_table_data <- data.frame(
  PFAS = factor(c("PFBA", "PFBA", "PFHxA", "PFOA", "PFOS"),
                levels = unique(tk_params$PFAS)),
  Species = factor(c("Rat", "Mouse", "Mouse", "Rat", "Rat"),
                   levels = c("Rat", "Mouse", "Human", "Monkey")),
  Sex = factor(c("Female", "Male", "Female", "Female", "Male")),
  # Route = factor(c("Intravenous", "Intragastric", "Oral", "Oral") ,
  #                levels = c("Intravenous", "Intragastric", "Intraperitoneal", "Oral", "Dermal")),
  Model_Type = factor(c("single-compartment", "two-compartment", "single-compartment", "two-compartment", "PBPK")),
  #Body_Weight_kg = c(0.3, 0.025, 0.3, 0.3),
  Dose_mg_per_kg = c(175, 35, 62.5, 50, 2.5),
  Interval_Hours = as.integer(c(24, 24, 24, 24, 24)),
  Exposure_Duration_Days = as.integer(c(17, 28, 28, 28, 28)),
  Time_Serum_Collected_hr = as.integer(c(432, 696, 696, 696, 696)),
  Serum_Concentration_mg_L = as.numeric(c(4.44, 86, 3.1, 9.326, 173.7)),
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

# Function to calculate AUC using the trapezoidal rule
calculate_auc <- function(time, concentration) {
  stopifnot(length(time) == length(concentration))
  time <- sort(time)
  
  auc <- sum((concentration[-1] + concentration[-length(concentration)]) / 2 * diff(time))
  return(auc)
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
                box(title = "Experimental Conditions", status = "primary", solidHeader = TRUE, width = 12,
                    p("Enter values for each species, PFAS chemical, model type, body weight, dose, interval, duration, and serum collection details."),
                    p("Use the interactive table below or upload a file to proceed."),
                    rHandsontableOutput("experiment_table"),
                    br(),
                    p("Optional: upload data below:"),
                    fileInput("upload_experiment_data", "Upload Experiment Data (.csv or .xlsx)", accept = c(".csv", ".xlsx"))
                ),
                br(),
                uiOutput("params_check_message"), # Add dynamic message below button
                br(),
                useShinyjs(),
                div(
                  style = "text-align: center;",
                  actionButton(
                    "run", 
                    label = HTML('<i class="fa fa-rocket"></i> Run Simulation'),
                    style = "color: #fff; background-color: #077336; border-color: #2e6da4; font-size: 20px; padding: 10px 15px;"
                  )
                )
              )
      ),
      ############ Parameters Tab Item #####
      tabItem(tabName = "parameters",
              box(status = "primary", width = 12, collapsible = T,
              fluidRow(
                tabBox(width = 12,
                       tabPanel(title = "Simple TK Models",
                                rHandsontableOutput("params_table"),  # Dynamically generate UI for each species' parameters
                                br(),
                                p("Optionally upload TK params data:"),
                                fileInput("upload_params_data", "Upload Parameters Data (.csv or .xlsx)", accept = c(".csv", ".xlsx"))
                                ),
                       tabPanel(title = "PBPK Models",
                                uiOutput("species_sex_params_grid")  # Dynamically generate UI for each species' parameters
                       ),
                       )
                )
              )
              ),
      ############ Results Tab Item #####
      tabItem(tabName = "results",
              box(status = "primary", width = 12, collapsible = T,
              fluidRow(
                tabBox(width = 12,
                       tabPanel(title = "Time-Series Results",
                                box(title = "Concentration Plot", status = "primary", solidHeader = TRUE, width = 12,
                                    plotlyOutput("concentrationPlot"),
                                    p("Download modeled time-series dataset:"),
                                    downloadButton("download_time_conc", "Download Concentration Data")
                                    ),
                                box(title = "Modeled vs. Measured Concentration", status = "primary", solidHeader = TRUE, width = 12,
                        plotlyOutput("scatterPlot")),
                        fluidRow(
                          h3("Modeled and Measured Concentrations"),
                          DTOutput("concentration_at_time")
                        ),
                        ),
                       tabPanel(title = "Summary Results",
                                fluidRow(
                                  h3("Summary Values"),
                                  p("Values below for maximum concentration (C_max), area under the curve (AUC), and time-weighted average concentraiton (C_TWA) are in mg/L for serum"),
                                  DTOutput("summary_table"),
                                  plotlyOutput("heatmap")
                                  )
                                )
                       )
                )
              )
              )
    )
  )
)


################# SERVER ############

server <- function(input, output, session) {
  
  
#################################################### Data entry and Upload #######################################
######## Experiment data table #####  
  # Reactive value to store experiment table data
  experiment_data <- reactiveVal(initial_table_data)
  
  # Upload and update experiment_data
  observeEvent(input$upload_experiment_data, {
    req(input$upload_experiment_data)
    file <- input$upload_experiment_data
    ext <- tools::file_ext(file$name)
    if (ext == "csv") {
      new_data <- read_csv(file$datapath, show_col_types = FALSE)
    } else if (ext == "xlsx") {
      new_data <- read_excel(file$datapath)
    } else {
      showNotification("Invalid file type. Please upload a .csv or .xlsx file.", type = "error")
      return(NULL)
    }
    
    # Validate and update experiment_data
    if (!all(c("PFAS", "Species", "Sex", 
               #"Route", 
               "Model_Type"
               #, "Body_Weight_kg"
               ) %in% colnames(new_data))) {
      showNotification("Invalid data format. Missing required columns.", type = "error")
      return(NULL)
    }
    experiment_data(new_data)
    showNotification("Experiment data successfully updated.", type = "message")
  })
  
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
  
  ########################### Parameters Data Table ############
  # initalize input params
  initial_params_data <- tk_params
  
  # Reactive value to store params table data
  params_data <- reactiveVal(initial_params_data)
  
  # Reactive value to store processed params table
  processed_params <- reactiveVal(NULL)
  
  # Observe changes in experiment_data() and update processed_params
  observeEvent(experiment_data(), {
    # Filter experimental data to non-PBPK rows
    exp_data <- experiment_data() %>%
      filter(Model_Type != "PBPK") %>%  # Exclude PBPK rows for this table
      distinct(PFAS, Species, Sex,
               #Route, 
               Model_Type) %>%
      mutate(selected = "selected")
    
    # Perform the join with tk_params
    params_processed <- exp_data %>%
      left_join(tk_params, by = c("PFAS", "Species", "Sex"
                                  #, "Route"
                                  )
                ) %>%
      filter(selected == "selected") %>%  # Keep only selected rows
      select(-selected) %>%  # Drop the helper column
      distinct() %>% 
      arrange(desc(PFAS), Species, Sex)  # Arrange for consistent display
    
    # Debugging: Inspect the processed data
    print("Processed params table:")
    print(params_processed)
    
    # Update reactive value
    processed_params(params_processed)
  })
  
  # Upload and update params_data (model parameters for species/sex combinations)
  observeEvent(input$upload_params_data, {
    shiny::req(input$upload_params_data)
    file <- input$upload_params_data
    ext <- tools::file_ext(file$name)
    if (ext == "csv") {
      new_params <- read_csv(file$datapath, show_col_types = FALSE)
    } else if (ext == "xlsx") {
      new_params <- read_excel(file$datapath)
    } else {
      showNotification("Invalid file type. Please upload a .csv or .xlsx file.", type = "error")
      return(NULL)
    }

    # Update reactive_params for each species/sex combination
    for (row in 1:nrow(new_params)) {
      species <- new_params[row, "Species"]
      sex <- new_params[row, "Sex"]
      params <- new_params[row, ]

      if (!is.null(species) && !is.null(sex)) {
        reactive_params[[paste0(species, "_", sex)]] <- params
      }
    }
    showNotification("Parameters data successfully updated.", type = "message")
  })

  # Render the RHandsontable using processed_params
  output$params_table <- renderRHandsontable({
    shiny::req(processed_params())  # Ensure processed_params is available
    
    # JavaScript renderer function to format numbers to 4 significant digits
    custom_renderer <- "
    function (instance, td, row, col, prop, value, cellProperties) {
      Handsontable.renderers.TextRenderer.apply(this, arguments);
      if (!isNaN(value) && typeof value === 'number') {
        td.innerHTML = value.toPrecision(4);
      }
    }
  "
    
    # Create rhandsontable
    rhandsontable(processed_params()) %>%
      hot_cols(strict = TRUE, allowInvalid = FALSE) %>%
      hot_col("PFAS", readOnly = TRUE) %>%
      hot_col("Species", readOnly = TRUE) %>%
      hot_col("Sex", readOnly = TRUE) %>%
     # hot_col("Route", readOnly = TRUE) %>%
      hot_col("Model_Type", readOnly = TRUE) %>% 
      hot_cols(renderer = custom_renderer)  # Apply the custom renderer to all columns
  })
  
 
  # Observe changes in the rhandsontable and update reactive value
  observeEvent(input$params_table, {
    # Get the updated data from the input
    updated_data <- hot_to_r(input$params_table)

    # Update the reactive value with the new data
    params_data(updated_data)
  })
  
 
  ############################### PBPK Params Table #############################
  # Reactive values for PBPK model parameters
  reactive_params <- reactiveValues()
  
  # Observe experiment table and dynamically update parameter tables for each species/sex
  observeEvent(experiment_data(), {
    shiny::req(experiment_data())  # Ensure experiment_data() is available
    
    table_data <- experiment_data()
    pbpk_data <- subset(table_data, Model_Type == "PBPK")
    selected_pfas <- unique(pbpk_data$PFAS)
    
    # Loop through each unique Species/Sex combination in pbpk_data
    for (sp in unique(pbpk_data$Species)) {
      for (sx in unique(pbpk_data$Sex)) {
        # Filter data for current species and sex
        species_sex_data <- subset(pbpk_data, Species == sp & Sex == sx)
        params <- if (!is.null(species_models[[sp]])) species_models[[sp]]$params else NULL
        
        # Initialize parameter table if it doesn't already exist
        if (is.null(reactive_params[[paste0(sp, "_", sx)]])) {
          reactive_params[[paste0(sp, "_", sx)]] <- data.frame(
            Parameter = if (!is.null(params)) names(params) else character(0),
            stringsAsFactors = FALSE
          )
        }
        
        # Populate columns for each selected PFAS
        for (pfas in selected_pfas) {
          is_supported_model <- pfas == "PFOS" && 
            any(species_sex_data$PFAS == "PFOS" & species_sex_data$Sex == sx)
          
          if (is_supported_model && !is.null(params)) {
            reactive_params[[paste0(sp, "_", sx)]][[pfas]] <- signif(sapply(params, exp), 4)
          } else {
            reactive_params[[paste0(sp, "_", sx)]][[pfas]] <- NA
          }
        }
      }
    }
  })
  
  # Render UI for species/sex parameter tables in a grid layout
  output$species_sex_params_grid <- renderUI({
    shiny::req(experiment_data())
    table_data <- experiment_data()
    pbpk_data <- subset(table_data, Model_Type == "PBPK")
    
    if (nrow(pbpk_data) == 0) return(NULL)
    
    species_sex_combinations <- unique(paste(pbpk_data$Species, pbpk_data$Sex, sep = "_"))
    
    # Create a list of fluidRows with each row containing two columns
    ui_elements <- list()
    for (i in seq(1, length(species_sex_combinations), by = 2)) {
      row_elements <- list()
      
      for (j in 0:1) {
        index <- i + j
        if (index <= length(species_sex_combinations)) {
          species_sex <- species_sex_combinations[index]
          sp_sex_split <- strsplit(species_sex, "_")[[1]]
          species <- sp_sex_split[1]
          sex <- sp_sex_split[2]
          
          row_elements[[j + 1]] <- column(
            width = 6,
            box(
              title = paste(species, sex, "Model Parameters"), 
              status = "primary", 
              solidHeader = TRUE, 
              width = 12,
              DTOutput(paste0("params_", species, "_", sex))
            )
          )
        }
      }
      
      # Add each pair of species/sex tables to a new row
      ui_elements[[length(ui_elements) + 1]] <- fluidRow(row_elements)
    }
    
    # Return the complete list of UI elements for grid layout
    ui_elements
  })
  
  # Render and update species/sex-specific parameter tables
  observe({
    for (sp in unique(initial_table_data$Species)) {
      for (sx in unique(initial_table_data$Sex)) {
        local({
          species <- sp
          sex <- sx
          output[[paste0("params_", species, "_", sex)]] <- renderDT({
            shiny::req(reactive_params[[paste0(species, "_", sex)]])
            datatable(
              reactive_params[[paste0(species, "_", sex)]], 
              editable = TRUE, 
              options = list(pageLength = 10, autoWidth = TRUE),
              rownames = FALSE
            )
          })
        })
      }
    }
  })
  
  # Update reactive parameters when the DataTable is edited
  observe({
    for (sp in unique(initial_table_data$Species)) {
      for (sx in unique(initial_table_data$Sex)) {
        local({
          species <- sp
          sex <- sx
          observeEvent(input[[paste0("params_", species, "_", sex, "_cell_edit")]], {
            info <- input[[paste0("params_", species, "_", sex, "_cell_edit")]]
            params <- reactive_params[[paste0(species, "_", sex)]]
            params[info$row, info$col] <- as.numeric(info$value)
            reactive_params[[paste0(species, "_", sex)]] <- params
          })
        })
      }
    }
  })
  
  
################ error output #####
  validate_params <- function(experiment_data, params_data, reactive_params, required_params) {
    # Check if experiment_data or params_data is empty
    if (is.null(experiment_data) || nrow(experiment_data) == 0) {
      return("No experiment data available.")
    }
    if (is.null(params_data) || nrow(params_data) == 0) {
      return("No parameters data available. Please upload or provide parameter data to run the simulation.")
    }
    
    # Convert factors to characters
    experiment_data[] <- lapply(experiment_data, function(x) if (is.factor(x)) as.character(x) else x)
    params_data[] <- lapply(params_data, function(x) if (is.factor(x)) as.character(x) else x)
    
    # Identify missing parameters or cells
    missing_rows <- lapply(seq_len(nrow(experiment_data)), function(row_idx) {
      row <- experiment_data[row_idx, ]
      model_type <- row$Model_Type
      pfas <- row$PFAS
      species <- row$Species
      sex <- row$Sex
     # route <- row$Route
      required <- required_params[[model_type]]
      
      if (model_type == "PBPK") {
        # Check PBPK model parameters in reactive_params
        param_table <- reactive_params[[paste0(species, "_", sex)]]
        if (is.null(param_table) || !pfas %in% colnames(param_table)) {
          return(paste(species, sex, pfas, "(PBPK parameters missing)"))
        }
      } else {
        # Check required parameters for non-PBPK models in params_data
        matching_rows <- params_data[
          params_data$PFAS == pfas &
            params_data$Species == species &
            params_data$Sex == sex,# &
          #  params_data$Route == route,
        ]
        
        if (nrow(matching_rows) == 0) {
          return(paste(species, sex, pfas, #route, 
                       "(No parameter data)"))
        }
        
        # Check if all required parameters are available as columns
        missing <- required[!required %in% colnames(matching_rows)]
        if (length(missing) > 0) {
          return(paste(species, sex, pfas, 
                       #route, 
                       "(Missing columns:", paste(missing, collapse = ", "), ")"))
        }
        
        # Check for missing values in the required columns
        for (param in required) {
          if (any(is.na(matching_rows[[param]]))) {
            return(paste(species, sex, pfas,
                         #route,
                         "(Missing values in:", param, ")"))
          }
        }
      }
      return(NULL)
    })
    
    # Flatten and return missing rows
    missing_rows <- unlist(Filter(Negate(is.null), missing_rows))
    
    if (length(missing_rows) == 0) {
      return(NULL)  # No issues
    } else {
      # Combine missing rows with <br> for line breaks
      return(paste("Missing parameter data for:<br>", paste(missing_rows, collapse = "<br>")))
    }
  }
  

  # Define dynamic required parameters
  required_params <- reactive({
    list(
      "PBPK" = NULL,
      "single-compartment" = c("Clearance_L_per_kg_d", "Volume_of_Distribution_L_per_kg", "Half_Life_hr"),
      "two-compartment" = c("Clearance_L_per_kg_d", "Volume_of_Distribution_L_per_kg", "Half_Life_hr", "Absorption_Coefficient_unitless")
    )
  })
  
  # Reactive value to store the validation message
  validation_message <- reactiveVal(NULL)
  
  # Update validation message dynamically
  observeEvent(list(experiment_data(), params_data()), {
    print("experiment_data or params_data changed")
    shiny::req(experiment_data(), params_data())
    
    # Compute validation result
    validation_result <- validate_params(
      experiment_data = experiment_data(),
      params_data = params_data(),
      reactive_params = reactive_params,
      required_params = required_params()  # Use dynamic required_params
    )
    
    # Set the message dynamically
    if (is.null(validation_result)) {
      validation_message(
        tags$p(
          "All parameters are sufficiently available to run the simulation.",
          style = "color: green; font-size: 16px; text-align: center;"
        )
      )
    } else {
      validation_message(
        tags$p(
          HTML(validation_result),  # Use HTML to interpret line breaks
          style = "color: red; font-size: 16px; text-align: center;"
        )
      )
    }
  })
  
  
  # Render the validation message in the UI
  output$params_check_message <- renderUI({
    shiny::req(validation_message())
    validation_message()
  })
  
    # Ensure the output is updated even when the tab is not active
  outputOptions(output, "params_check_message", suspendWhenHidden = FALSE)
  
  
  ### Observe and Force UI Re-rendering on Data Changes
  observeEvent(list(experiment_data(), params_data()), {
    print("Triggered re-rendering due to data change")  # Debugging message
    output$params_check_message <- renderUI({
      shiny::req(experiment_data(), params_data())
      
      # Validate parameters
      validation_result <- validate_params(
        experiment_data = experiment_data(),
        params_data = params_data(),
        reactive_params = reactive_params,
        required_params = required_params()  # Use dynamic required_params
      )
      
      # Display the appropriate message
      if (is.null(validation_result)) {
        tags$p(
          "All parameters are sufficiently available to run the simulation.",
          style = "color: green; font-size: 16px; text-align: center;"
        )
      } else {
        tags$p(
          HTML(validation_result),  # Use HTML to interpret line breaks
          style = "color: red; font-size: 16px; text-align: center;"
        )
      }
    })
  })
  

 #########################################################################################################
############################################### Run Simulation ##########################################
 #####################################################################################################
  
  simulation_results <- eventReactive(input$run, {
    shiny::req(params_data(), experiment_data())
    
    print("Running simulation...")
    
    # Ensure consistent data types before joining
    exp_data <- experiment_data() %>%
      mutate(
        Species = as.character(Species),
        PFAS = as.character(PFAS),
        Sex = as.character(Sex)
      )
    
    params <- params_data() %>%
      mutate(
        Species = as.character(Species),
        PFAS = as.character(PFAS),
        Sex = as.character(Sex)
      )
    
    # Join experiment_data with params_data
    experiment_with_params <- exp_data %>%
      left_join(params, by = c("Species", "Sex", "PFAS"))
    
    # Define consistent column structure
    result_cols <- c(
      "Species", "PFAS", "Sex", "Model_Type", "Dose_mg_per_kg",
      "Exposure_Duration_Days", "Interval_Hours", "Time", "Concentration"
    )
    
  
    # Run the simulation for each row in the joined data
    results <- lapply(1:nrow(experiment_with_params), function(i) {
      row <- experiment_with_params[i, ]
      pfas <- row$PFAS
      species <- row$Species
      sex <- row$Sex
      model_type <- row$Model_Type
    #  bw <- row$Body_Weight_kg
      dose <- row$Dose_mg_per_kg
      interval <- row$Interval_Hours
      dose_frequency_per_day <- 24 / interval
      exposure_duration_days <- row$Exposure_Duration_Days
      exposure_duration_hrs <- exposure_duration_days * 24
      
      cat("Debugging inputs:\n")
      cat("Species:", species, "\n")
      cat("PFAS:", pfas, "\n")
      cat("Model Type:", as.character(model_type), "\n")
      cat("Dose per kg:", dose, "mg/kg\n")
      cat("Exposure Duration:", exposure_duration_hrs, "hours\n")
      
      if (model_type == "PBPK") {
        # Extract default kinetic parameters from species_models
        species_model <- species_models[[species]]
        pars <- species_model$params
        model <- species_model$model
        default_bw <- species_model$default_bw
        # For PBPK model, dose calculations
        dose_mg <- default_bw * dose
        interval_between_doses <- 24 / dose_frequency_per_day
        total_doses <- exposure_duration_days * dose_frequency_per_day
        
        cat("Dose:", dose_mg, "mg\n")
        cat("Body Weight:", default_bw, "kg\n")
        cat("Interval:", interval, "hours\n")
        
        ex <- ev(ID = 1, amt = dose_mg, ii = interval_between_doses, addl = total_doses - 1, cmt = "AST")
        tgrid <- tgrid(0, exposure_duration_hrs + 96, 1)
      
        
        output <- tryCatch({
          model %>%
            param(10 ^ pars) %>% #get out of log10 space
            Req(Plasma) %>%
            update(atol = 1E-8, maxsteps = 10000) %>%
            mrgsim_d(data = ex, tgrid = tgrid)
        }, error = function(e) {
          message(paste("Error for Species:", species, pfas, "Message:", e$message))
          return(NULL)
        })
        
        if (!is.null(output)) {
          data.frame(
            Species = species,
            PFAS = pfas,
            Sex = sex,
            Model_Type = model_type,
            Dose_mg_per_kg = dose,
            Exposure_Duration_Days = exposure_duration_days,
            Interval_Hours = interval,
            Time = output$time / 24,
            Concentration = signif(output$Plasma, 4)
          )
        } else {
          # Return NA for PBPK model failures
          data.frame(
            Species = species,
            PFAS = pfas,
            Sex = sex,
            Model_Type = model_type,
            Dose_mg_per_kg = dose,
            Exposure_Duration_Days = exposure_duration_days,
            Interval_Hours = interval,
            Time = NA,
            Concentration = NA
          )
        }
      } else {
        # For simplified model (one or two-compartment)
        ke <- log(2) / row$Half_Life_hr
        ka <- if (model_type == "two-compartment") row$Absorption_Coefficient_unitless else NULL
        vd <- row$Volume_of_Distribution_L_per_kg
        n_doses <- floor(exposure_duration_hrs / interval)
        times <- seq(0, exposure_duration_hrs + 96, by = 1)
        
        conc <- if (!is.null(ka)) {
          map_dbl(times, ~ full_conc(dose, vd, ke, ka, interval, n_doses, .x, exposure_duration_hrs))
        } else {
          map_dbl(times, ~ simplified_conc(dose, vd, ke, interval, n_doses, .x, exposure_duration_hrs))
        }
        
        data.frame(
          Species = species,
          PFAS = pfas,
          Sex = sex,
          Model_Type = model_type,
          Dose_mg_per_kg = dose,
          Exposure_Duration_Days = exposure_duration_days,
          Interval_Hours = interval,
          Time = times / 24,
          Concentration = conc
        )
      }
    })
    
    # Ensure all results have the same columns and structure
    results <- lapply(results, function(res) {
      res[result_cols[!result_cols %in% names(res)]] <- NA
      res[result_cols]
    })
    
    do.call(rbind, results)
  })
  
  
  
  
  
  
  

  ###################################################### Results ###########################################  
  ############################### Concentration-time Plot ##############
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
  
  ############################ Scatterplot of measured vs. modeled #############
  # Render scatter plot for modeled vs. measured concentrations
  output$scatterPlot <- renderPlotly({
    sim_results <- simulation_results()
    exp_data <- experiment_data() %>% filter(!is.na(Time_Serum_Collected_hr))
    
    measured_conc <- exp_data %>%
      mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
        subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
      })) %>%
      filter(!is.na(modeled_concentration))
    
    # Calculate dynamic limits based on the data range
    data_range <- range(
      c(measured_conc$modeled_concentration, measured_conc$Serum_Concentration_mg_L),
      na.rm = TRUE
    )
    
    # Expand the range slightly for better visualization
    expanded_limits <- c(
      10^(floor(log10(min(data_range)))),
      10^(ceiling(log10(max(data_range))))
    )
    
    #make scatterplot
    p <- ggplot(measured_conc, aes(x = modeled_concentration, 
                                   y = Serum_Concentration_mg_L,
                                   color = paste(PFAS, Species, #route,
                                                 sep = " - "),
                                   shape = Sex,
                                   text = paste("Chemical:", PFAS, "<br>",
                                                "Species:", Species, "<br>",
                                                "Sex:", Sex, "<br>",
                                                "Model:", Model_Type, "<br>",
                                                "Dose:", Dose_mg_per_kg, "mg/kg-d", "<br>",
                                                "Exposure Duration:", Exposure_Duration_Days, "days", "<br>",
                                                "Measured Serum Concentration:", Serum_Concentration_mg_L, "mg/L", "<br>",
                                                "Modeled Serum Concentration:", modeled_concentration, "mg/L", "<br>"))) +
      geom_point(size = 2, alpha = 0.9) +
      scale_x_log10(limits = expanded_limits) +
      scale_y_log10(limits = expanded_limits) +
      scale_color_discrete_c4a_cat(name = "") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") + # Add 10x and 0.1x deviation lines
      geom_abline(slope = 1, intercept = log10(10), linetype = "dotted", color = "gray70") +
      geom_abline(slope = 1, intercept = log10(0.1), linetype = "dotted", color = "gray70") +
      labs(title = "Modeled vs. Measured Concentration",
           x = "Modeled Concentration (mg/L)",
           y = "Measured Concentration (mg/L)") +
      theme_minimal() +
      theme(legend.title = element_blank())
    
    ggplotly(p, tooltip = "text")
  })
  
  ###### Excel Sheet for downloading time-concentration data
  output$download_time_conc <- downloadHandler(
    filename = function() {
      paste("time_series_concentration_data", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      # Load required library
      library(openxlsx)
      
      # load data
      sim_results <- simulation_results()
      
      # Prepare the data (assuming `sim_results` contains your time-series data)
      time_series_data <- sim_results %>%
        group_by(Species, Sex, PFAS, Model_Type) %>%
        group_split()
      
      # Create a new workbook
      wb <- createWorkbook()
      
      # Iterate over each group and add a sheet
      for (group_data in time_series_data) {
        # Define the sheet name (adjust as necessary to avoid overly long names)
        sheet_name <- paste0(
          group_data$Species[1], "-", 
          group_data$Sex[1], "-", 
          group_data$PFAS[1], "-", 
          group_data$Model_Type[1]
        )
        
        # Shorten sheet name if it exceeds Excel's limit
        sheet_name <- substr(sheet_name, 1, 31)
        
        # Add the data to the workbook
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, group_data)
      }
      
      # Save the workbook
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  
  ################################### Datatable of Time Point ###############
  # Display concentration at specified time points in a table
  output$concentration_at_time <- renderDT({
    sim_results <- simulation_results()
    exp_data <- experiment_data() %>% filter(!is.na(Time_Serum_Collected_hr))
    
    data <- exp_data %>%
      mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
        subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
      })) %>%
      select(PFAS, Species, Dose_mg_per_kg, Time_Serum_Collected_hr, Serum_Concentration_mg_L, modeled_concentration, Model_Type) #%>% 
      # left_join(summary_stats(),
      #           by = c("PFAS", "Species", "Sex", "Model_Type", "Dose_mg_per_kg", #"Interval_Hours",
      #                  "Exposure_Duration_Days"))
    
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

########################################## Calculate Summary Statistics ############  
  # Calculate summary statistics (C_max, C_TWA, AUC) for a specified time period
  summary_stats <- reactive({
    print("Calculating summary statistics...")
    
    # Manually check if the simulation results exist and are valid
    sim_results <- simulation_results()
    if (is.null(sim_results) || nrow(sim_results) == 0) {
      return(NULL)
    }
    
    sim_results %>%
      #get experimental input data
      left_join(experiment_data(),
                by = c("PFAS", "Species", "Sex", "Model_Type", "Dose_mg_per_kg", "Interval_Hours", "Exposure_Duration_Days")) %>% 
      group_by(Species, PFAS, Sex, Model_Type, Dose_mg_per_kg, Exposure_Duration_Days) %>%
      summarize(
        C_max = max(Concentration),
        AUC = {
          time_subset <- Time >= 0 & Time <= Exposure_Duration_Days
          calculate_auc(Time[time_subset], Concentration[time_subset])
        },
        C_TWA = AUC / (Exposure_Duration_Days - 0)
      ) %>%
      ungroup() %>% 
      distinct()
  })

  ##### Summary Stats Table #####
  # Render the summary statistics table with export options
  output$summary_table <- renderDT({
    print("Rendering summary table...")

    # Get the summary statistics
    summary <- summary_stats()

    # Manually check if the summary statistics exist and are valid
    if (is.null(summary) || nrow(summary) == 0) {
      return(NULL)
    }

    datatable(summary,
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
      ) %>%
      formatSignif(columns = c("C_max", "AUC", "C_TWA"), digits = 4)
  })

#### Summary stats heatmap ####
  output$heatmap <- renderPlotly({
    
    params_heat <- summary_stats() %>%
      pivot_longer(
        cols = c(C_max, C_TWA), # Replace these with the actual column names in `summary_table`
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      group_by(Metric) %>%
      mutate(
        normalized_value = (Value - min(Value, na.rm = TRUE)) /
          (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))
      ) %>%
      ungroup() %>%
      arrange(PFAS) %>%
      ggplot(aes(
        x = Metric,
        y = paste(PFAS, Sex, Species, Model_Type, sep = " - "), 
        fill = normalized_value + 0.001,
        text = paste0(
          "Chemical: ", PFAS, "<br>",
          "Species: ", Species, "<br>",
          "Model Type: ", Model_Type, "<br>",
          "Sex: ", Sex, "<br>",
          "Computed Metric: ", Metric, "<br>",
          "Value: ", Value, "<br>"
        )
      )) +
      geom_tile(color = "white") +
      scale_fill_gradient(
        trans = scales::log_trans(base = 10),
        low = "#56B1F7",
        high = "red4",
        space = "Lab",
        na.value = "grey50"
      ) +
      labs(
        x = "Toxicokinetic Parameters",
        y = "Chemical - Sex - Species - Model Type",
        title = "Comparison of Computed Metrics",
        subtitle = "Identical Exposure Duration (28 d) and Dose (1 mg/kg-d)"
      ) +
      theme_minimal(base_size = 15) +
      theme(
        axis.text.y = element_text(size = 14),  # Adjust the y-axis text size
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"
      )
    
    
    params_heat
  })
  
  
} # close server

shinyApp(ui = ui, server = server)

