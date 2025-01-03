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
library(shinythemes)


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
tk_params <- readRDS("Additional files/Datasets/general/tk_params.rds")

# Define species models and parameters
species_models <- list(
  Rat = list(PFAS = "PFOS", Sex = "Male", model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
  Mouse = list(PFAS = "PFOS", Sex = "Male",model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
  Human = list(PFAS = "PFOS", Sex = "Male",model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
  Monkey = list(PFAS = "PFOS", Sex = "Male",model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
)


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
# Function for single-compartment model
simplified_conc <- function(D_per_dose, Vd, ke, tau, n, t, exposure_duration_hr) {
  # Calculate concentration at each time point
  concentration <- numeric(length(t))
  for (i in 1:n) {
    dose_time <- (i - 1) * tau
    concentration <- concentration + (D_per_dose / Vd) * exp(-ke * (t - dose_time)) * (t >= dose_time)
  }
  return(concentration)
}

# Function for two-compartment model
full_conc <- function(D_per_dose, Vd, ke, ka, tau, n, t, exposure_duration_hr) {
  # Calculate concentration at each time point
  concentration <- numeric(length(t))
  # Calculate concentrations for each dose administered
  for (i in 1:n) {  # Loop over each dose (from 1 to n)
    dose_time <- (i - 1) * tau  # Calculate the time at which the current dose is administered
    for (j in 1:length(t)) {  # Loop over each time point in the defined sequence
      if (t[j] >= dose_time) {  # Check if the current time point is greater than or equal to the dose time
        if (ka != ke) {  # Ensure that the absorption rate constant is not equal to the elimination rate constant
          exp_ke <- exp(-ke * (t[j] - dose_time))  # Calculate the exponential term for elimination
          exp_ka <- exp(-ka * (t[j] - dose_time))  # Calculate the exponential term for absorption
          if (is.finite(exp_ke) && is.finite(exp_ka)) {  # Check if both exponential calculations are finite
            # Update the concentration at the current time point using the pharmacokinetic formula
            concentration[j] <- concentration[j] + (D_per_dose / Vd) * (ka / (ka - ke)) *
              (exp_ke - exp_ka)  # Add the contribution of the current dose to the concentration
          }
        } else {
        # If ka equals ke, handle it here
        concentration[j] <- concentration[j] + (D_per_dose / Vd) * ka * (t[j] - dose_time) * exp(-ke * (t[j] - dose_time))
      }
    }
    }
  }
  return(concentration)
}



##### ORIGINAL FORMULA ######
# simplified_conc <- function(D_per_dose, Vd, ke, tau, n, t, exposure_duration_hr) {
#   if (t == 0) {
#     # Explicitly set concentration to 0 at the start of the simulation
#     C <- 0
#   } else if (t <= exposure_duration_hr) {
#     # Calculate the cumulative concentration over time (transient phase)
#     doses_given <- floor(t / tau)
#     C <- 0  # Initialize concentration
#     for (i in 0:(doses_given - 1)) {
#       dose_time <- i * tau
#       C <- C + (D_per_dose / Vd) * exp(-ke * (t - dose_time))
#     }
#   } else {
#     # After the exposure duration, calculate residual concentration from all doses
#     doses_given <- floor(exposure_duration_hr / tau)
#     C <- 0  # Initialize concentration
#     for (i in 0:(doses_given - 1)) {
#       dose_time <- i * tau
#       C <- C + (D_per_dose / Vd) * exp(-ke * (t - dose_time))
#     }
#   }
#   return(C)
# }
# 
# full_conc <- function(D_per_dose, Vd, ke, ka, tau, n, t, exposure_duration_hr) {
#   if (t == 0) {
#     # Explicitly set concentration to 0 at the start of the simulation
#     C <- 0
#   } else if (t <= exposure_duration_hr) {
#     # Calculate the cumulative concentration over time (transient phase)
#     doses_given <- floor(t / tau)
#     C <- 0  # Initialize concentration
#     for (i in 0:(doses_given - 1)) {
#       dose_time <- i * tau
#       C <- C + (D_per_dose * ka) / (Vd * (ka - ke)) * (
#         exp(-ke * (t - dose_time)) - exp(-ka * (t - dose_time))
#       )
#     }
#   } else {
#     # After the exposure duration, calculate residual concentration from all doses
#     doses_given <- floor(exposure_duration_hr / tau)
#     C <- 0  # Initialize concentration
#     for (i in 0:(doses_given - 1)) {
#       dose_time <- i * tau
#       C <- C + (D_per_dose * ka) / (Vd * (ka - ke)) * (
#         exp(-ke * (t - dose_time)) - exp(-ka * (t - dose_time))
#       )
#     }
#   }
#   return(C)
# }

# Function to calculate AUC using the trapezoidal rule
calculate_auc <- function(time, concentration) {
  stopifnot(length(time) == length(concentration))
  time <- sort(time)
  
  auc <- sum((concentration[-1] + concentration[-length(concentration)]) / 2 * diff(time))
  return(auc)
}

########################## USER INTERFACE ################

ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "OEHHA PFAS PK Modelling Application", titleWidth = 400),
  dashboardSidebar(
    # Logo Image
    tags$img(src="main_logo_drop.png", width = "75%", height = "75%", style = 'display: block; margin-left: auto; margin-right: auto;'),
    tags$div("Logo Copyright OEHHA (2024)", align = 'center', style = 'font-size: 12px; display: block; margin-left: auto; margin-right: auto;'), 
    tags$div("Contact: Scott.Coffin@oehha.ca.gov", align = 'center', style = 'font-size: 12px; display: block; margin-left: auto; margin-right: auto;'), 
    br(),
    br(),
    sidebarMenu(
      menuItem("Experiment Inputs", tabName = "inputs", icon = icon("flask")),
      menuItem("Model Parameters", tabName = "parameters", icon = icon("sliders-h")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Additional OEHHA Apps", tabName = "additional_apps", icon = icon("cubes")),
      # GitHub Link at the Bottom
      div(style = "position: absolute; bottom: 12px; width: 100%;",
          menuItem(
            "GitHub Repository",
            href = "https://github.com/ScottCoffin/PFOS_PBPK_fork",
            icon = icon("github"),
            newtab = TRUE # Open link in a new tab
          )
      )
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      ############ Input Tab Item #####
      tabItem(tabName = "inputs",
              fluidRow(
                box(title = "Experimental Conditions", status = "primary", solidHeader = TRUE, width = 12,
                    p("Enter values for each species, PFAS chemical, model type, body weight, dose, interval, duration, and serum collection details."),
                    p(strong("Use the interactive table below or upload a file using the button below.")),
                    rHandsontableOutput("experiment_table"),
                    br(),
                    div(style = "display: flex; justify-content: center; align-items: center; height: 100%;",
                        box(status = "primary", width = 8,
                            div(style = "text-align: center;",
                                p(strong("Optionally")," you may upload experimental conditions data. Ensure that all of the above column names are present, and that all fields are valid (i.e., numeric values in numeric columns, factors in factor columns). A sample spreadsheet may be downloaded here that serves as a guide."),
                                downloadButton("example_experiment_data", "Download Example Spreadsheet"),
                                br(),
                                br(),
                                fileInput("upload_experiment_data", "Upload Experimental Conditions Table (.csv or .xlsx)", accept = c(".csv", ".xlsx"))
                            )
                        )
                    ),
                    p("Once experimental conditions data are entered above, please proceed by clicking on the", strong("Model Parameters tab"), "on the left side of the app.")
                )
              )
      ),
      ############ Parameters Tab Item #####
      tabItem(tabName = "parameters",
              fluidRow(
                box(title = "PK Model Parameter Customization", status = "primary", solidHeader = TRUE, width = 12,
                    tabBox(width = 12,
                       tabPanel(title = "Simple TK Models",
                                p("Toxicokinetic parameter data for PFAS are from a curated large variety of data sources. Whenever possible, TK parameters from authoritative sources (e.g., ATSDR, USEPA) are used when multiple are available for a given chemical-species-sex combination. The full, compiled TK dataset may be viewed in the interactive datatable below."),
                                rHandsontableOutput("params_table"),  # Dynamically generate UI for each species' parameters
                                p("Please note that `Absorption Coefficient` refers to the absorption phase kinetic coefficient (commonly referred to as `Ka`) - and is not the same as the absorption fraction (commonly referred to as `F`."),
                                br(),
                                p("Optionally upload TK params data:"),
                                fileInput("upload_params_data", "Upload Parameters Data (.csv or .xlsx)", accept = c(".csv", ".xlsx")),
                                br(),
                                uiOutput("params_check_message"), # Add dynamic message below button
                                useShinyjs(),
                                div(
                                  style = "text-align: center;",
                                  actionButton(
                                    "run", 
                                    label = HTML('<i class="fa fa-rocket"></i> Run Simulation'),
                                    style = "color: #fff; background-color: #077336; border-color: #2e6da4; font-size: 20px; padding: 10px 15px;"
                                  ),
                                  actionButton(
                                    "reset", 
                                    label = HTML('<i class="fa fa-refresh"></i> Reset'),
                                    style = "color: #fff; background-color: #d9534f; border-color: #d43f3a; font-size: 16px; padding: 10px 15px;"
                                  ),
                                  br(),
                                  br(),
                                  p("To view simulation results, go to the ", strong("Results"), "tab on the left side of the app."), 
                                ),
                                br(),
                                h3("Explore the full toxicokinetic dataset below:"),
                                box(title = "Full TK Dataset", width = 12, collapsible = T, collapsed = T,
                                    DTOutput("Full_TK_Datatable"),
                                    ),
                                ),
                       tabPanel(title = "PBPK Models",
                                p("This RShiny app allows users to run the optimized models from", a("Chou & Lin et al. (2019). Bayesian evaluation of a physiologically based pharmacokinetic (PBPK) model for perfluorooctane sulfonate (PFOS) to characterize the interspecies uncertainty between mice, rats, monkeys, and humans: Development and performance verification",
                                                                                                     href = "https://www.sciencedirect.com/science/article/pii/S016041201930203X")),
                                p("The model was only optimized for", strong("males"), "and via an oral route of exposure. The model consists of four organ compartments (plasma, liver, kidney, and rest of body). For more information, please refer to their publication. A schematic of the PBPK model framework from the publication is shown below:"),
                                tags$a(
                                  href = "https://www.sciencedirect.com/science/article/pii/S016041201930203X",
                                  target = "_blank",
                                  tags$img(src = "fig1.PNG", width = "70%", alt = "PBPK Framework")
                                ),
                                br(),
                                br(),
                                h3("Customize PBPK Parameters"),
                                p("This PBPK model is fully customizable (however not recommended)."),
                                uiOutput("species_sex_params_grid"),  # Dynamically generate UI for each species' parameters
                       )
                )
              )
              )
              ),
      ############ Results Tab Item #####
      tabItem(tabName = "results",
              box(title = "Results", status = "primary", width = 12, collapsible = T,solidHeader = TRUE,
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
              ),
      ################## Additional Apps Item #################
      tabItem(
        tabName = "additional_apps",
        fluidRow(
          box(
            title = "Explore Additional OEHHA Applications",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("OEHHA offers several interactive applications to explore and analyze environmental and public health data. Below, you can find links, descriptions, and GitHub repositories for some of these tools."),
            br(),
            
            # Application 1
            box(
              title = "OEHHA Chemical Data Explorer Application",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              p("The OEHHA Chemical Data Explorer Tool helps users rapidly identify, visualize, and access chemical hazard, production, and exposure data."),
              tags$a(
                href = "https://oehha.shinyapps.io/MasterChemicalList/", 
                target = "_blank", 
                tags$span(
                  icon("link"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff;",
                  " Visit the Chemical Data Explorer Tool"
                )
              ),
              tags$span(
                tags$a(
                  href = "https://github.com/ScottCoffin/MasterChemicalList", 
                  target = "_blank",
                  icon("github"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff; margin-left: 10px;",
                  " View GitHub Repository"
                )
              ),
              br(),
              tags$a(
                href = "https://oehha.shinyapps.io/MasterChemicalList/",
                target = "_blank",
                tags$img(src = "chemical_app_screenshot.png", width = "70%", alt = "Screenshot of OEHHA Chemical Data Explorer")
              ),
              br(),
              br()
            ),
            
            # Application 2
            box(
              title = "CalEnviroScreen",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              p("The CalEnviroScreen 4.0 Tool provides data and maps on pollution burden and population characteristics across California to identify disadvantaged communities."),
              tags$a(
                href = "https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-40", 
                target = "_blank", 
                tags$span(
                  icon("link"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff;",
                  " Visit the CalEnviroScreen 4.0 Tool"
                )
              ),
              # tags$span(
              #   tags$a(
              #     href = "https://github.com/oehha/calenviroscreen", 
              #     target = "_blank",
              #     icon("github"),
              #     style = "font-size: 2rem; font-weight: bold; color: #007bff; margin-left: 10px;",
              #     " View GitHub Repository"
              #   )
              # ),
              br(),
              tags$a(
                href = "https://oehha.ca.gov/calenviroscreen/report/calenviroscreen-40",
                target = "_blank",
                tags$img(src = "calenviroscreen_screenshot.png", width = "70%", alt = "Screenshot of CalEnviroScreen 4.0")
              ),
              br(),
              br()
            ),
            
            # Application 3
            box(
              title = "OEHHA Lead Leggett+ PBPK Modelling App",
              status = "info",
              solidHeader = TRUE,
              width = 12,
              p("The Leggett+ PBPK Modelling App allows for user-friendly pharmacokinetic modelling of lead for various exposure scenarios."),
              p("App under development!"),
              # tags$a(
              #   href = "https://oehha.ca.gov/lead-pbpk-app", 
              #   target = "_blank", 
              #   tags$span(
              #     icon("link"),
              #     style = "font-size: 2rem; font-weight: bold; color: #007bff;",
              #     " Visit the Lead PBPK Modelling App"
              #   )
              # ),
              tags$span(
                tags$a(
                  href = "https://github.com/ScottCoffin/Leggett-_lead_PBPK", 
                  target = "_blank",
                  icon("github"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff; margin-left: 10px;",
                  " View GitHub Repository"
                )
              ),
              br(),
              # tags$a(
              #   href = "https://oehha.ca.gov/lead-pbpk-app",
              #   target = "_blank",
              #   tags$img(src = "lead_pbpk_app_screenshot.png", width = "70%", alt = "Screenshot of OEHHA Lead PBPK Modelling App")
              # ),
              # br(),
              br()
            )
          )
        )
      )
      
      
    ) #
  )
)


################# SERVER ############

server <- function(input, output, session) {
  
##### RESET BUTTON #####
  observeEvent(input$reset, {
    # Reset experiment data
    experiment_data(initial_table_data)
    # Reset params data
    params_data(initial_params_data)
    # Reset reactive_params (clear or reinitialize as needed)
    reactive_params <- reactiveValues()  # Clears the reactiveValues object
    # Reset validation message
    validation_message(NULL)
    
    # Reset file inputs (using JavaScript)
    runjs("
    document.getElementById('upload_experiment_data').value = null;
    document.getElementById('upload_params_data').value = null;
  ")
    
    # Reset rendered handsontable
    output$experiment_table <- renderRHandsontable({
      rhandsontable(initial_table_data) %>% 
        hot_cols(strict = TRUE, allowInvalid = FALSE)
    })
    
    # Clear other outputs if necessary (example)
    output$params_table <- renderTable(NULL)
    
    # Optionally reset other UI elements
    # Reset button styles
    runjs("document.getElementById('run').style.backgroundColor = '#077336';")
  })
  
  # Download handler for the example CSV
  output$"example_experiment_data" <- downloadHandler(
    filename = function() {
      "example_experiment_data.csv"
    },
    content = function(file) {
      write.csv(initial_table_data, file, row.names = FALSE)
    }
  )
  
#################################################### Data entry and Upload #######################################
######## Experiment data table #####  
  # Reactive value to store experiment table data
  experiment_data <- reactiveVal(initial_table_data)
  
  # Upload and update experiment_data
  observeEvent(input$upload_experiment_data, {
 
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
      print("Uploaded data column names:")
      print(colnames(new_data))
      showNotification("Invalid data format. Missing required columns.", type = "error")
      return(NULL)
    }
    experiment_data(new_data)
    showNotification("Experiment data successfully updated.", type = "message")
    print("Experiment data successfully updated.")
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
    
    print("Updating experiment table based on new input...")
    
    # Reapply original data types and levels
    for (col in names(updated_data)) {
      # Check if the original type is not NA
      if (!is.na(original_types[col])) {
        # Use tryCatch to handle potential errors
        tryCatch({
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
        }, error = function(e) {
          # Handle the error gracefully
          print(paste("Error converting column", col, ":", e$message))
          # Optionally, you can set the column to NA or keep it unchanged
          updated_data[[col]] <- NA  # or leave it unchanged
        })
      } else {
        # Debugging output for missing type
        print(paste("Warning: Original type for column", col, "is NA. Skipping conversion."))
      }
    }
    
    # Update the reactive value with the new data
    experiment_data(updated_data)
  })
  
  ########################### Parameters Data Table ############
  # initalize input params
  # Reactive value to store params table data
  params_data <- reactiveVal(NULL)
  
  # Reactive value to store processed params table
  processed_params <- reactiveVal(NULL)
  
  # Observe changes in experiment_data() and update processed_params
  observeEvent(experiment_data(), {

    shiny::req(experiment_data())
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
      #calculate half-life from Clearance and Volume of Distribution if value is missing
      mutate(
        Half_Life_hr = case_when(
        is.na(Half_Life_hr) & !is.na(Clearance_L_per_kg_d) & !is.na(Volume_of_Distribution_L_per_kg) ~
          (log(2) * Volume_of_Distribution_L_per_kg) / Clearance_L_per_kg_d, #Ke = Vd / Cl
        T ~ Half_Life_hr
        ),
        Volume_of_Distribution_L_per_kg = case_when(
          is.na(Volume_of_Distribution_L_per_kg) & !is.na(Clearance_L_per_kg_d) & !is.na(Half_Life_hr) ~
            Clearance_L_per_kg_d / (log(2) / Half_Life_hr), #Vd  = Ke / Cl
          T ~ Volume_of_Distribution_L_per_kg
          )
      ) %>% 
      arrange(desc(PFAS), Species, Sex)  # Arrange for consistent display

    # Debugging: Inspect the processed data
    print("Processed params table:")
    print(params_processed)

    # Update reactive value
    processed_params(params_processed)
  })

  # Upload and update params_data (model parameters for species/sex combinations)
  observeEvent(input$upload_params_data, {

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

  #### Create data object for plotting
  TK_data <- reactive({
    
    full_tk <- readRDS("Additional files/Datasets/general/tk_df.rds") 
    tk_df <- full_tk %>% 
      select(chem, cas, species_name, sex, strain, tissue, route, standard_endpoint, standard_value, standard_unit,
             pubmed_id, authors, year, source, `Preferred value`) %>% 
      rename(PFAS = chem,
             CAS = cas,
             Species = species_name,
             Sex = sex,
             Strain = strain,
             Tissue = tissue,
             Route = route,
             Endpoint = standard_endpoint,
             Value = standard_value,
             Units = standard_unit,
             PubMedID = pubmed_id,
             Authors = authors,
             Year = year,
             Source = source
      ) %>% 
      filter(Endpoint %in% c("CL", "F", "HLe_invivo", "VDss", "kabs")) %>% 
      mutate(Endpoint = case_when(
        Endpoint == "CL" ~ "Clearance",
        Endpoint == "F" ~ "Bioavailable Fraction",
        Endpoint == "HLe_invivo" ~ "Elimination Half-Life",
        Endpoint == "VDss" ~ "Volume of Distribution",
        Endpoint == "kabs" ~ "Absorption Kinetic Constant"
      ))  %>% 
      mutate(across(c(PFAS, CAS, Species, Strain, Route, Endpoint, Units, `Preferred value`), as.factor)) %>% 
      droplevels() 
    
    tk_df
  })
  
  #### Full TK DataTable ####
  output$Full_TK_Datatable <- renderDT(server = F,{
    tk_df <- TK_data()
        print("Rendering TK datatable...")
        
        
        # Define the color palette based on your existing scheme
        color_palette <- c(
          "#fdb80b",  # Yellow
          "#0097ab",  # Teal
          "#0a4531",  # Dark Green
          "#017a3d",  # Medium Green
          "#f7a600",  # Gold (analogous to yellow)
          "#00758f",  # Deep teal
          "#044d34"   # Very dark green
        )
        # Define a corresponding text color palette (white for dark colors)
        text_palette <- c(
          "black",    # Yellow
          "white",    # Teal
          "white",    # Dark Green
          "white",    # Medium Green
          "black",    # Gold
          "white",    # Deep Teal
          "white"     # Very Dark Green
        )
        
    dt <- datatable(tk_df,
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
        # Target the "Species" column to determine colors
        columns = "Species",
        target = "row",
        # Use a custom JavaScript function to assign row-level colors
        backgroundColor = styleEqual(
          unique(tk_df$Species),
          color_palette          # Use the defined color palette
        ),
        color = styleEqual(
          unique(tk_df$Species),  # Map unique levels to text colors
          text_palette           # Use the defined text color palette
      )
      )
    
    dt
    
  })
  
 
  ############################### PBPK Params Table #############################
  # Reactive values for PBPK model parameters
  reactive_params <- reactiveValues()
  
  observeEvent(experiment_data(), {
  shiny::req(experiment_data())

  # Filter PBPK rows from experiment_data
  table_data <- experiment_data()
  pbpk_data <- subset(table_data, Model_Type == "PBPK")

  print("Processing PBPK data...")
  print(head(pbpk_data))  # Debugging: Inspect PBPK data

  for (sp in unique(pbpk_data$Species)) {
    for (sx in unique(pbpk_data$Sex)) {
      # Get species-specific model details
      model_info <- species_models[[sp]]

      # Skip processing if model_info is NULL or doesn't match Sex
      if (is.null(model_info) || model_info$Sex != sx) {
        print("Skipping processing for PBPK since model_info is NULL or doesn't match sex.")
        next
      }

      # Initialize parameter table if it doesn't exist
      if (is.null(reactive_params[[paste0(sp, "_", sx)]])) {
        reactive_params[[paste0(sp, "_", sx)]] <- data.frame(
          Parameter = if (!is.null(model_info$params)) names(model_info$params) else character(0),
          stringsAsFactors = FALSE
        )
        print("Initializing parameter table for PBPK model")
      }

      # Populate columns dynamically for the supported PFAS
      for (pfas in unique(pbpk_data$PFAS)) {
        if (pfas == model_info$PFAS) {
          reactive_params[[paste0(sp, "_", sx)]][[pfas]] <- signif(sapply(model_info$params, exp), 4)
          print("Reactive parameters:")
          print(reactive_params)
        } else {
          reactive_params[[paste0(sp, "_", sx)]][[pfas]] <- NA
        }
      }
    }
  }
})


  #just for debugging
  observeEvent(experiment_data(), {
    shiny::req(experiment_data())

    table_data <- experiment_data()
    pbpk_data <- subset(table_data, Model_Type == "PBPK")

    for (sp in unique(pbpk_data$Species)) {
      for (sx in unique(pbpk_data$Sex)) {
        model_info <- species_models[[sp]]

        if (is.null(model_info)) {
          message("Model info is NULL for Species =", sp)
          next
        }

        if (is.na(model_info$Sex)) {
          message("Sex is NA for Species =", sp, "Sex =", sx)
          next
        }

        if (model_info$Sex != sx) {
          message("Skipping unsupported combination: Species =", sp, "Sex =", sx)
          next
        }

        message("Processing Species =", sp, "Sex =", sx, "Supported PFAS =", model_info$PFAS)
      }
    }

    if (is.null(species_models) || length(species_models) == 0) {
      showNotification("Species models data is missing or invalid.", type = "error")
      return(NULL)
    }
    message("Species processing complete")
  })
  
  
  # Render UI for species/sex parameter tables in a grid layout
  output$species_sex_params_grid <- renderUI({
    message("rendering species_sex_params_grid")
    shiny::req(experiment_data())
    table_data <- experiment_data()
    pbpk_data <- subset(table_data, Model_Type == "PBPK")

    if (nrow(pbpk_data) == 0) return(NULL)

    species_sex_combinations <- unique(paste(pbpk_data$Species, pbpk_data$Sex, sep = "_"))

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

          if (!is.null(reactive_params[[paste0(species, "_", sex)]])) {
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
      }

      # Add each pair of species/sex tables to a new row
      ui_elements[[length(ui_elements) + 1]] <- fluidRow(row_elements)
    }

    ui_elements
  })


  # Render and update species/sex-specific parameter tables
  observe({
    for (sp in unique(experiment_data()$Species)) {
      for (sx in unique(experiment_data()$Sex)) {
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
    for (sp in unique(experiment_data()$Species)) {
      for (sx in unique(experiment_data()$Sex)) {
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
          HTML("All parameters are sufficiently available to run the simulation! <br> Please hit `Run Simulation` button to proceed."),
          style = "color: green; font-size: 20px; text-align: center;"
        )
      )
    } else {
      validation_message(
        tags$p(
          HTML(validation_result),  # Use HTML to interpret line breaks
          style = "color: red; font-size: 20px; text-align: center;"
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
          HTML("All parameters are sufficiently available to run the simulation! <br> Please hit `Run Simulation` button below to proceed."),
          style = "color: green; font-size: 20px; text-align: center;"
        )
      } else {
        tags$p(
          HTML(validation_result),  # Use HTML to interpret line breaks
          style = "color: red; font-size: 20px; text-align: center;"
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
        Sex = as.character(Sex),
        Model_Type = as.character(Model_Type)
      )
    
    print(head(params_data()))
    
    params <- params_data() %>%
      mutate(
        Species = as.character(Species),
        PFAS = as.character(PFAS),
        Sex = as.character(Sex),
        Model_Type = as.character(Model_Type)
      )
    
    # Join experiment_data with params_data
    experiment_with_params <- exp_data %>%
      left_join(params, by = c("Species", "Sex", "PFAS", "Model_Type"))
    
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
        tgrid <- tgrid(0, exposure_duration_hrs + 96, by = 1) #model every 1 hour
      
        
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
        times <- seq(0, exposure_duration_hrs + 96, by = 1) # model every 1 hr
        
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
    
    # Create a new column for custom legend labels
    sim_results <- sim_results %>%
      mutate(Legend_Label = paste(Species, PFAS, Sex, Model_Type, paste0(Dose_mg_per_kg, " mg/kg"), sep = " | "))
    
    # Determine the number of unique legend labels
    num_colors_needed <- length(unique(sim_results$Legend_Label))
    custom_colors <- hcl(seq(15, 375, length.out = num_colors_needed), 100, 65)  # Adjust parameters as needed
    
    # Plot with the custom legend labels
    p <- ggplot(sim_results, aes(x = Time, y = Concentration, color = Legend_Label,
                                 group = interaction(Species, PFAS, Sex, Model_Type, Exposure_Duration_Days, Interval_Hours, Dose_mg_per_kg), # Explicit grouping
                                 text = paste("Chemical:", PFAS, "<br>",
                                              "Species:", Species, "<br>",
                                              "Sex:", Sex, "<br>",
                                              "Model:", Model_Type, "<br>",
                                              "Dose:", Dose_mg_per_kg, "mg/kg-d", "<br>",
                                              "Exposure Duration:", Exposure_Duration_Days, "days", "<br>",
                                              "Dosing Interval:", Interval_Hours, "hrs", "<br>",
                                              "Time:", signif(Time, 4), "Days", "<br>",
                                              "Modeled Serum Concentration:", signif(Concentration,3), "mg/L", "<br>")
                                 )) +
      geom_line(size = 1) +
      labs(title = "Plasma Concentration Over Time",
           x = "Time (days)", 
           y = "Concentration (ug/ml)",
           color = "") +
      #scale_color_discrete_c4a_cat("hcl.dark3") +
      scale_color_manual(values = custom_colors) +
      xlim(0, max(sim_results$Time, na.rm = TRUE) + 4) +
      theme_minimal(base_size = 13) + 
      theme(legend.title = element_blank())
    
    # Convert to plotly
    ggplotly(p, tooltip = "text")
  })
  
  ############################ Scatterplot of measured vs. modeled #############
  # Render scatter plot for modeled vs. measured concentrations
  output$scatterPlot <- renderPlotly({
    sim_results <- simulation_results()
    exp_data <- experiment_data() %>% filter(!is.na(Time_Serum_Collected_hr))
    
    # Pre-calculate time in days for joining
    exp_data <- exp_data %>%
      mutate(Time_in_days = Time_Serum_Collected_hr / 24)
    
    # Join with sim_results
    measured_conc <- exp_data %>%
      left_join(sim_results, by = c("Time_in_days" = "Time", "Species" = "Species", "PFAS", "Sex", "Model_Type", "Dose_mg_per_kg", "Exposure_Duration_Days", "Interval_Hours")) %>%
      # rename("Predicted Concentration (mg/L)" = Concentration,
      #        "Model" = Model_Type,
      #        "Dose (mg/kg)" = Dose_mg_per_kg,
      #        "Dosing Interval" = Interval_Hours,
      #        "Exposure (Days)" = Exposure_Duration_Days,
      #        "Hr Serum Collected" = Time_Serum_Collected_hr,
      #        "Measured Concentration (mg/L)" = Serum_Concentration_mg_L
      # ) %>% 
   #   select(-Time_in_days) %>% 
    #  mutate(across(c(PFAS, Species, Sex, Model), as.factor)) %>% 
      rename(modeled_concentration = Concentration) %>% 
      mutate(Legend_Label = paste(Species, PFAS, Model_Type, paste0(Dose_mg_per_kg, " mg/kg"), sep = " | "))
    
    #legacy approach
    # measured_conc <- exp_data %>%
    #   mutate(modeled_concentration = map2_dbl(Time_Serum_Collected_hr, Species, function(time, species) {
    #     subset(sim_results, Time == time / 24 & Species == species)$Concentration[1]
    #   })) %>%
    #   filter(!is.na(modeled_concentration)) %>% 
    #   # custom legend for plot
    #   mutate(Legend_Label = paste(Species, PFAS, Model_Type, paste0(Dose_mg_per_kg, " mg/kg"), sep = " | "))
    
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
    
    # Determine the number of unique legend labels
    num_colors_needed <- length(unique(measured_conc$Legend_Label))
    custom_colors <- hcl(seq(15, 375, length.out = num_colors_needed), 100, 65)  # Adjust parameters as needed
    
    #make scatterplot
    p <- ggplot(measured_conc, aes(y = modeled_concentration, 
                                   x = Serum_Concentration_mg_L,
                                   color = Legend_Label,
                                   shape = Sex,
                                   group = interaction(Species, PFAS, Sex, Model_Type, Exposure_Duration_Days, Interval_Hours, Dose_mg_per_kg), # Explicit grouping
                                   text = paste("Chemical:", PFAS, "<br>",
                                                "Species:", Species, "<br>",
                                                "Sex:", Sex, "<br>",
                                                "Dose:", Dose_mg_per_kg, "mg/kg-d", "<br>",
                                                "Exposure Duration:", Exposure_Duration_Days, "days", "<br>",
                                                "Dosing Interval:", Interval_Hours, "hrs", "<br>",
                                                "Measured Serum Concentration:", Serum_Concentration_mg_L, "mg/L", "<br>",
                                                "Model:", Model_Type, "<br>",
                                                "Modeled Serum Concentration:", signif(modeled_concentration,3), "mg/L", "<br>"))) +
      geom_point(size = 2, alpha = 0.9) +
      scale_x_log10(limits = expanded_limits) +
      scale_y_log10(limits = expanded_limits) +
      #scale_color_discrete_c4a_cat("hcl.dark3", name = "") +
      scale_color_manual(values = custom_colors) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") + # Add 10x and 0.1x deviation lines
      geom_abline(slope = 1, intercept = log10(10), linetype = "dotted", color = "gray70") +
      geom_abline(slope = 1, intercept = log10(0.1), linetype = "dotted", color = "gray70") +
      labs(title = "Modeled vs. Measured Concentration",
           color = "",
           shape = "",
           x = "Modeled Concentration (mg/L)",
           y = "Measured Concentration (mg/L)") +
      theme_minimal(base_size = 15) +
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
    
    # Pre-calculate time in days for joining
    exp_data <- exp_data %>%
      mutate(Time_in_days = Time_Serum_Collected_hr / 24)
    
    # Join with sim_results
    data <- exp_data %>%
      left_join(sim_results, by = c("Time_in_days" = "Time", "Species" = "Species", "PFAS", "Sex", "Model_Type", "Dose_mg_per_kg", "Exposure_Duration_Days", "Interval_Hours")) %>%
      rename("Predicted Concentration (mg/L)" = Concentration,
             "Model" = Model_Type,
             "Dose (mg/kg)" = Dose_mg_per_kg,
             "Dosing Interval" = Interval_Hours,
             "Exposure (Days)" = Exposure_Duration_Days,
             "Hr Serum Collected" = Time_Serum_Collected_hr,
             "Measured Concentration (mg/L)" = Serum_Concentration_mg_L
             ) %>% 
      select(-Time_in_days) %>% 
      mutate(across(c(PFAS, Species, Sex, Model), as.factor))
    
    #build datatable
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
        columns = "Model"  # Specifies to base coloring on the data_available column
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
      group_by(Species, PFAS, Sex, Model_Type, Dose_mg_per_kg, Interval_Hours, Exposure_Duration_Days) %>%
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
    summary <- summary_stats() %>% 
      rename("Model" = Model_Type,
             "Dose (mg/kg)" = Dose_mg_per_kg,
             "Exposure (Days)" = Exposure_Duration_Days,
             "Maximum Serum Concentration (mg/L)" = C_max,
             "Area Under Curve Serum Concentration (mg/L)" = AUC,
             "Time-Weighted Average Serum Concentraiton (mg/L)" = C_TWA)

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
        columns = "Model"  # Specifies to base coloring on the data_available column
      ) %>%
      formatSignif(columns = c("Maximum Serum Concentration (mg/L)", "Area Under Curve Serum Concentration (mg/L)", "Time-Weighted Average Serum Concentraiton (mg/L)"),
                   digits = 4)
  })

#### Summary stats heatmap ####
  output$heatmap <- renderPlotly({
    
    params_heat <- summary_stats() %>%
      pivot_longer(
        cols = c("C_max", "C_TWA"),
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      mutate(Metric = case_when(
        Metric == "C_max" ~ "Maximum Serum Concentration (mg/L)",
        Metric == "C_TWA" ~ "Time-Weighted Average Serum Concentration (mg/L)")) %>% 
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
        group = interaction(Species, PFAS, Sex, Model_Type, Exposure_Duration_Days, Interval_Hours, Dose_mg_per_kg), # Explicit grouping
        fill = normalized_value + 0.001,
        text = paste0(
          "Chemical: ", PFAS, "<br>",
          "Species: ", Species, "<br>",
          "Sex: ", Sex, "<br>",
          "Exposure Duration: ", Exposure_Duration_Days, " days", "<br>",
          "Dose: ", Dose_mg_per_kg, " mg/kg", "<br>",
          "Model Type: ", Model_Type, "<br>",
          "Computed Metric: ", Metric, "<br>",
          "Value: ", signif(Value,2), " mg/L"
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
        axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none"
      ) +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 15))  # Wrap x-axis text
    
    
    ggplotly(params_heat, tooltip = "text")
  })
  
  
} # close server

shinyApp(ui = ui, server = server)

