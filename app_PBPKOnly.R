library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(mrgsolve)
library(dplyr)
library(purrr)

#### PBPK Models ####
# Load PBPK models and parameters for different species
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

# Define species models and parameters
species_models <- list(
  Rat = list(model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
  Mouse = list(model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
  Human = list(model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
  Monkey = list(model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
)


#### Simple PK Models ####
# Default parameters for the simpler PK models
default_params <- list(
  Dose_mg_kg_d = 100,   # Total daily dose (mg/kg)
  HLe_invivo = 6,       # Half-life in hours
  tau_hr = 12,          # Dosing interval in hours
  Vd = 50,              # Volume of distribution (L/kg)
  kabs = 0.5,           # Absorption rate constant (1/h)
  exposure_duration_hr = 72,  # Exposure duration in hours
  n_extra_time = 48     # Extra time to simulate after exposure (hours)
)

# Simplified PK models
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

############################# UI ###################################
ui <- dashboardPage(
  dashboardHeader(title = "PBPK Model Simulation (PFOS, Males Only, Oral)"),
  
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
                box(title = "PFOS PBPK Model Info", status = "primary", solidHeader = TRUE, width = 12,
                    p("This RShiny app allows users to run the optimized models from", a("Chou & Lin et al. (2019). Bayesian evaluation of a physiologically based pharmacokinetic (PBPK) model for perfluorooctane sulfonate (PFOS) to characterize the interspecies uncertainty between mice, rats, monkeys, and humans: Development and performance verification",
                                                                                         href = "https://www.sciencedirect.com/science/article/pii/S016041201930203X")),
                    p("The model was only optimized for males and an oral route of exposure. The model consists of four organ compartments (plasma, liver, kidney, and rest of body). For more information, please refer to their publication."),
                    ),
                box(title = "Input Parameters", status = "primary", solidHeader = TRUE, width = 12,
                    selectInput("species", "Select Species", 
                                choices = c("Rat", "Mouse", "Human", "Monkey"), 
                                selected = "Rat", multiple = TRUE),
                    uiOutput("body_weight_inputs"),  # Dynamically generate body weight inputs
                    numericInput("dose_mg_per_kg", "Dose (mg/kg)", value = 50),
                    numericInput("interval_hours", "Interval between doses (hours)", value = 12),
                    numericInput("exposure_duration_days", "Exposure Duration (days)", value = 30),
                    numericInput("dose_fraction", "Dose Fraction (e.g., 0.5 for half-dose)", value = 0.5),
                    numericInput("dose_frequency_per_day", "Dose Frequency per Day", value = 2),
                    numericInput("auc_start", "AUC Calculation Start Time (days)", value = 0),
                    numericInput("auc_end", "AUC Calculation End Time (days)", value = 30),
                    actionButton("run", "Run Simulation")
                )
              )
      ),
      
      tabItem(tabName = "results",
              fluidRow(
                box(title = "Concentration Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("concentrationPlot")),
                br(),
                h3("Estimate serum concentration at a specific time point:"),
                numericInput("time_point", "Time Point (hours)", value = NULL),
                tableOutput("concentration_at_time")
                ),
              fluidRow(
                box(title = "Summary Statistics (mg/L)", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("summary_table"))
              ),
      ),
      
      tabItem(tabName = "parameters",
              fluidRow(
                uiOutput("species_params")  # Dynamically generate UI for each species' parameters
              )
      )
    )
  )
)

#################################### SERVER #######################
server <- function(input, output, session) {
  
  # Initialize reactive values for model parameters and body weights for each species
  reactive_params <- reactiveValues()
  reactive_bw <- reactiveValues(
    Rat = 0.3,
    Mouse = 0.025,
    Human = 82.3,
    Monkey = 3.5
  )
  
  observe({
    for (sp in input$species) {
      if (is.null(reactive_params[[sp]])) {
        species <- species_models[[sp]]
        reactive_params[[sp]] <- data.frame(
          Parameter = names(species$params),
          Value = signif(sapply(species$params, exp), 4)
        )
        reactive_bw[[sp]] <- species$default_bw
      }
    }
  })
  
  # Corrected renderUI for body weight inputs with debugging
  output$body_weight_inputs <- renderUI({
    print("Rendering body weight inputs...")
    print(input$species)  # Debugging output
    
    # Manually check if input$species is valid
    if (is.null(input$species) || length(input$species) == 0) {
      return(NULL)
    }
    
    # If input$species is valid, create the UI elements
    lapply(input$species, function(sp) {
      print(paste("Creating input for species:", sp))  # Debugging output
      numericInput(
        inputId = paste0("bw_", sp), 
        label = paste(sp, "Body Weight (kg)"), 
        value = reactive_bw[[sp]]  # Ensure `reactive_bw` is available as a reactive list
      )
    })
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
  
  # Function to run the PBPK model simulation for a single species
  run_simulation <- function(species, bw, dose_mg_per_kg, interval_hours, exposure_duration_days, dose_fraction, dose_frequency_per_day) {
    pars <- setNames(as.list(reactive_params[[species]]$Value), reactive_params[[species]]$Parameter)
    model <- species_models[[species]]$model
    
    dose_mg <- dose_mg_per_kg * bw * dose_fraction
    interval_between_doses <- 24 / dose_frequency_per_day
    total_doses <- exposure_duration_days * dose_frequency_per_day
    
    ex <- ev(ID = 1, amt = dose_mg, ii = interval_between_doses, addl = total_doses - 1, cmt = "AST", replicate = FALSE)
    tgrid <- tgrid(0, interval_hours * (exposure_duration_days - 1) + 24 * 100, 1)
    
    print(head(pars))
    
    
    output <- model %>%
      param(pars) %>%
      Req(Plasma) %>%
      update(atol = 1E-8, maxsteps = 10000) %>%
      mrgsim_d(data = ex, tgrid = tgrid)
    
    cat("Debugging inputs:\n")
    cat("Dose:", dose_mg, "mg/kg\n")
    cat("Body Weight:", bw, "kg\n")
    cat("Interval:", interval_hours, "hours\n")
    cat("Exposure Duration:", exposure_duration_days, "hours\n")
    # cat("ex:", ex, "\n")
    # cat("tgrid:", tgrid, "\n")
    
    data.frame(Species = species, Time = output$time / 24, Concentration = signif(output$Plasma, 4))
  }
  
  # Run the simulation for all selected species when the button is clicked
  simulation_results <- eventReactive(input$run, {
    do.call(rbind, lapply(input$species, function(sp) {
      bw <- input[[paste0("bw_", sp)]]
      run_simulation(
        sp,
        bw,
        input$dose_mg_per_kg,
        input$interval_hours,
        input$exposure_duration_days,
        input$dose_fraction,
        input$dose_frequency_per_day
      )
    }))
  })
  
  #update time point for observed serum level to be 24 hr following last dose
  observe({
    # Update the time_point input when exposure_duration_days is changed
    updateNumericInput(session, "time_point",
                       value = (input$exposure_duration_days * 24) + 24)
  })
  
  
  # Extract concentration at specific time point
  concentration_at_time <- reactive({
    sim_results <- simulation_results()
    if (is.null(sim_results) || nrow(sim_results) == 0) {
      return(NULL)
    }
    
    # Filter the simulation result at the specified time point
    sim_results %>%
      filter(Time == input$time_point / 24) %>%
      select(Species, Concentration)
  })

### Concentration-time-point table ####
  # Display concentration at the specified time point
  output$concentration_at_time <- renderTable({
    concentration_at_time()
  })
  
  # Calculate summary statistics (C_max, C_TWA, AUC) for a specified time period
  summary_stats <- reactive({
    print("Calculating summary statistics...")
    
    # Manually check if the simulation results exist and are valid
    sim_results <- simulation_results()
    if (is.null(sim_results) || nrow(sim_results) == 0) {
      return(NULL)
    }
    
    sim_results %>%
      group_by(Species) %>%
      summarize(
        C_max = max(Concentration),
        AUC = {
          time_subset <- Time >= input$auc_start & Time <= input$auc_end
          calculate_auc(Time[time_subset], Concentration[time_subset])
        },
        C_TWA = AUC / (input$auc_end - input$auc_start)
      ) %>%
      ungroup()
  })
  
  
  # Function to calculate AUC using the trapezoidal rule
  calculate_auc <- function(time, concentration) {
    stopifnot(length(time) == length(concentration))
    time <- sort(time)
    
    auc <- sum((concentration[-1] + concentration[-length(concentration)]) / 2 * diff(time))
    return(auc)
  }

###### Concentration-Time Plot #####
  # Render the concentration plot
  output$concentrationPlot <- renderPlotly({
    print("Rendering concentration plot...")
    
    # Manually check if the simulation results exist and are valid
    sim_results <- simulation_results()
    if (is.null(sim_results) || nrow(sim_results) == 0) {
      return(NULL)
    }
    
    plot_data <- sim_results
    
    p <- ggplot(plot_data, aes(x = Time, y = Concentration, color = Species)) +
      geom_line(size = 1) +
      labs(title = "Plasma Concentration Over Time",
           x = "Time (days)", 
           y = "Concentration (ug/ml)") +
      theme_minimal()
    
    ggplotly(p, tooltip = c("x", "y", "color"))
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
              extensions = 'Buttons',
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv'),
                columnDefs = list(list(targets = "_all", className = "dt-center"))
              )
    ) %>%
      formatSignif(columns = c("C_max", "AUC", "C_TWA"), digits = 4)
  })
  
} # close the app





# Run the app
shinyApp(ui = ui, server = server)
