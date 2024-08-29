library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(ggplot2)
library(mrgsolve)
library(magrittr)
library(dplyr)

# load models
RatPBPK.code <- readRDS("../Additional files/Results/Workplace/ratPBPK.RDS")
ratpbpk <- mcode("Ratpbpk", RatPBPK.code)

MicePBPK.code <- readRDS("../Additional files/Results/Workplace/micePBPK.RDS")
micepbpk <- mcode("Micepbpk", MicePBPK.code)

humanPBPK.code <- readRDS("../Additional files/Results/Workplace/humanPBPK.RDS")
humanpbpk <- mcode("Humanpbpk", humanPBPK.code)

monkeyPBPK.code <- readRDS("../Additional files/Results/Workplace/monkeyPBPK.RDS")
monkeypbpk <- mcode("Monkeypbpk", monkeyPBPK.code)

# load best fit parameters
rat_params <- readRDS("../Additional files/Results/Workplace/rat.MCMC.rds")
rat_best <- rat_params$bestpar
mouse_params <- readRDS("../Additional files/Results/Workplace/mouse.MCMC.rds")
mouse_best <- mouse_params$bestpar
human_params <- readRDS("../Additional files/Results/Workplace/human.MCMC.rds")
human_best <- human_params$bestpar
monkey_params <- readRDS("../Additional files/Results/Workplace/monkey.MCMC.rds")
monkey_best <- monkey_params$bestpar


######################### UI #######################
# Define UI
# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "PBPK Model Simulation"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Experiment Inputs", tabName = "inputs", icon = icon("flask")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Model Parameters", tabName = "parameters", icon = icon("sliders-h"))
    ),
    
    collapsed = TRUE  # Start with the sidebar collapsed
    
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "inputs",
              fluidRow(
                box(title = "Input Parameters", status = "primary", solidHeader = TRUE, width = 12,
                    selectInput("species", "Select Species", 
                                choices = c("Rat", "Mouse", "Human", "Monkey"), 
                                selected = "Rat", multiple = TRUE),
                    numericInput("dose_mg_per_kg", "Dose (mg/kg)", value = 50),
                    numericInput("interval_hours", "Interval between doses (hours)", value = 12),
                    numericInput("exposure_duration_days", "Exposure Duration (days)", value = 30),
                    numericInput("dose_fraction", "Dose Fraction (e.g., 0.5 for half-dose)", value = 0.5),
                    numericInput("dose_frequency_per_day", "Dose Frequency per Day", value = 2),
                    
                    actionButton("run", "Run Simulation")
                )
              )
      ),
      
      tabItem(tabName = "results",
              fluidRow(
                box(title = "Concentration Plot", status = "primary", solidHeader = TRUE, width = 12,
                    plotlyOutput("concentrationPlot")),
                
                box(title = "Summary Statistics", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("summary_table"))
              ),
              
              fluidRow(
                box(title = "Download Data", status = "primary", solidHeader = TRUE, width = 6,
                    downloadButton("downloadData", "Download Simulation Data")),
                
                box(title = "Download Plots", status = "primary", solidHeader = TRUE, width = 6,
                    downloadButton("downloadPlots", "Download Plots"))
              )
      ),
      
      tabItem(tabName = "parameters",
              fluidRow(
                box(title = "Editable Model Parameters", status = "primary", solidHeader = TRUE, width = 12,
                    DTOutput("model_parameters"))
              )
      )
    )
  )
)


server <- function(input, output, session) {
  
  species_models <- list(
    Rat = list(model = ratpbpk, params = rat_best),
    Mouse = list(model = micepbpk, params = mouse_best),
    Human = list(model = humanpbpk, params = human_best),
    Monkey = list(model = monkeypbpk, params = monkey_best)
  )
  
  # Define data_input function
  data_input <- reactive({
    if (!is.null(input$file)) {
      # Batch mode: read CSV
      read.csv(input$file$datapath)
    } else {
      # Single experiment mode: use inputs from the UI
      data.frame(
        dose_mg_per_kg = input$dose_mg_per_kg,
        interval_hours = input$interval_hours,
        exposure_duration_days = input$exposure_duration_days,
        dose_fraction = input$dose_fraction,
        dose_frequency_per_day = input$dose_frequency_per_day
      )
    }
  })
  
  reactive_params <- reactiveVal({
    data.frame(
      Parameter = names(rat_best),  # Assuming 'rat_best' is the base parameters used
      Value = signif(sapply(rat_best, exp), 4)
    )
  })
  
  output$model_parameters <- renderDT({
    datatable(reactive_params(), editable = TRUE, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  
  
  run_simulation <- function(species, dose_mg_per_kg, interval_hours, exposure_duration_days, dose_fraction, dose_frequency_per_day) {
    species_info <- species_models[[species]]
    model <- species_info$model
    pars <- species_info$params
    
    BW <- 0.3
    dose_mg <- dose_mg_per_kg * BW * dose_fraction
    interval_between_doses <- 24 / dose_frequency_per_day
    total_doses <- exposure_duration_days * dose_frequency_per_day
    
    ex <- ev(ID = 1, amt = dose_mg, ii = interval_between_doses, addl = total_doses - 1, cmt = "AST", replicate = FALSE)
    tgrid <- tgrid(0, interval_hours * (exposure_duration_days - 1) + 24 * 100, 1)
    
    output <- model %>%
      param(pars) %>%
      Req(Plasma) %>%
      update(atol = 1E-8, maxsteps = 10000) %>%
      mrgsim_d(data = ex, tgrid = tgrid)
    
    data.frame(Species = species, Time = output$time / 24, Concentration = signif(output$Plasma, 4))
  }
  
  simulation_results <- eventReactive(input$run, {
    req(input$species)
    data <- data_input()
    species_results <- lapply(input$species, function(species) {
      run_simulation(
        species,
        data$dose_mg_per_kg[1],
        data$interval_hours[1],
        data$exposure_duration_days[1],
        data$dose_fraction[1],
        data$dose_frequency_per_day[1]
      )
    })
    bind_rows(species_results)
  })
  
  output$concentrationPlot <- renderPlotly({
    req(simulation_results())
    
    plot_data <- simulation_results()
    
    p <- ggplot(plot_data, aes(x = Time, y = Concentration, color = Species)) +
      geom_line(size = 1) +
      labs(title = "Plasma Concentration Over Time",
           x = "Time (days)", 
           y = "Concentration (ug/ml)") +
      theme_minimal()
    
    ggplotly(p, tooltip = c("x", "y", "color"))
  })
  
  output$summary_table <- renderDT({
    summary_stats <- simulation_results() %>%
      group_by(Species) %>%
      summarize(Cmax = max(Concentration), 
                AUC = sum(Concentration) * diff(Time[1:2]), 
                TWAC = AUC / max(Time))
    
    datatable(summary_stats)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulation_data", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(simulation_results(), file)
    }
  )
  
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste("simulation_plots", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(ggplotly(ggplot()), file)
    }
  )
}



# Run the app
shinyApp(ui = ui, server = server)
