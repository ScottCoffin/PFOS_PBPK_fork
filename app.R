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
library(shinycssloaders) #shows spinner when loading


#### run this before deploying to update date ###
# Record the date of app deployment.
record_deployment_date <-  function(deployment_history_file = "deployment_history.txt") {
  if (!file.exists(deployment_history_file)) {file.create(deployment_history_file)}
  # record the time
  deployment_time <- Sys.Date()
  cat(paste0(deployment_time, "\n"), file = deployment_history_file, append = TRUE)}
# be sure to run this to update text file before deployment. 
#record_deployment_date() #keep commented out for deployment, but execute when updating

# Return the last recorded deployment date of the application.
load_deployment_date <-
  function(deployment_history_file = "deployment_history.txt") {
    deployment_history <- readLines(deployment_history_file)
    
    # return the most recent line
    deployment_history[[length(deployment_history)]]
  }

############### STATIC DATA ############

############## PFOS PBPK Models and Parameters ######
############ - code and models from Chou & Lin 2019 (https://www.sciencedirect.com/science/article/pii/S016041201930203X) --
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

################# PFOA Male Rat PBPK Model ##########
##### Code and data from Fan et al. (2025) ####
#load PBPK models for different PFAS


################# PFOS, PFOA, PFHxS Male Mouse Models ######
######--- code and data from Zhu et al. (2023) https://ehp.niehs.nih.gov/doi/10.1289/EHP11969
### their code is incomplete!!

##### Fabian Fischer's 2025 Active Transporter PBPK Model ####
# load Fischer model function
source("Additional files/Code/Fischer 2025/Fischer_PBPK.R")
# loads "Fischer_PBPK_mouse" function
## ---- Fischer (MassTransferPBPK) parameter catalog (used for validation & UI) ----
# Path already used by your server call:
FISCHER_PARAM_PATH <- "Additional files/Datasets/Fischer/pfas_parameters.csv"

# Read once (be tolerant if file is missing at dev-time)
fischer_params <- tryCatch({
  if (file.exists(FISCHER_PARAM_PATH)) {
    readr::read_csv(FISCHER_PARAM_PATH, show_col_types = FALSE)
  } else {
    NULL
  }
}, error = function(e) NULL)

# What PFAS are supported by the MassTransferPBPK?
fischer_supported_pfas <- if (!is.null(fischer_params) && "PFAS" %in% names(fischer_params)) {
  sort(unique(fischer_params$PFAS))
} else character(0)

# (Current Fischer function is mouse/male only)
fischer_supported_species <- "Mouse"
fischer_supported_sex     <- "Male"

# helper to validate a single MassTransferPBPK combo
# tweak your check to use the reactive list
# MassTransferPBPK: Mouse/Male only + PFAS present in Fischer param file
is_mass_transfer_supported <- function(pfas, species, sex,
                                       pfas_pool       = fischer_supported_pfas,
                                       allowed_species = fischer_supported_species,
                                       allowed_sex     = fischer_supported_sex) {
  n <- length(pfas)
  if (length(pfas_pool) == 0) return(rep(FALSE, n))
  (pfas %in% pfas_pool) & (species == allowed_species) & (sex == allowed_sex)
}


## ---- Palettes used across plots ----
# Model-type colors (stable, named)
model_type_colors <- c(
  "PBPK"               = "#3A9AB2",
  "single-compartment" = "#6FB2C1",
  "two-compartment"    = "#91BAB6",
  "biphasic"           = "#91BAD6",
  "MassTransferPBPK"   = "#1E88E5"  # Fischer
)

## ---- Model preference (used everywhere) ----
# Assumption: biphasic ≈ two-compartment (place it between 2- and 1-comp)
model_preference <- c("MassTransferPBPK","PBPK","two-compartment","biphasic","single-compartment")
pref_rank <- function(x) match(x, model_preference, nomatch = NA_integer_)


# Compartment linetypes (stable, named)
compartment_linetypes <- c(
  "Plasma"  = "solid",
  "Blood"   = "solid",
  "Liver"   = "longdash",
  "Kidneys" = "dotdash",
  "Gut"     = "twodash",
  "Rest"    = "dotted"
)

# Compartment colors if/when you want color by compartment
compartment_colors <- c(
  "Plasma"  = "#1f38f7",
  "Blood"   = "#1f78b4",
  "Liver"   = "#33a02c",
  "Kidneys" = "#6a3d9a",
  "Gut"     = "#ff7f00",
  "Rest"    = "#b15928"
)

##-------------------------

# general TK params for non PBPK-models
# This data is prepared in the TK_data_prep.R script. The values in this are imported from the PFHpA Repo (https://github.com/ScottCoffin/PFHpA.git) in output/data/dose_to_serum_data.xlsx
tk_params <- readRDS("Additional files/Datasets/general/tk_params.rds")

#imprt entire EAS-E Suite + additional data to get sources
tk_df <- readRDS("Additional files/Datasets/general/tk_df.rds") %>% 
  mutate(source = paste(authors, year)) %>% 
  select(source, pubmed_id, doi) %>% 
  distinct(source, .keep_all = T) 
  

## load 'mismatch data' for serving pup interactive plotly ##
# this is generated in the dose_to_serum.Rmd script in the PFHpA repo!
mismatch_data <- readRDS("Additional files/Datasets/general/mismatch_data_for_shiny.rds") %>% 
  mutate(Sex = case_when(
    sex == "F" ~ "Female",
    sex == "M" ~ "Male"
  )) %>% 
  mutate(PFAS = case_when(
    chem == "PFHPA" ~ "PFHpA",
    chem == "PFHXA" ~ "PFHxA",
    chem == "PFHXS" ~ "PFHxS",
    chem == "PFPRA" ~ "PFPrA",
    chem == "TFA" ~ "TRIFLUOROACETIC ACID",
    T ~ chem
  )) %>% 
  # join with tk_params to get doi based on source
  left_join(tk_df, by = "source")

# Define species models and parameters
species_models <- list(
  Rat = list(PFAS = "PFOS", Sex = "Male", model = ratpbpk, params = rat_best, default_bw = 0.3),   # 0.3 kg for rat
  Mouse = list(PFAS = "PFOS", Sex = "Male",model = micepbpk, params = mouse_best, default_bw = 0.025), # 0.025 kg for mouse
  Human = list(PFAS = "PFOS", Sex = "Male",model = humanpbpk, params = human_best, default_bw = 82.3),  # 82.3 kg for human
  Monkey = list(PFAS = "PFOS", Sex = "Male",model = monkeypbpk, params = monkey_best, default_bw = 3.5)  # 3.5 kg for monkey
)
 
# PFAS levels available  (ALPHABETICAL)
pfas_levels <- sort(unique(c(as.character(tk_params$PFAS), fischer_supported_pfas)))


# Initialize editable table data with unrestricted input fields
initial_table_data <- data.frame(
  PFAS = factor(c("PFBA", "PFBA", "PFHxA", "PFOA", "PFOA", "PFOS", "GENX", "PFBS", "PFOA"),
                levels = pfas_levels), # ensure Fischer model in there too!
  Species = factor(c("Rat", "Mouse", "Mouse", "Rat", "Rat", "Rat", "Rat", "Mouse", "Mouse"),
                   levels = c("Rat", "Mouse", "Human", "Monkey")),
  Sex = factor(c("Female", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male")),
  # Route = factor(c("Intravenous", "Intragastric", "Oral", "Oral") ,
  #                levels = c("Intravenous", "Intragastric", "Intraperitoneal", "Oral", "Dermal")),
  Model_Type = factor(c("single-compartment", "two-compartment", "single-compartment", "two-compartment", "biphasic", "PBPK", "biphasic", "MassTransferPBPK", "MassTransferPBPK"),
                      levels = model_preference),
  #Body_Weight_kg = c(0.3, 0.025, 0.3, 0.3),
  Dose_mg_per_kg = c(175, 35, 62.5, 50, 30, 2.5, 30, 5, 5),
  Interval_Hours = as.integer(c(24, 24, 24, 24, 24, 24, 24, 24, 24)),
  Exposure_Duration_Days = as.integer(c(17, 28, 28, 28, 28, 28, 28, 30, 30)),
  Time_Serum_Collected_hr = as.integer(c(432, 696, 696, 696, 696, 696, 696, 720, 720)),
  Serum_Concentration_mg_L = as.numeric(c(4.44, 86, 3.1, 9.326, 51.65, 173.7, NA, 0.0895, 0.370)),
  stringsAsFactors = FALSE
)

# validation data (created in PFHpA repo)
validation_data <- read_excel("Additional files/Datasets/general/validation_data.xlsx")

# Capture the original column types and levels
original_types <- sapply(initial_table_data, class)
original_levels <- lapply(initial_table_data, levels)


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

##The biphasic model for PFAAs that are known to follow such kinetic profile (e.g., GenX, PFOA, PFBS, PFOS) is as follows (doi.org/10.1016/j.envint.2018.01.011):
two_phase <- function(D_per_dose, 
                               VDc, # VDc: beta phase volume of distribution (steady-state)
                               #VD_ss in Gomis et al. (2018), but VDc in EASESUITE - beta phase Vd
                               k_beta, # k_beta: elimination rate constant (beta phase)
                               #ke in Gomis et al. (2018), but k_beta in EASESUITE
                               kabs, # kabs: absorption rate constant
                               VD2, # VD2: alpha phase volume of distribution (initial distribution)
                               #VD0 in gomis et al. (2018), but VD2 in EASESuite
                               k_alpha, # k_alpha: distribution rate constant (alpha phase)
                               #"k_d" in Gomis et al. (2018, but k_alpha in EASESUITE)
                               tau, # tau: dosing interval
                               n_doses, # n_doses: number of doses
                               t, # t: time points
                               exposure_duration_hr) {
  # Initialize concentration vector
  concentration <- numeric(length(t))
  
  # Dynamically calculate t_alpha_end
  t_alpha_end <- if (!is.na(k_alpha) && !is.na(k_beta) && k_alpha > k_beta) {
    (log(k_alpha) - log(k_beta)) / (k_alpha - k_beta)
  } else {
    0  # If k_alpha or k_beta is invalid, set t_alpha_end to 0
  }
  
  # Loop over each dose administered
  for (i in 1:n_doses) {
    dose_time <- (i - 1) * tau  # Time at which the current dose is administered
    
    # Loop over each time point
    for (j in 1:length(t)) {
      if (t[j] >= dose_time) {  # Only calculate for time points after the dose is administered
        if (t[j] <= dose_time + t_alpha_end) {  # Alpha phase
          if (!is.na(kabs) && !is.na(k_alpha) && kabs != k_alpha) {
            exp_k_alpha <- exp(-k_alpha * (t[j] - dose_time))
            exp_kabs <- exp(-kabs * (t[j] - dose_time))
            concentration[j] <- concentration[j] + 
              (D_per_dose / VD2) * (kabs / (kabs - k_alpha)) * (exp_k_alpha - exp_kabs)
          }
        } else {  # Beta phase
          if (!is.na(k_beta)) {
            # Calculate concentration at the end of the alpha phase
            C_alpha_end <- if (!is.na(kabs) && !is.na(k_alpha) && kabs != k_alpha) {
              (D_per_dose / VD2) * (kabs / (kabs - k_alpha)) * 
                (exp(-k_alpha * t_alpha_end) - exp(-kabs * t_alpha_end))
            } else {
              0
            }
            # Beta phase concentration
            concentration[j] <- concentration[j] + 
              C_alpha_end * exp(-k_beta * (t[j] - dose_time - t_alpha_end))
          }
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
  
  # Sort time and concentration together
  sorted_indices <- order(time)
  time <- time[sorted_indices]
  concentration <- concentration[sorted_indices]
  
  # Apply trapezoidal rule
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
    uiOutput("deploymentDate"),
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
                                downloadButton("validation_data", "Full Validation Dataset"), 
                                br(),
                                br(),
                                fileInput("upload_experiment_data", "Upload Experimental Conditions Table (.csv or .xlsx)", accept = c(".csv", ".xlsx"))
                            )
                        )
                    ),
                    h4("Once experimental conditions data are entered above, please proceed by clicking on the", strong("Model Parameters tab"), "on the left side of the app.")
                )
              ),
              tabPanel(title = "Model Availability",
                       fluidRow(
                         box(title = "Searchable availability table", width = 12, status = "primary", solidHeader = TRUE,
                             DTOutput("availability_table")
                         )
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
                                  p("To view simulation results, go to the ", strong("Results"), "tab on the left side of the app.") 
                                ),
                                br(),
                                h3("View TK data for selected chemical-species-sex combinations"),
                                box(title = "TK Data for selections", width = 12, collapsible = T, collapsed = F,
                                    #plotlyOutput("TK_plotly")
                                    withSpinner(uiOutput("TK_plot.ui"))
                                    ),
                                box(width = 12, 
                                     div(style = "display: flex; justify-content: center; align-items: center;", # Center the button
                                column(width = 2, downloadButton("download_TKPlotly_widget", "Download Plotly", icon("download"), style=" font-size: 15px; padding: 10px 18px; color: #fff; background-color: #337ab7; border-color: #2e6da4"))
                                )
                                ),
                                br(),
                                h3("An interactive datatable of the parameters selected in your models are below:"),
                                box(title = "Selected TK Parameter Dataset", width = 12, collapsible = T,
                                    DTOutput("TK_used_DT")
                                    ),
                                br(),
                                h4("Expand the below table to view sources"),
                                box(title = "Selected TK Parameter Dataset (sources)", width = 12, collapsible = T, collapsed = T,
                                    DTOutput("TK_sources_used_DT")
                                ),
                                h3("Explore the full toxicokinetic dataset below:"),
                                box(title = "Full TK Dataset", width = 12, collapsible = T, collapsed = T,
                                    DTOutput("Full_TK_Datatable")
                                    )
                                ),
                       tabPanel(title = "PBPK Models",
                                p("This RShiny app allows users to run the optimized models from", a("Chou & Lin et al. (2019). Bayesian evaluation of a physiologically based pharmacokinetic (PBPK) model for perfluorooctane sulfonate (PFOS) to characterize the interspecies uncertainty between mice, rats, monkeys, and humans: Development and performance verification",
                                                                                                     href = "https://www.sciencedirect.com/science/article/pii/S016041201930203X"),
                                  "and", a("Cheng and Ng (2017). A Permeability-Limited Physiologically Based Pharmacokinetic (PBPK) Model for Perfluorooctanoic acid (PFOA) in Male Rats", href = "https://pubs.acs.org/doi/10.1021/acs.est.7b02602")),
                                p("Both of these model have only been optimized for", strong("males"), "and via an IV and oral route of exposure (PFOS model is oral-only). The model consists of four organ compartments (plasma, liver, kidney, and rest of body). For more information, please refer to their publications.
                                  Schematic of the PBPK model frameworks from these publications are shown below:"),
                                h3("PFOS PBPK Male Rat Model"),
                                br(),
                                tags$a(
                                  href = "https://www.sciencedirect.com/science/article/pii/S016041201930203X",
                                  target = "_blank",
                                  tags$img(src = "fig1.PNG", width = "70%", alt = "PBPK Framework")
                                ),
                                br(),
                                
                                h2("Customize PFOS PBPK Parameters"),
                                p("The PFOS PBPK model is fully customizable (however not recommended)."),
                                uiOutput("species_sex_params_grid"),  # Dynamically generate UI for each species' parameters
                                
                                br(),
                                h3("PFOA PBPK Male Rat Model"),
                                br(),
                                tags$a(
                                  href = "https://pubs.acs.org/doi/10.1021/acs.est.7b02602",
                                  target = "_blank",
                                  tags$img(src = "PFOA_PBPK.jpeg", width = "70%", alt = "PBPK Framework")
                                ),
                                p("Rat model structure.There are seven tissues including blood (B), liver (L), gut (G), kidney (K), muscle (M), adipose (A), and “rest of
body” (R). All tissues except blood contain both a vascular space (e.g., KB for kidney) and tissue space, the latter of which can be further divided into
two subcompartments: interstitial fluid (e.g., KF for kidney) and tissue (e.g., KT for kidney). Liver, gut, and kidney also contain bile, gut lumen, and
filtrate, respectively. Blood flow rate for each tissue is indicated (e.g., QBK is the blood flow rate to kidney). In compartments with albumin, α2μ-
globulin, or liver-type fatty acid binding protein (L-FABP), PFOA binding occurs. PFOA exchange among connected compartments can occur via
passive diffusion (all compartments) and active transport (kidney tissue (KT) and liver tissue (LT))."),
                                br(),
                                br(),
                                
                                # p("Additional PBPK models are available for mice for PFOS, PFOA, and PFHxS via the oral, nasal, dermal, and IV pathways based on the models described in", a("Zhu et al. (2023). Exploring Route-Specific Pharmacokinetics of PFAS in Mice by Coupling In Vivo tests and Physiologically Based Toxicokinetic Models.", href = "https://ehp.niehs.nih.gov/doi/10.1289/EHP11969")),
                                # h3("PFOS, PFOA, PFHxS Male Mouse Models"),
                                # tags$a(
                                #   href = "https://ehp.niehs.nih.gov/doi/10.1289/EHP11969",
                                #   target = "_blank",
                                #   tags$img(src = "ehp11969_f1.jpg", width = "70%", alt = "PBPK Framework")
                                # ),
                                # br(),
                                # br(),
                                box(
                                  title = "MassTransferPBPK (Fischer 2025) — Schematic",
                                  status = "warning", solidHeader = TRUE, width = 12,
                                  uiOutput("fischer_schematic")
                                ),
                                
                                box(
                                  title = "MassTransferPBPK (Fischer 2025) — Model synopsis",
                                  status = "warning", solidHeader = TRUE, width = 12,
                                  uiOutput("fischer_synopsis")
                                ),
                                
                                # Fishcer Model
                                box(
                                  title = "MassTransferPBPK (Fischer 2025) — Model Information",
                                  status = "warning",
                                  solidHeader = TRUE,
                                  width = 12,
                                  p("This permeability/active-transport PBPK model (i.e., Fischer et al., 2025 https://pubs.acs.org/doi/10.1021/acs.est.5c05473) is currently implemented for",
                                    tags$b("Male Mouse"),
                                    "and supports the PFAS listed below as found in the app’s Fischer parameter file."),
                                  p("To run this model, choose Model Type = ", tags$code("MassTransferPBPK"),
                                    ", Species = ", tags$code("Mouse"),
                                    ", Sex = ", tags$code("Male"),
                                    ". If a PFAS is not listed here, the app will flag it during validation."),
                                  p("Edit parameter values below. The model will use these edited values when you click ",
                                    tags$b("Run Simulation"), ". If no PFAS are listed, add a MassTransferPBPK row for Mouse/Male in the Experiment table."),
                                  DTOutput("massTransfer_info_table"),
                                  p("Unsupported selections will be flagged during validation or when you click Run.")
                                )
                                
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
                                
                                # --- Filters ------------------------------------------------------------
                                box(
                                  title = "Filter visuals",
                                  status = "primary",
                                  solidHeader = TRUE,
                                  width = 12,
                                  collapsible = TRUE,
                                  collapsed = FALSE,
                                  # Hide filters until we have results (prevents empty dropdowns)
                                  conditionalPanel("output.hasResults",
                                                   fluidRow(
                                                     column(3, selectizeInput(
                                                       "flt_species", "Species", choices = NULL, multiple = TRUE,
                                                       options = list(placeholder = "Select species…", plugins = list("remove_button"))
                                                     )),
                                                     column(3, selectizeInput(
                                                       "flt_pfas", "PFAS", choices = NULL, multiple = TRUE,
                                                       options = list(placeholder = "Select PFAS…", plugins = list("remove_button"))
                                                     )),
                                                     column(3, selectizeInput(
                                                       "flt_model", "Model", choices = NULL, multiple = TRUE,
                                                       options = list(placeholder = "Select model(s)…", plugins = list("remove_button"))
                                                     )),
                                                     column(3, selectizeInput(
                                                       "flt_comp", "Compartment", choices = NULL, multiple = TRUE,
                                                       options = list(placeholder = "Select compartment(s)…", plugins = list("remove_button"))
                                                     ))
                                                   ),
                                                   div(style = "margin-top:6px; color:#666;",
                                                       tags$small("Tip: start typing to search; use Backspace to remove a tag."))
                                  ),
                                  conditionalPanel("!output.hasResults",
                                                   tags$em("Run a simulation to enable filters.")
                                  )
                                ),
                                
                                # --- Concentration Plot -------------------------------------------------
                                box(
                                  title = "Concentration over time",
                                  status = "primary",
                                  solidHeader = TRUE,
                                  width = 12,
                                  # right-aligned download button
                                  div(style = "display:flex; justify-content:flex-end; gap:10px; margin-bottom:8px;",
                                      downloadButton("download_time_conc", "Download Concentration Data")
                                  ),
                                  withSpinner(plotlyOutput("concentrationPlot"))
                                ),
                                
                                # --- Modeled vs Measured -----------------------------------------------
                                box(
                                  title = "Modeled vs. measured concentration",
                                  status = "primary",
                                  solidHeader = TRUE,
                                  width = 12,
                                  div(style = "display:flex; justify-content:flex-end; gap:10px; margin-bottom:8px;",
                                      downloadButton("download_scatterPlot_widget", "Download Interactive Plotly",
                                                     icon = icon("download"),
                                                     style="font-size:15px; padding:10px 18px; color:#fff; background-color:#FFA500; border-color:#2e6da4")
                                  ),
                                  withSpinner(plotlyOutput("scatterPlot"))
                                ),
                                
                                # --- Data Table ---------------------------------------------------------
                                box(
                                  title = "Modeled and measured concentrations",
                                  status = "primary",
                                  solidHeader = TRUE,
                                  width = 12,
                                  withSpinner(DTOutput("concentration_at_time"))
                                )
                       ),
                       tabPanel(title = "Summary Results",
                                fluidRow(
                                  h3("Summary Values"),
                                  p("Values below for maximum concentration (C_max), area under the curve (AUC), and time-weighted average concentraiton (C_TWA) are in mg/L for serum"),
                                  withSpinner(DTOutput("summary_table")),
                                  withSpinner(plotlyOutput("heatmap"))
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
                href = "https://oehha.shinyapps.io/OEHHA-Data-Explorer/", 
                target = "_blank", 
                tags$span(
                  icon("link"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff;",
                  " Visit the Chemical Data Explorer Tool"
                )
              ),
              tags$span(
                tags$a(
                  href = "https://github.com/ScottCoffin/OEHHA-Data-Explorer/", 
                  target = "_blank",
                  icon("github"),
                  style = "font-size: 2rem; font-weight: bold; color: #007bff; margin-left: 10px;",
                  " View GitHub Repository"
                )
              ),
              br(),
              tags$a(
                href = "https://oehha.shinyapps.io/OEHHA-Data-Explorer/",
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
  
  # report dynamic deployment date
  output$deploymentDate <- renderUI({
    tags$h4(
      paste("App updated on", load_deployment_date()),
      style = "text-align: center; color: #D3D3D3;"
    )
  })
  
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
  
  # Download handler for the full validation CSV
  output$"validation_data" <- downloadHandler(
    filename = function() {
      "validation_data.csv"
    },
    content = function(file) {
      write.csv(validation_data, file, row.names = FALSE)
    }
  )
  
  # helper for UI
  output$hasResults <- shiny::reactive({
    df <- simulation_results()
    !is.null(df) && nrow(df) > 0
  })
  outputOptions(output, "hasResults", suspendWhenHidden = FALSE)
  
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
    
    # After reading `new_data` from csv/xlsx, normalize types
    num_keys <- c("Dose_mg_per_kg","Interval_Hours","Exposure_Duration_Days",
                  "Time_Serum_Collected_hr","Serum_Concentration_mg_L")
    
    new_data <- new_data %>%
      mutate(
        PFAS       = trimws(as.character(PFAS)),
        Species    = trimws(as.character(Species)),
        Sex        = trimws(as.character(Sex)),
        Model_Type = trimws(as.character(Model_Type))
      ) %>%
      mutate(across(all_of(num_keys), ~ suppressWarnings(as.numeric(.))))
    
    # remove any user-supplied columns that will collide with modeled outputs
    new_data <- new_data %>% select(-any_of(c("Compartment", "Concentration", "Time")))
    
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
    
    # prevent crashing with NA values
    updated_data <- drop_empty_experiment_rows(updated_data)   
    # Update the reactive value with the new data
    experiment_data(updated_data)
  })
  
  # Notify only: show unsupported choices and a suggested alternative; DO NOT mutate the table
  # Notify only: suggest a preferred model if (and only if) it's actually available for that combo
  observeEvent(experiment_data(), {
    ed <- drop_empty_experiment_rows(experiment_data())
    if (!nrow(ed)) return()
    
    avail <- isolate(model_availability())
    
    # Best available model per PFAS/Species/Sex, using your preference order
    best_tbl <- avail %>%
      dplyr::filter(Available, Model_Type %in% model_preference) %>%
      dplyr::group_by(PFAS, Species, Sex) %>%
      dplyr::summarise(
        best = Model_Type[which.min(pref_rank(Model_Type))],
        .groups = "drop"
      )
    
    # Join suggestions to the user’s rows; only message when the suggested model differs
    nudges <- ed %>%
      dplyr::left_join(best_tbl, by = c("PFAS","Species","Sex")) %>%
      dplyr::filter(!is.na(best), best != Model_Type) %>%
      dplyr::mutate(msg = sprintf("%s/%s/%s — consider %s (chosen: %s)",
                                  Species, Sex, PFAS, best, Model_Type)) %>%
      dplyr::pull(msg)
    
    if (length(nudges)) {
      showNotification(
        paste("More preferred model available:", paste(nudges, collapse = " • ")),
        type = "message", duration = 8
      )
    }
  })
  
  
  
  
  ########################### Model Availability #########
  # Make a full PFAS × Species × Sex grid and compute availability per model
  model_availability <- reactive({
    combos <- expand.grid(
      PFAS    = as.character(pfas_levels),
      Species = c("Rat","Mouse","Human","Monkey"),
      Sex     = c("Male","Female"),
      stringsAsFactors = FALSE
    )
    
    combos$`MassTransferPBPK` <- is_mass_transfer_supported(combos$PFAS, combos$Species, combos$Sex)
    combos$PBPK               <- is_pbpk_supported(combos$PFAS, combos$Species, combos$Sex)
    combos$`two-compartment`  <- is_simple_supported(combos$PFAS, combos$Species, combos$Sex, "two-compartment")
    combos$`single-compartment` <- is_simple_supported(combos$PFAS, combos$Species, combos$Sex, "single-compartment")
    # We intentionally do NOT rank "biphasic" because your stated order excludes it.
    # If you want to show availability only:
    combos$biphasic <- is_simple_supported(combos$PFAS, combos$Species, combos$Sex, "two-compartment") &
      !is.na(tk_params$K_alpha)[match(paste(combos$PFAS,combos$Species,combos$Sex,sep="|"),
                                      paste(tk_params$PFAS,tk_params$Species,tk_params$Sex,sep="|"))]
    
    long <- combos |>
      tidyr::pivot_longer(
        cols = c("MassTransferPBPK","PBPK","two-compartment","single-compartment","biphasic"),
        names_to = "Model_Type", values_to = "Available"
      )
    
    # Preferred = best-ranked among the four in model_preference only
    long <- long |>
      dplyr::group_by(PFAS, Species, Sex) |>
      dplyr::mutate(
        .best_rank = if (any(Available & Model_Type %in% model_preference)) {
          min(pref_rank(Model_Type[Available & Model_Type %in% model_preference]), na.rm = TRUE)
        } else NA_real_,
        Preferred = Available & (pref_rank(Model_Type) == .best_rank) & (Model_Type %in% model_preference)
      ) |>
      dplyr::ungroup() |>
      dplyr::select(-.best_rank)
    
    long
  })
  
  
  
  ### availability plotly ##
  output$availability_matrix_plot <- renderPlotly({
    df <- model_availability()
    shiny::req(nrow(df) > 0)
    
    # order PFAS for readability
    df <- df %>% arrange(PFAS)
    df$Model_Type <- factor(df$Model_Type, levels = model_preference)
    
    # base heat (available vs not)
    base <- ggplot(df, aes(x = Model_Type,
                           y = paste(PFAS, Sex, Species, sep=" • "),
                           fill = ifelse(Available, "yes", "no"),
                           text = paste0(
                             "PFAS: ", PFAS, "<br>",
                             "Species: ", Species, " &nbsp; Sex: ", Sex, "<br>",
                             "Model: ", Model_Type, "<br>",
                             "Available: ", ifelse(Available,"Yes","No"),
                             ifelse(Preferred & Available, "<br><b>Preferred ✅</b>", "")
                           ))) +
      geom_tile(color = "white", height = 0.9) +
      scale_fill_manual(values = c("no" = "#e9ecef", "yes" = "#a9dfbf"), guide = "none") +
      labs(x = NULL, y = "PFAS • Sex • Species") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    
    stars <- df %>% dplyr::filter(Preferred & Available)
    p <- base + geom_point(data = stars,
                           aes(x = Model_Type, y = paste(PFAS, Sex, Species, sep=" • ")),
                           shape = 8, size = 3, inherit.aes = FALSE)
    
    ggplotly(p, tooltip = "text")
  })
  
  ### availability DT ###
  output$availability_table <- renderDT({
    df <- model_availability()
    shiny::req(nrow(df) > 0)
    
    tbl <- df %>%
      dplyr::group_by(PFAS, Species, Sex) %>%
      dplyr::summarise(
        Preferred = {
          avail <- Model_Type[Available]
          if (length(avail)) avail[which.min(pref_rank(avail))] else NA_character_
        },
        Available_Models = {
          avail <- Model_Type[Available]
          if (length(avail)) paste(avail[order(pref_rank(avail))], collapse = "  •  ")
          else "—"
        },
        .groups = "drop"
      ) %>%
      dplyr::arrange(PFAS, Species, Sex) %>% 
      mutate(PFAS = as.factor(PFAS),
             Species = as.factor(Species),
             Sex = as.factor(Sex))
    
    # pretty badges using your model_type_colors
    badge_for <- function(label) {
      if (is.na(label) || label == "—") return("—")
      col <- unname(model_type_colors[label])
      sprintf('<span style="background:%s;color:#fff;border-radius:10px;padding:2px 8px; margin-right:4px; font-size:12px;">%s</span>', col, label)
    }
    tbl$Preferred <- vapply(tbl$Preferred, badge_for, character(1))
    tbl$Available_Models <- vapply(tbl$Available_Models, function(s) {
      if (s == "—") return("—")
      labs <- strsplit(s, "  •  ", fixed = TRUE)[[1]]
      paste(vapply(labs, badge_for, character(1)), collapse = " ")
    }, character(1))
    
    datatable(tbl, escape = FALSE, rownames = FALSE,
              filter = "top",
              options = list(pageLength=20, scrollX=TRUE,
                             dom='Blrtip', buttons=list(I('colvis'),'copy')))
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
    
    exp_data <- experiment_data() %>%
      dplyr::filter(Model_Type != "PBPK") %>%
      dplyr::distinct(PFAS, Species, Sex, Model_Type) %>%
      dplyr::mutate(dplyr::across(c(PFAS,Species,Sex,Model_Type), as.character))
    
    if (!nrow(exp_data)) {
      processed_params(exp_data)
      params_data(exp_data)   # clears the table
      return()
    }
    
    params_processed <- exp_data %>%
      dplyr::left_join(
        tk_params %>% dplyr::mutate(dplyr::across(c(PFAS,Species,Sex), as.character)),
        by = c("PFAS","Species","Sex")
      ) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        Half_Life_hr = dplyr::case_when(
          is.na(Half_Life_hr) & !is.na(Clearance_L_per_kg_d) & !is.na(Volume_of_Distribution_L_per_kg) ~
            (log(2) * Volume_of_Distribution_L_per_kg) / Clearance_L_per_kg_d,
          TRUE ~ Half_Life_hr
        ),
        Volume_of_Distribution_L_per_kg = dplyr::case_when(
          is.na(Volume_of_Distribution_L_per_kg) & !is.na(Clearance_L_per_kg_d) & !is.na(Half_Life_hr) ~
            Clearance_L_per_kg_d / (log(2) / Half_Life_hr),
          TRUE ~ Volume_of_Distribution_L_per_kg
        )
      ) %>%
      dplyr::arrange(dplyr::desc(PFAS), Species, Sex)
    
    # prune to current rows but keep user edits that still match
    params_data(reconcile_params(params_processed, params_data()))
    
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
    pd <- params_data()
    if (is.null(pd) || !nrow(pd)) {
      # show a simple note instead of leaving the old table around
      return(rhandsontable(data.frame(Note = "No simple-TK rows in the Experiment table."),
                           readOnly = TRUE))
    }
    
    # JS renderer unchanged if you want it:
    custom_renderer <- "
    function (instance, td, row, col, prop, value, cellProperties) {
      Handsontable.renderers.TextRenderer.apply(this, arguments);
      if (!isNaN(value) && typeof value === 'number') td.innerHTML = Number(value).toPrecision(4);
    }"
    
    rhandsontable(pd) %>%
      hot_cols(strict = TRUE, allowInvalid = FALSE, renderer = custom_renderer) %>%
      hot_col("PFAS",       readOnly = TRUE) %>%
      hot_col("Species",    readOnly = TRUE) %>%
      hot_col("Sex",        readOnly = TRUE) %>%
      hot_col("Model_Type", readOnly = TRUE)
  })
  
 
  # Observe changes in the rhandsontable and update reactive value
  observeEvent(input$params_table, {
    upd <- hot_to_r(input$params_table)
    key_cols <- c("PFAS","Species","Sex","Model_Type")
    if (nrow(upd)) upd <- dplyr::mutate(upd, dplyr::across(dplyr::all_of(key_cols), ~ as.character(.)))
    params_data(upd)
  })
  
  
  # create reactive df used in TK plotly and DT
  mismatch_data_norm_reactive <- reactive({
    # get combos in use for filtering plot #
    experiment_data <- experiment_data()
    
    needed_vals <- experiment_data %>% 
      distinct(PFAS, Species, Sex) %>% 
      mutate(combo = paste0(PFAS, Species, Sex))
    
    # Normalize the value within each variable (column) group
    mismatch_data_norm <- mismatch_data %>%
      mutate(Units = case_when(
        variable == "Volume of Distribution" ~ "L/kg",
        variable == "Clearance" ~ "mL hr^-1 kg^-1",
        variable == "Elimination Half-Life" ~ "hr",
        variable == "K_absorption" ~ "hr^-1",
        variable == "K (alpha phase)" ~ "hr^-1",
        variable == "K (beta phase)" ~ "hr^-1",
        variable == "Volume of Distribution (t=0)" ~ "L/kg",
        variable == "Volume of Distribution (beta phase)" ~ "L/kg"
      )) %>% 
      mutate(Species = species_name) %>% 
      mutate(combo = paste0(PFAS, Species, Sex)) %>% 
      filter(combo %in% needed_vals$combo) %>% 
      group_by(variable) %>%
      ungroup() %>%
      arrange(chem)  # Arrange by chemical
  
    
    #print
    mismatch_data_norm
  })
  
  #### create interactive DT for TK data used in selection ###
  output$TK_used_DT <- renderDT(server = F,{
    
    tk_df <- mismatch_data_norm_reactive() %>% 
      select(PFAS, Species, Sex, variable, value, Units) %>% 
      # concatenate variable with Units using parentheses
      mutate(variable = paste(variable, " (", Units, ")", sep = "")) %>%
      select(-Units) %>% 
      distinct() %>% 
      # pivot wider on variable with value
      pivot_wider(names_from = variable, values_from = value)
      print("Rendering Used TK datatable...")
      
    
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
    
    # Get unique species
    species_levels <- unique(tk_df$Species)
    
    # Dynamically truncate or recycle the color palette
    color_palette_final <- rep(color_palette, length.out = length(species_levels))
    text_palette_final <- rep(text_palette, length.out = length(species_levels))
    
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
    ) %>%  formatStyle(
        columns = "Species",
        target = "row",
        backgroundColor = styleEqual(species_levels, color_palette_final),
        color = styleEqual(species_levels, text_palette_final)
      )

    dt

         })
  
  #### create interactive DT for TK data (sources only) used in selection ###
  output$TK_sources_used_DT <- renderDT(server = F,{
    
    print(head(mismatch_data_norm_reactive()))
    
    tk_df <- mismatch_data_norm_reactive() %>% 
      select(PFAS, Species, Sex, variable, source, pubmedID, doi, Units) %>% 
      # concatenate variable with Units using parentheses
      mutate(variable = paste(variable, " (", Units, ")", sep = "")) %>%
      mutate(source = paste(source, " (PMID:", pubmedID , "; DOI:", doi, ")", sep = "")) %>% 
      select(-Units, - pubmedID, - doi) %>% 
      distinct() %>% 
      # pivot wider on variable with value
      pivot_wider(names_from = variable, values_from = source)
    print("Rendering Used TK Sources datatable...")
    
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
    
    # Get unique species
    species_levels <- unique(tk_df$Species)
    
    # Dynamically truncate or recycle the color palette
    color_palette_final <- rep(color_palette, length.out = length(species_levels))
    text_palette_final <- rep(text_palette, length.out = length(species_levels))
    
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
    ) %>%  formatStyle(
      columns = "Species",
      target = "row",
      backgroundColor = styleEqual(species_levels, color_palette_final),
      color = styleEqual(species_levels, text_palette_final)
    )
    
    dt
    
  })
  
  #### Plotly for TK Data that's actually powering the app! ###
  TK_plot <- reactive({
   
    #import reactive df from above
    mismatch_data_norm <- mismatch_data_norm_reactive()
      
    # Create the heatmap with the normalized value
    tk_heat_ggplot <- ggplot(mismatch_data_norm, 
                      aes(x = variable, 
                          y = paste(chem, sex, species_name, sep = " - "), 
                                              fill = value,
                                                      #normalized_value + 1e-12,
                                              text = paste("Chemical:", chem,"<br>",
                                                           "Species:", species_name, "<br>",
                                                           "Sex:", sex, "<br>",
                                                           "Variable:", variable, "<br>",
                                                           "Value:", value, Units, "<br>",
                                                           #"Match: ", mismatch_info, "<br>",
                                                           "Source:", source,"<br>",
                                                           "PubmedID:", pubmedID
                                              ))) +
      geom_tile(color = "white") +
      
      # Use the normalized value for the color scale
      scale_fill_gradient(trans = scales::log_trans(base = 10),
                          low = "red",
                          high = "#56B1F7",
                          space = "Lab",
                          na.value = "grey50"
      ) +
      
      # Labels and theme
      labs(x = "Toxicokinetic Parameters", 
           y = "Chemical - Sex - Species") +
      theme_minimal(base_size = 15) +
      theme(axis.text.y = element_text(size = 14),  # Adjust the y-axis text size
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            plot.title = element_blank(),
            legend.position = "none")
    
    
    # Normalize the data by variable (column)
    mismatch_data_norm <- mismatch_data_norm %>%
      group_by(variable) %>%
      mutate(
        normalized_value = (value - min(value, na.rm = TRUE)) / 
          (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
      ) %>%
      ungroup()
    
    # Define a custom color scale with three colors
    custom_colorscale <- list(
      list(0, "#2CA02C"),    # Low values
      list(0.5, "#1F77B4"), # Mid-range values
      list(1, "#E73F74")    # High values
    )

    # Create the heatmap using plot_ly
    tk_heat <- plot_ly(
      data = mismatch_data_norm,
      x = ~variable,  # Columns (x-axis)
      y = ~paste(chem, sex, species_name, sep = " - "),  # Rows (y-axis)
      z = ~normalized_value,  # Normalized values for the color scale
      type = "heatmap",
      colorscale = custom_colorscale,  # Apply the custom color scale
      text = ~paste(
        "Chemical:", chem, "<br>",
        "Species:", species_name, "<br>",
        "Sex:", sex, "<br>",
        "Variable:", variable, "<br>",
        "Value:", value, Units, "<br>",
        "Source:", source, "<br>",
        "PubmedID:", pubmedID
      ),
      hoverinfo = "text"  # Use the custom hover text
    )
    
    # Customize layout
    tk_heat <- tk_heat %>%
      layout(
        title = list(title = "Toxicokinetic Parameters Heatmap",
                     tickfont = list(size = 18)  # Increase x-axis text size
                     ),
        xaxis = list(title = "Toxicokinetic Parameters", tickangle = 45,
                     tickfont = list(size = 16)  # Increase x-axis text size
                     ),
        yaxis = list(title = "Chemical - Sex - Species",
                     tickfont = list(size = 16)  # Increase x-axis text size
                     ),
        colorbar = list(title = "Relative Value")
      )
    
    
    # Display the heatmap
    print(tk_heat)
  })
  
  # Create plotly 
#  output$TK_plotly <- renderPlotly({TK_plotly_react()})
  
  session_store <- reactiveValues()
  
  # Render TK heatmaply     
  output$TK_plotly <- renderPlotly({
    TK_plot <- TK_plot()
    
    # session_store$TK_plotly <- ggplotly(TK_plot, tooltip = "text") %>% 
    #   layout(legend = list(orientation = "h",   # Horizontal legend
    #                        x = 0.5,              # Center legend horizontally
    #                        xanchor = "center",   # Center the legend's position
    #                        y = 1.3))            # Position above the plot
    
    session_store$TK_plotly <- TK_plot
      
    print(session_store$TK_plotly)
  })
  
  
  # calculate TK heatmap chemical count for dynamic plotting
  TK_plot_count <- reactive({
    filtered_data <- experiment_data()
    
    needed_vals <- filtered_data %>% 
      distinct(PFAS, Species, Sex) %>% 
      mutate(combo = paste0(PFAS, Species, Sex))
    
    as.numeric(n_distinct(needed_vals$combo))
  })
  
  #assign  value for plot height based on column count (pixels)
  TK_plot_height <- reactive(100 * TK_plot_count())
  
  #render UI (dynamic height)
  output$TK_plot.ui <- renderUI({
    plotHeight <- TK_plot_height()
    plotlyOutput("TK_plotly", height = plotHeight, width = 1000)
  })
  
  #create download button for in TK heatmap
  output$downloadTKPlot <- downloadHandler(
    filename = function() {
      paste('TKHeatmap', Sys.Date(), '.', input$plotDevice_TK, sep='')
    },
    content = function(file) {
      ggsave(file, plot = TK_plot(), 
             width = input$plotWidth_silico, height = input$plotHeight_TK, device = input$plotDevice_TK)
    }
  )
  
  #download button for silico scatterplotly widget
  output$download_TKPlotly_widget <- downloadHandler(
    filename = function() {
      paste("TKHeatmaply-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$TK_plotly), file, selfcontained = TRUE)
    }
  )
  

  #### Create data object for plotting all TK data available (included data not used)
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
      filter(Endpoint %in% c("CL", "F", "HLe_invivo", "VDss", "kabs", "VD2", "VDc", "k_alpha", "k_beta")) %>% 
      mutate(Endpoint = case_when(
        Endpoint == "CL" ~ "Clearance",
        Endpoint == "F" ~ "Bioavailable Fraction",
        Endpoint == "HLe_invivo" ~ "Elimination Half-Life",
        Endpoint == "VDss" ~ "Volume of Distribution",
        Endpoint == "kabs" ~ "Absorption Kinetic Constant",
        Endpoint == "VD2" ~ "Volume of Distribution (alpha)",
        Endpoint == "VDc" ~ "Volume of Distribution (beta)",
        Endpoint == "k_alpha" ~ "Alpha-phase elimination rate constant",
        Endpoint == "k_beta" ~ "Beta-phase elimination rate constant"
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
  # --- Which PFAS are selected for MassTransferPBPK right now?
  mass_transfer_selected_pfas <- reactive({
    shiny::req(experiment_data())
    experiment_data() %>%
      dplyr::filter(Model_Type == "MassTransferPBPK") %>%
      # current Fischer function is Mouse/Male only
      dplyr::filter(Species == fischer_supported_species,
                    Sex     == fischer_supported_sex) %>%
      dplyr::pull(PFAS) %>%
      unique()
  })
  
  # --- Base param rows from the Fischer file for those PFAS
  mass_transfer_params_base <- reactive({
    if (is.null(fischer_params)) return(NULL)
    pf <- mass_transfer_selected_pfas()
    if (!length(pf)) return(NULL)
    df <- fischer_params %>% dplyr::filter(PFAS %in% pf)
    if (!nrow(df)) return(NULL)
    df
  })
  
  # --- Editable copy the user can change
  mass_transfer_params_edit <- reactiveVal(NULL)
  
  observeEvent(mass_transfer_params_base(), {
    mass_transfer_params_edit(mass_transfer_params_base())
  })
  
  current_supported_pfas <- reactive({
    df <- mass_transfer_params_edit()
    if (!is.null(df) && "PFAS" %in% names(df)) unique(df$PFAS) else fischer_supported_pfas
  })
  
  output$massTransfer_info_table <- renderDT({
    df <- mass_transfer_params_edit()
    if (is.null(df) || !nrow(df)) {
      return(DT::datatable(
        data.frame(
          Note = "Select at least one MassTransferPBPK row in the Experiment table (Mouse/Male only) to view/edit parameters."
        ),
        rownames = FALSE, options = list(dom = 't', pageLength = 5)
      ))
    }
    
    # lock the PFAS column from editing (we'll enforce it in the edit handler)
    DT::datatable(
      df,
      rownames = FALSE,
      editable = TRUE,
      options = list(scrollX = TRUE, pageLength = 20)
    )
  })
  
  # Apply user edits back into the reactiveVal
  observeEvent(input$massTransfer_info_table_cell_edit, {
    info <- input$massTransfer_info_table_cell_edit
    df <- isolate(mass_transfer_params_edit())
    
    colname <- colnames(df)[info$col + 1]  # DT is 0-indexed
    if (identical(colname, "PFAS")) return()  # don't allow changing PFAS labels
    
    # Coerce to numeric if the target column is numeric
    new_val <- info$value
    if (is.numeric(df[[colname]])) {
      new_val <- suppressWarnings(as.numeric(new_val))
    }
    df[info$row, colname] <- new_val
    mass_transfer_params_edit(df)
  })
  
  
  
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
    for (sp in na.omit(unique(experiment_data()$Species))) {
      for (sx in na.omit(unique(experiment_data()$Sex))) {
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
    if (is.null(experiment_data) || nrow(experiment_data) == 0) {
      return("No experiment data available.")
    }
    if (is.null(params_data) || nrow(params_data) == 0) {
      return("No parameters data available. Please upload or provide parameter data to run the simulation.")
    }
    
    experiment_data[] <- lapply(experiment_data, function(x) if (is.factor(x)) as.character(x) else x)
    params_data[]     <- lapply(params_data,     function(x) if (is.factor(x)) as.character(x) else x)
    
    # drop rows where ALL key fields are blank (belt-and-suspenders)
    req_keys <- c("Model_Type","PFAS","Species","Sex")
    experiment_data <- experiment_data %>%
      dplyr::filter(!if_all(all_of(req_keys), ~ is.na(.) | . == ""))
    
    safe_eq <- function(a, b) isTRUE(a == b)
    `%||%`   <- function(x, y) if (is.null(x)) y else x
    
    missing_rows <- lapply(seq_len(nrow(experiment_data)), function(row_idx) {
      row <- experiment_data[row_idx, ]
      
      # skip rows with any missing key (user is still editing)
      if (any(is.na(row[req_keys])) || any(row[req_keys] == "")) return(NULL)
      
      model_type <- row$Model_Type
      pfas  <- row$PFAS
      species <- row$Species
      sex     <- row$Sex
      required <- required_params[[model_type]] %||% character(0)
      
      if (safe_eq(model_type, "PBPK")) {
        if (!identical(pfas, "PFOS")) {
          return(paste(species, sex, pfas, "(PBPK supports PFOS only — choose PFOS or change model type)"))
        }
        param_table <- reactive_params[[paste0(species, "_", sex)]]
        if (is.null(param_table) || !pfas %in% colnames(param_table)) {
          return(paste(species, sex, pfas, "(PBPK parameter table missing for this species/sex)"))
        }
        return(NULL)
      }
      
      if (safe_eq(model_type, "MassTransferPBPK")) {
        if (!is_mass_transfer_supported(pfas, species, sex)) {
          return(paste0(
            species, " ", sex, " ", pfas,
            " (MassTransferPBPK allowed — PFAS: ", paste(fischer_supported_pfas, collapse = ", "),
            "; Species: ", fischer_supported_species,
            "; Sex: ", fischer_supported_sex, ")"
          ))
        }
        return(NULL)
      }
      
      # simple TK families
      matching_rows <- params_data[
        params_data$PFAS == pfas &
          params_data$Species == species &
          params_data$Sex == sex, , drop = FALSE
      ]
      if (nrow(matching_rows) == 0) {
        return(paste(species, sex, pfas, "(No parameter data)"))
      }
      missing_cols <- required[!required %in% colnames(matching_rows)]
      if (length(missing_cols) > 0) {
        return(paste(species, sex, pfas, "(Missing columns:", paste(missing_cols, collapse = ", "), ")"))
      }
      for (param in required) {
        if (any(is.na(matching_rows[[param]]))) {
          return(paste(species, sex, pfas, "(Missing values in:", param, ")"))
        }
      }
      NULL
    })
    
    missing_rows <- unlist(Filter(Negate(is.null), missing_rows))
    if (length(missing_rows) == 0) NULL else paste("Missing parameter data for:<br>", paste(missing_rows, collapse = "<br>"))
  }
  
  
  

  # Define dynamic required parameters
  required_params <- reactive({
    list(
      "PBPK"               = NULL,
      "MassTransferPBPK"   = NULL,  # params handled internally from Fischer file
      "single-compartment" = c("Volume_of_Distribution_L_per_kg", "Half_Life_hr"),
      "two-compartment"    = c("Volume_of_Distribution_L_per_kg", "Half_Life_hr", "Absorption_Coefficient_unitless"),
      "biphasic"           = c("Volume_of_Distribution_alpha_L_per_kg", "Volume_of_Distribution_beta_L_per_kg", "K_alpha", "K_beta", "Absorption_Coefficient_unitless")
    )
  })
  
  
  # Reactive value to store the validation message
  validation_message <- reactiveVal(NULL)
  
  # Update validation message dynamically
  observeEvent(list(experiment_data(), params_data()), {
    print("experiment_data or params_data changed")
    ed <- drop_empty_experiment_rows(experiment_data())
    pd <- params_data()
    
    
    # Compute validation result
    validation_result <- validate_params(
      experiment_data = ed,
      params_data = pd,
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
  
  #### Notifier to use preferred model
  # # Notify only: suggest a preferred model if (and only if) it's actually available for that combo
  # observeEvent(experiment_data(), {
  #   ed <- drop_empty_experiment_rows(experiment_data())
  #   if (!nrow(ed)) return()
  #   
  #   avail <- isolate(model_availability())
  #   
  #   # Best available model per PFAS/Species/Sex, using your preference order
  #   best_tbl <- avail %>%
  #     dplyr::filter(Available, Model_Type %in% model_preference) %>%
  #     dplyr::group_by(PFAS, Species, Sex) %>%
  #     dplyr::summarise(
  #       best = Model_Type[which.min(pref_rank(Model_Type))],
  #       .groups = "drop"
  #     )
  #   
  #   # Join suggestions to the user’s rows; only message when the suggested model differs
  #   nudges <- ed %>%
  #     dplyr::left_join(best_tbl, by = c("PFAS","Species","Sex")) %>%
  #     dplyr::filter(!is.na(best), best != Model_Type) %>%
  #     dplyr::mutate(msg = sprintf("%s/%s/%s — consider %s (chosen: %s)",
  #                                 Species, Sex, PFAS, best, Model_Type)) %>%
  #     dplyr::pull(msg)
  #   
  #   if (length(nudges)) {
  #     showNotification(
  #       paste("More preferred model available:", paste(nudges, collapse = " • ")),
  #       type = "message", duration = 8
  #     )
  #   }
  # })
  
  
  
  
  # Render the validation message in the UI
  output$params_check_message <- renderUI({
    ed <- drop_empty_experiment_rows(experiment_data())
    if (!nrow(ed)) return(NULL)
    
    chosen_ok <- with(ed, dplyr::case_when(
      Model_Type == "MassTransferPBPK" ~ is_mass_transfer_supported(PFAS, Species, Sex),
      Model_Type == "PBPK"             ~ is_pbpk_supported(PFAS, Species, Sex),
      Model_Type == "two-compartment"  ~ is_simple_supported(PFAS, Species, Sex, "two-compartment"),
      Model_Type == "single-compartment" ~ is_simple_supported(PFAS, Species, Sex, "single-compartment"),
      Model_Type == "biphasic"         ~ isTRUE(is_simple_supported(PFAS, Species, Sex, "two-compartment")), # conservative
      TRUE ~ FALSE
    ))
    
    if (all(chosen_ok, na.rm = TRUE)) {
      tags$p(HTML("All selected models are supported for those PFAS/species/sex."), style = "color: green; font-size: 16px; text-align: center;")
    } else {
      bad <- ed[!chosen_ok %in% TRUE, c("PFAS","Species","Sex","Model_Type")]
      msg <- paste(apply(bad, 1, paste, collapse = " | "), collapse = "<br>")
      tags$p(HTML(paste0("Unsupported selections (won't run):<br>", msg)),
             style = "color: red; font-size: 16px; text-align: center;")
    }
  })
  
    # Ensure the output is updated even when the tab is not active
  outputOptions(output, "params_check_message", suspendWhenHidden = FALSE)
  
  
  ### Observe and Force UI Re-rendering on Data Changes
  observeEvent(list(experiment_data(), params_data()), {
    output$params_check_message <- renderUI({
      ed <- drop_empty_experiment_rows(experiment_data())
      pd <- params_data()
      validation_result <- validate_params(
        experiment_data = ed,
        params_data     = pd,
        reactive_params = reactive_params,
        required_params = required_params()
      )
      if (is.null(validation_result)) {
        tags$p(HTML("All parameters are sufficiently available to run the simulation! <br> Please hit `Run Simulation` button below to proceed."),
               style = "color: green; font-size: 20px; text-align: center;")
      } else {
        tags$p(HTML(validation_result),
               style = "color: red; font-size: 20px; text-align: center;")
      }
    })
  })
  
  
  
  # Show the Fischer schematic
  output$fischer_schematic <- renderUI({
    # If you have permission to host a local copy, save it as: www/fischer2025_fig1.png
    # Otherwise, leave only the link (the <img> can be commented out).
    tags$figure(
      # Link to the article
      tags$a(
        href   = "https://doi.org/10.1021/acs.est.5c05473",
        target = "_blank",
        # Show a local image
        tags$img(
          src  = "acsest5c05473.JPG",  
          alt  = "Permeability-/transport-limited PBPK schematic with Blood/Plasma, Liver, Kidneys, Gut, and Rest compartments connected by passive PS and active transport pathways.",
          width = "80%"
        )
      ),
      tags$figcaption(
        HTML("Figure 1. Schematic of the MassTransferPBPK framework (adapted from Fischer et&nbsp;al., 2025, <i>Environmental Science &amp; Technology</i>, DOI: 10.1021/acs.est.5c05473).")
      )
    )
  })
  
  output$fischer_synopsis <- renderUI({
    # This uses your already-loaded `fischer_supported_pfas` vector
    pfas_list <- if (length(fischer_supported_pfas)) {
      paste(sort(fischer_supported_pfas), collapse = ", ")
    } else {
      "—"
    }
    
    HTML(paste0(
      "<ul style='margin-left:1em;'>",
      "<li><b>Model type:</b> permeability-/transport-limited PBPK (‘MassTransferPBPK’). ",
      "Compartments implemented in the app: Blood/Plasma, Liver, Kidneys, Gut, and Rest.</li>",
      "<li><b>Processes:</b> passive diffusion (permeability–surface terms) and active transport clearances ",
      "in relevant tissues; renal filtration and potential tubular reabsorption are represented; ",
      "oral dosing enters via gut; systemic mixing represented in blood/plasma.</li>",
      "<li><b>Inputs (in this app):</b> PFAS, dose (mg/kg), dosing interval (h), exposure duration (days). ",
      "Supported species/sex: Mouse (Male). Parameters are read from the Fischer parameter file and ",
      "can be edited in the table above.</li>",
      "<li><b>Outputs returned by <code>Fischer_PBPK_mouse()</code>:</b> Time series with ",
      "<code>time_h</code> (hours) and concentrations (mg/L) for <code>C_plasma</code>, ",
      "<code>C_blood</code>, <code>C_liver</code>, <code>C_kidneys</code>, <code>C_gut</code>, and <code>C_rest</code>. ",
      "The app maps these to ‘Plasma’, ‘Blood’, ‘Liver’, ‘Kidneys’, ‘Gut’, and ‘Rest’.</li>",
      "<li><b>Supported PFAS in your current parameter file:</b> ", pfas_list, ".</li>",
      "<li><b>Citation:</b> Fischer et&nbsp;al., 2025, <i>Environmental Science &amp; Technology</i>, ",
      "DOI: <a target='_blank' href='https://doi.org/10.1021/acs.est.5c05473'>10.1021/acs.est.5c05473</a>.</li>",
      "</ul>"
    ))
  })
  

#########################################################################################################
############################################### Run Simulation ##########################################
#####################################################################################################
  ## Helper for TK_params interactive table to update in real-time
  # helper: merge defaults with existing edits but only for currently selected combos
  reconcile_params <- function(defaults, edits) {
    keys <- c("PFAS","Species","Sex","Model_Type")
    
    to_char <- function(df) {
      if (is.null(df)) return(df)
      # ensure keys exist & are character
      for (k in keys) if (!k %in% names(df)) df[[k]] <- NA_character_
      dplyr::mutate(df, dplyr::across(dplyr::all_of(keys), ~ as.character(.)))
    }
    
    defaults <- to_char(defaults)
    edits    <- to_char(edits)
    
    if (is.null(defaults) || !nrow(defaults)) return(defaults)  # nothing selected → empty
    if (is.null(edits)    || !nrow(edits))    return(defaults)  # no edits yet → just defaults
    
    # only keep columns that exist in defaults (plus keys)
    edits <- dplyr::select(edits, dplyr::all_of(intersect(names(edits), names(defaults))))
    
    merged <- dplyr::left_join(defaults, edits, by = keys, suffix = c("", ".edit"))
    
    # coalesce edited values over defaults
    nonkey <- setdiff(intersect(names(defaults), names(edits)), keys)
    for (c in nonkey) {
      e <- paste0(c, ".edit")
      if (e %in% names(merged)) {
        merged[[c]] <- dplyr::coalesce(merged[[e]], merged[[c]])
        merged[[e]] <- NULL
      }
    }
    merged
  }
  
  observeEvent(experiment_data(), {
    ed <- drop_empty_experiment_rows(experiment_data())
    
    # simple-TK rows only (exclude PBPK)
    exp_keys <- ed %>%
      dplyr::filter(Model_Type != "PBPK") %>%
      dplyr::distinct(PFAS, Species, Sex, Model_Type) %>%
      dplyr::mutate(dplyr::across(c(PFAS,Species,Sex,Model_Type), as.character))
    
    if (!nrow(exp_keys)) {
      processed_params(exp_keys)  # empty
      params_data(exp_keys)       # <-- clear the handson table
      return()
    }
    
    # defaults from tk_params for current combos
    defaults <- exp_keys %>%
      dplyr::left_join(
        tk_params %>% dplyr::mutate(dplyr::across(c(PFAS,Species,Sex), as.character)),
        by = c("PFAS","Species","Sex")
      ) %>%
      dplyr::distinct() %>%
      dplyr::mutate(
        Half_Life_hr = dplyr::case_when(
          is.na(Half_Life_hr) & !is.na(Clearance_L_per_kg_d) & !is.na(Volume_of_Distribution_L_per_kg) ~
            (log(2) * Volume_of_Distribution_L_per_kg) / Clearance_L_per_kg_d,
          TRUE ~ Half_Life_hr
        ),
        Volume_of_Distribution_L_per_kg = dplyr::case_when(
          is.na(Volume_of_Distribution_L_per_kg) & !is.na(Clearance_L_per_kg_d) & !is.na(Half_Life_hr) ~
            Clearance_L_per_kg_d / (log(2) / Half_Life_hr),
          TRUE ~ Volume_of_Distribution_L_per_kg
        )
      ) %>%
      dplyr::arrange(dplyr::desc(PFAS), Species, Sex)
    
    processed_params(defaults)
    
    # merge in any previous edits, but prune to current combos
    params_data(reconcile_params(defaults, params_data()))
  })
  
  
  
  #### Helper for experiment table to handle NAs without crashing ###
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  drop_empty_experiment_rows <- function(df) {
    req <- c("PFAS","Species","Sex","Model_Type")
    if (!all(req %in% names(df))) return(df)
    df %>%
      dplyr::mutate(across(all_of(req), ~ if (is.factor(.)) as.character(.) else .)) %>%
      dplyr::filter(!if_all(all_of(req), ~ is.na(.) | . == "")) %>%  # remove rows with all 4 keys blank
      dplyr::mutate(
        PFAS       = factor(PFAS, levels = pfas_levels),
        Species    = factor(Species, levels = c("Rat","Mouse","Human","Monkey")),
        Sex        = factor(Sex, levels = c("Female","Male")),
        Model_Type = factor(Model_Type, levels = model_preference)
      )
  }
  
  
  
  #### Helper Function ####
  guess_time_days <- function(df, exposure_duration_days) {
    # prefer explicit hour columns
    if ("time_h" %in% names(df))      return(df$time_h / 24)
    if ("time_hours" %in% names(df))  return(df$time_hours / 24)
    if ("t_hr" %in% names(df))        return(df$t_hr / 24)
    # fallbacks
    if ("time" %in% names(df)) {
      t <- df$time
      # heuristic: if it runs far beyond exposure days, assume hours
      if (max(t, na.rm = TRUE) > exposure_duration_days * 1.5) return(t / 24)
      return(t)
    }
    if ("t" %in% names(df)) {
      tt <- df$t
      if (max(tt, na.rm = TRUE) > exposure_duration_days * 1.5) return(tt / 24)
      return(tt)
    }
    rep(NA_real_, nrow(df))
  }
  
  
  # 1) Make sure simulate_pfas() is sourced/available above this server code
  # source("R/simulate_pfas.R")  # if you saved it there
  
  simulation_results <- eventReactive(input$run, {
    shiny::req(experiment_data())
    
    # Prefer the edited table, but fall back to the processed defaults
    params <- params_data()
    if (is.null(params) || !nrow(params)) {
      params <- processed_params()
    }
    shiny::req(params)
    
    print("Running simulation...")
    
    exp_data <- experiment_data() %>%
      dplyr::mutate(
        Species = as.character(Species),
        PFAS = as.character(PFAS),
        Sex = as.character(Sex),
        Model_Type = as.character(Model_Type)
      )
    
    # Warn and exclude PBPK rows that are not PFOS (Chou and Lin model)
    bad_pbpk <- exp_data %>%
      dplyr::filter(Model_Type == "PBPK", PFAS != "PFOS")
    
    if (nrow(bad_pbpk) > 0) {
      shiny::showNotification(
        "PBPK (Chou & Lin 2019) supports PFOS only. Non-PFOS PBPK rows were skipped.",
        type = "warning", duration = 8
      )
      exp_data <- dplyr::anti_join(exp_data, bad_pbpk,
                                   by = c("PFAS","Species","Sex","Model_Type","Dose_mg_per_kg",
                                          "Interval_Hours","Exposure_Duration_Days","Time_Serum_Collected_hr",
                                          "Serum_Concentration_mg_L"))
    }
    
    # Warn & drop unsupported MassTransferPBPK rows before simulation
    bad_mass <- exp_data %>%
      dplyr::filter(Model_Type == "MassTransferPBPK") %>%
      dplyr::mutate(.ok = is_mass_transfer_supported(PFAS, Species, Sex)) %>%
      dplyr::filter(!.ok)
    
    if (nrow(bad_mass) > 0) {
      combos <- bad_mass %>%
        dplyr::mutate(label = paste0(Species, "/", Sex, "/", PFAS,
                                     " — allowed PFAS: ", paste(fischer_supported_pfas, collapse = ", "),
                                     "; species: ", fischer_supported_species,
                                     "; sex: ", fischer_supported_sex)) %>%
        dplyr::pull(label) %>% unique()
      
      shiny::showNotification(
        paste0("MassTransferPBPK rows skipped: ", paste(combos, collapse = " • ")),
        type = "warning", duration = 10
      )
      
      exp_data <- dplyr::anti_join(
        exp_data, bad_mass %>% dplyr::select(-.ok),
        by = c("PFAS","Species","Sex","Model_Type","Dose_mg_per_kg",
               "Interval_Hours","Exposure_Duration_Days","Time_Serum_Collected_hr",
               "Serum_Concentration_mg_L")
      )
    }
    
    
    
    print(head(params_data()))
    
    params <- params_data() %>%
      dplyr::mutate(
        Species = as.character(Species),
        PFAS = as.character(PFAS),
        Sex = as.character(Sex),
        Model_Type = as.character(Model_Type)
      )
    
    # Join experiment_data with params_data
    experiment_with_params <- exp_data %>%
      dplyr::left_join(params, by = c("Species", "Sex", "PFAS", "Model_Type"))
    
    # 2) Add "Compartment" so we can return multiple outputs per timepoint
    result_cols <- c(
      "Species", "PFAS", "Sex", "Model_Type", "Dose_mg_per_kg",
      "Exposure_Duration_Days", "Interval_Hours", "Compartment", "Time", "Concentration"
    )
    
    # Build a custom parameter CSV for MassTransferPBPK if the user edited anything
    custom_fischer_path <- NULL
    mt_overrides <- mass_transfer_params_edit()
    if (!is.null(mt_overrides) && nrow(mt_overrides)) {
      custom_fischer_path <- tempfile(fileext = ".csv")
      readr::write_csv(mt_overrides, custom_fischer_path)
    }
    
    
    # Run the simulation for each row in the joined data
    results <- lapply(seq_len(nrow(experiment_with_params)), function(i) {
      row <- experiment_with_params[i, ]
      pfas <- row$PFAS
      species <- row$Species
      sex <- row$Sex
      model_type <- row$Model_Type
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
      
      # MrGSolve styles
      if (model_type == "PBPK") {
        species_model <- species_models[[species]]
        pars <- species_model$params
        model <- species_model$model
        default_bw <- species_model$default_bw
        
        dose_mg <- default_bw * dose
        interval_between_doses <- 24 / dose_frequency_per_day
        total_doses <- exposure_duration_days * dose_frequency_per_day
        
        cat("Dose:", dose_mg, "mg\n")
        cat("Body Weight:", default_bw, "kg\n")
        cat("Interval:", interval, "hours\n")
        
        ex <- ev(ID = 1, amt = dose_mg, ii = interval_between_doses, addl = total_doses - 1, cmt = "AST")
        tgrid <- tgrid(0, exposure_duration_hrs + 96, by = 1) # model every 1 hour
        
        output <- tryCatch({
          model %>%
            param(10 ^ pars) %>% # get out of log10 space
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
            Compartment = "Plasma",
            Time = output$time / 24,
            Concentration = signif(output$Plasma, 4),
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            Species = species,
            PFAS = pfas,
            Sex = sex,
            Model_Type = model_type,
            Dose_mg_per_kg = dose,
            Exposure_Duration_Days = exposure_duration_days,
            Interval_Hours = interval,
            Compartment = "Plasma",
            Time = NA,
            Concentration = NA,
            stringsAsFactors = FALSE
          )
        }
        
        # Fischer Model
      } else if (model_type == "MassTransferPBPK") {
        
        # decide which param file to use
        param_file <- if (!is.null(custom_fischer_path)) custom_fischer_path else FISCHER_PARAM_PATH
        
        sim_df <- tryCatch({
          Fischer_PBPK_mouse(
            PFAS = pfas,
            dose_mg_per_kg = dose,
            exposure_duration_days = exposure_duration_days,
            interval_hours = interval,
            pfas_param_path = param_file  # <-- edited values flow through here
          )
        }, error = function(e) {
          message(paste("Fischer_PBPK_mouse error for Species:", species, pfas, "Message:", e$message))
          return(NULL)
        })
        
        if (is.null(sim_df)) {
          return(data.frame(
            Species = species, PFAS = pfas, Sex = sex, Model_Type = model_type,
            Dose_mg_per_kg = dose, Exposure_Duration_Days = exposure_duration_days,
            Interval_Hours = interval, Compartment = NA, Time = NA, Concentration = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        # --- guarantee a plasma column ---
        if (!"C_plasma" %in% names(sim_df)) {
          if ("Plasma" %in% names(sim_df)) {
            sim_df$C_plasma <- sim_df$Plasma
          } else if ("C_blood" %in% names(sim_df)) {
            # "central == plasma" per app decision
            sim_df$C_plasma <- sim_df$C_blood
          } else {
            sim_df$C_plasma <- NA_real_
          }
        }
        
        # ensure a days column for joining/plotting
        sim_df$Time_days <- guess_time_days(sim_df, exposure_duration_days)
        
        keep_cols <- intersect(
          c("C_blood", "C_plasma", "C_liver", "C_kidneys", "C_gut", "C_rest"),
          names(sim_df)
        )
        
        long <- tidyr::pivot_longer(
          sim_df,
          cols      = keep_cols,
          names_to  = "Compartment",
          values_to = "Concentration"
        )
        
        long$Compartment <- dplyr::recode(
          long$Compartment,
          C_blood = "Blood", C_plasma = "Plasma", C_liver = "Liver",
          C_kidneys = "Kidneys", C_gut = "Gut", C_rest = "Rest",
          .default = long$Compartment
        )
        
        out <- dplyr::tibble(
          Species = species,
          PFAS = pfas,
          Sex = sex,
          Model_Type = model_type,
          Dose_mg_per_kg = dose,
          Exposure_Duration_Days = exposure_duration_days,
          Interval_Hours = interval,
          Compartment = long$Compartment,
          Time = long$Time_days,                  # <- use robust days column
          Concentration = long$Concentration
        )
        
        out
        
      } else {
        # ----- simplified PK branches -----
        ke <- log(2) / row$Half_Life_hr
        ka <-  row$Absorption_Coefficient_unitless
        vd <-  row$Volume_of_Distribution_L_per_kg
        VD2 <- row$Volume_of_Distribution_alpha_L_per_kg
        VDc <- row$Volume_of_Distribution_beta_L_per_kg
        k_alpha <- row$K_alpha
        k_beta <- row$K_beta
        n_doses <- floor(exposure_duration_hrs / interval)
        times <- seq(0, exposure_duration_hrs + 96, by = 1)
        
        conc <- if (model_type == "biphasic" & !is.na(ka) & !is.na(VD2) & !is.na(VDc) & !is.na(k_alpha) & !is.na(k_beta)){
          purrr::map_dbl(times, ~ two_phase(
            D_per_dose = dose, VD2 = VD2, k_alpha = k_alpha, kabs = ka,
            k_beta = k_beta, VDc = VDc, tau = interval, n_doses = n_doses, t = .x,
            exposure_duration_hr = exposure_duration_hrs))
        } else if (model_type == "two-compartment") {
          purrr::map_dbl(times, ~ full_conc(dose, vd, ke, ka, interval, n_doses, .x, exposure_duration_hrs))
        } else {
          purrr::map_dbl(times, ~ simplified_conc(dose, vd, ke, interval, n_doses, .x, exposure_duration_hrs))
        }
        
        data.frame(
          Species = species,
          PFAS = pfas,
          Sex = sex,
          Model_Type = model_type,
          Dose_mg_per_kg = dose,
          Exposure_Duration_Days = exposure_duration_days,
          Interval_Hours = interval,
          Compartment = "Plasma",
          Time = times / 24,
          Concentration = conc,
          stringsAsFactors = FALSE
        )
      }
    })
    
    # Ensure uniform structure
    results <- lapply(results, function(res) {
      res[result_cols[!result_cols %in% names(res)]] <- NA
      res[result_cols]
    })
    
    do.call(rbind, results)
  })
  
  

  ###################################################### Results ###########################################  
  ####### HELPERS #####
  # PBPK (Chou & Lin) availability = PFOS only; species/sex must match your species_models list
  # PBPK: supports PFOS-only and sex/species exactly as defined in species_models
  is_pbpk_supported <- function(pfas, species, sex) {
    n <- length(pfas)
    out <- logical(n)
    sm_names <- names(species_models)
    for (i in seq_len(n)) {
      sp <- species[i]; sx <- sex[i]; pf <- pfas[i]
      if (is.na(sp) || is.na(sx) || is.na(pf) || !(sp %in% sm_names)) { out[i] <- FALSE; next }
      sm <- species_models[[sp]]
      out[i] <- identical(pf, sm$PFAS) && identical(sx, sm$Sex)  # PFOS + sex=Male in your list
    }
    out
  }
  
  # Simple TK availability: all required columns exist & are non-NA for the combo
  simple_required <- list(
    `single-compartment` = c("Volume_of_Distribution_L_per_kg","Half_Life_hr"),
    `two-compartment`    = c("Volume_of_Distribution_L_per_kg","Half_Life_hr","Absorption_Coefficient_unitless")
  )
  is_simple_supported <- function(pfas, species, sex, model) {
    req_cols <- simple_required[[model]]
    if (is.null(req_cols) || !all(req_cols %in% names(tk_params))) return(rep(FALSE, length(pfas)))
    keys_ok <- tk_params |>
      dplyr::filter(dplyr::if_all(dplyr::all_of(req_cols), ~ !is.na(.))) |>
      dplyr::mutate(.key = paste(PFAS, Species, Sex, sep="|")) |>
      dplyr::distinct(.key)
    (paste(pfas, species, sex, sep="|") %in% keys_ok$.key)
  }
  
  # Generic TK availability (= do we have the required, non-NA params in tk_params?)
  has_params_for <- function(pfas, species, sex, needs) {
    rows <- tk_params %>% dplyr::filter(PFAS==pfas, Species==species, Sex==sex)
    nrow(rows) > 0 && all(needs %in% names(rows)) && all(stats::complete.cases(rows[1, needs]))
  }
  
  
  # Robust normalizers
  normalize_text <- function(x) {
    x <- as.character(x)
    x <- gsub("[\u00A0\u2007\u202F]", " ", x, perl = TRUE)   # non-breaking spaces → space
    x <- trimws(x)
    x
  }
  
  normalize_pfas <- function(x) {
    x <- normalize_text(x)
    x <- toupper(x)
    x <- gsub("[^A-Z0-9+.-]", "", x)  # drop stray punctuation/Unicode
    # unify common synonyms (extend as needed)
    x <- dplyr::recode(x,
                       "HFPO-DA" = "GENX",
                       "HFPO_DA" = "GENX",
                       "HEXAFLUOROPROPYLENEOXIDEDIMERACID" = "GENX",
                       .default = x
    )
    x
  }
  
  norm_keys <- function(df) {
    force(df)
    chr <- intersect(c("Species","PFAS","Sex","Model_Type","Compartment"), names(df))
    for (c in chr) df[[c]] <- normalize_text(df[[c]])
    if ("PFAS" %in% names(df)) df$PFAS <- normalize_pfas(df$PFAS)
    df
  }
  
  time_key_days <- function(t_days) round(t_days * 24, 6) / 24  # nearest hour
  
  closest_join_multi <- function(exp_df, sim_df, comps = NULL, tol_hours = Inf) {
    gvars   <- c("Species","PFAS","Sex","Model_Type")
    numvars <- c("Dose_mg_per_kg","Exposure_Duration_Days","Interval_Hours")
    
    # sanity
    must <- c(gvars, "Time_key", "Compartment", "Time", "Concentration", numvars)
    miss <- setdiff(c("Time_key","Compartment","Time","Concentration"), names(sim_df))
    if (length(miss)) stop("sim_df missing columns: ", paste(miss, collapse = ", "))
    
    # which compartments to try
    if (is.null(comps)) comps <- sort(unique(sim_df$Compartment))
    
    exp_df <- exp_df %>% dplyr::mutate(exp_id = dplyr::row_number())
    
    purrr::map_dfr(seq_len(nrow(exp_df)), function(i) {
      e <- exp_df[i, ]
      
      # group by categories
      sgrp <- sim_df %>% dplyr::semi_join(e[gvars], by = gvars)
      if (!nrow(sgrp)) {
        # return one NA row per requested compartment
        return(purrr::map_dfr(comps, ~ dplyr::bind_cols(
          e,
          tibble::tibble(Compartment = .x, Time = NA_real_, Concentration = NA_real_)
        )))
      }
      
      # exact numeric keys (when present)
      exact <- sgrp
      for (v in numvars) {
        if (!is.na(e[[v]])) exact <- exact %>% dplyr::filter(.data[[v]] == e[[v]])
      }
      cand0 <- if (nrow(exact)) exact else sgrp
      
      # do one pick per compartment
      comps_avail <- intersect(comps, unique(cand0$Compartment))
      purrr::map_dfr(comps_avail, function(cp) {
        cand <- cand0 %>% dplyr::filter(Compartment == cp)
        if (!nrow(cand)) {
          return(dplyr::bind_cols(
            e, tibble::tibble(Compartment = cp, Time = NA_real_, Concentration = NA_real_)
          ))
        }
        
        dist <- abs(cand$Time_key - e$Time_key)
        if (is.finite(tol_hours)) {
          bad <- abs(cand$Time_key - e$Time_key) > tol_hours/24
          dist[bad] <- Inf
        }
        # add numeric-key distances (0 if exact/missing)
        for (v in numvars) {
          add <- if (!is.na(e[[v]])) abs(cand[[v]] - e[[v]]) else 0
          dist <- dist + add
        }
        
        k <- which.min(dist)
        if (!length(k) || is.infinite(dist[k])) {
          dplyr::bind_cols(
            e, tibble::tibble(Compartment = cp, Time = NA_real_, Concentration = NA_real_)
          )
        } else {
          dplyr::bind_cols(e, cand[k, c("Compartment","Time","Concentration")])
        }
      })
    }) %>% dplyr::select(-exp_id)
  }
  
  
  ############################### Concentration-time Plot ##############
  # Render concentration-time plot with dynamic xlim
  output$concentrationPlot <- renderPlotly({
    sim_results <- simulation_results()
    
    # guard: if no results yet
    if (is.null(sim_results) || !nrow(sim_results)) return(NULL)
    
    # apply the multi-select filters (treat empty as "all")
    species_sel <- if (is.null(input$flt_species) || !length(input$flt_species)) unique(sim_results$Species) else input$flt_species
    pfas_sel    <- if (is.null(input$flt_pfas)    || !length(input$flt_pfas))    unique(sim_results$PFAS)    else input$flt_pfas
    model_sel   <- if (is.null(input$flt_model)   || !length(input$flt_model))   unique(sim_results$Model_Type) else input$flt_model
    comp_sel    <- if (is.null(input$flt_comp)    || !length(input$flt_comp))    unique(sim_results$Compartment) else input$flt_comp
    
    sim_f <- sim_results %>%
      dplyr::filter(
        Species     %in% species_sel,
        PFAS        %in% pfas_sel,
        Model_Type  %in% model_sel,
        Compartment %in% comp_sel
      ) %>%
      dplyr::mutate(
        PFAS = as.factor(PFAS),
        Compartment = dplyr::recode(Compartment,
                                    Plasma="Plasma", Blood="Blood", 
                                    Liver="Liver", Kidneys="Kidneys", Gut="Gut", Rest="Rest",
                                    .default = as.character(Compartment)
        ),
        combo_line = interaction(Model_Type, Compartment, sep = " • ")
      )
    
    # build a stable linetype set (so legend order doesn’t jump)
    all_combo_levels <- sort(unique(interaction(
      c("PBPK","single-compartment","two-compartment","biphasic","MassTransferPBPK"),
      c("Plasma", "Blood","Liver","Kidneys","Gut","Rest"),
      sep = " • "
    )))
    base_lts <- c("solid","longdash","dotdash","twodash","dotted")
    combo_linetypes <- setNames(rep(base_lts, length.out = length(all_combo_levels)), all_combo_levels)
    
    p <- ggplot(
      sim_f,
      aes(
        x = Time, y = Concentration,
        color = PFAS,
        linetype = combo_line,
        group = interaction(Species, PFAS, Sex, Model_Type, Compartment,
                            Exposure_Duration_Days, Interval_Hours, Dose_mg_per_kg),
        text = paste(
          "Chemical:", PFAS, "<br>",
          "Species:", Species, "<br>",
          "Sex:", Sex, "<br>",
          "Model:", Model_Type, "<br>",
          "Compartment:", Compartment, "<br>",
          "Dose:", Dose_mg_per_kg, "mg/kg-d", "<br>",
          "Exposure:", Exposure_Duration_Days, "d<br>",
          "Interval:", Interval_Hours, "hrs<br>",
          "Time:", signif(Time, 4), "d<br>",
          "Conc:", signif(Concentration,3), "mg/L"
        )
      )
    ) +
      geom_line(size = 1) +
      scale_linetype_manual(values = combo_linetypes, drop = FALSE) +
      labs(
        title = "Concentration Over Time",
        x = "Time (days)", y = "Concentration (mg/L)",
        color = "PFAS", linetype = "Model • Compartment"
      ) +
      xlim(0, max(sim_f$Time, na.rm = TRUE) + 4) +
      theme_minimal(base_size = 13)
    
    plt <- ggplotly(p, tooltip = "text")
    session_store$concentrationPlot <- plt
    plt
  })
  
  
  
  ############################ Scatterplot of measured vs. modeled #############
  # Initialize choices from data after a run
  observeEvent(simulation_results(), {
    sim <- simulation_results()
    updateSelectizeInput(session, "flt_species", choices = sort(unique(sim$Species)), server = TRUE,
                         selected = sort(unique(sim$Species)))
    updateSelectizeInput(session, "flt_pfas",    choices = sort(unique(sim$PFAS)),    server = TRUE,
                         selected = sort(unique(sim$PFAS)))
    updateSelectizeInput(session, "flt_model",   choices = sort(unique(sim$Model_Type)), server = TRUE,
                         selected = sort(unique(sim$Model_Type)))
    updateSelectizeInput(session, "flt_comp",    choices = sort(unique(sim$Compartment)), server = TRUE,
                         selected = sort(unique(sim$Compartment)))
  })
  
  
  output$scatterPlot <- renderPlotly({
    sim_results <- simulation_results()
    if (is.null(sim_results) || !nrow(sim_results)) return(NULL)
    
    # --- fallbacks like your time-series plot ---
    species_sel <- if (is.null(input$flt_species) || !length(input$flt_species)) unique(sim_results$Species) else input$flt_species
    pfas_sel    <- if (is.null(input$flt_pfas)    || !length(input$flt_pfas))    unique(sim_results$PFAS)    else input$flt_pfas
    model_sel   <- if (is.null(input$flt_model)   || !length(input$flt_model))   unique(sim_results$Model_Type) else input$flt_model
    comp_sel    <- if (is.null(input$flt_comp)    || !length(input$flt_comp))    unique(sim_results$Compartment) else input$flt_comp
    
    # norm_keys() (which uppercases PFAS)
    pfas_sel <- normalize_pfas(pfas_sel)
    
    num_keys <- c("Dose_mg_per_kg","Exposure_Duration_Days","Interval_Hours")
    
  
    sim_f <- sim_results %>%
      norm_keys() %>%
      filter(
        Species %in% species_sel,
        PFAS %in% pfas_sel,
        Model_Type %in% model_sel
      ) %>%
      mutate(Time_key = time_key_days(Time)) %>%
      filter(Compartment == "Plasma") %>%
      mutate(across(all_of(num_keys), as.double))
    
    exp_data <- experiment_data() %>%
      filter(!is.na(Time_Serum_Collected_hr)) %>%
      norm_keys() %>%
      mutate(
        Time_in_days = Time_Serum_Collected_hr / 24,
        Time_key     = time_key_days(Time_in_days)
      ) %>%
      mutate(across(all_of(num_keys), as.double)) %>%
      # drop any uploaded columns that can collide with modeled ones
      select(-any_of(c("Compartment","Concentration","Time")))
    
    measured_conc <- closest_join_multi(exp_data, sim_f, comps = "Plasma", tol_hours = 24) %>%
      dplyr::mutate(
        modeled_concentration = Concentration,
        Legend_Label = paste(Species, PFAS, Model_Type, sep = " | ")
      )
    
    ### TEMP DEBUG
    print(
      measured_conc %>%
        dplyr::count(Species, Sex, PFAS, Model_Type, Compartment, wt = !is.na(Time)) %>%
        dplyr::arrange(PFAS, Model_Type, Compartment)
    )
    
    
    # guard: nothing matched
    if (!nrow(measured_conc) || all(!is.finite(measured_conc$modeled_concentration))) return(NULL)
    
    vals <- c(measured_conc$modeled_concentration, measured_conc$Serum_Concentration_mg_L)
    vals <- vals[is.finite(vals) & vals > 0]
    if (!length(vals)) return(NULL)
    
    expanded_limits <- c(
      10^(floor(log10(min(vals)))),
      10^(ceiling(log10(max(vals))))
    )
    
    
    
    p <- ggplot(measured_conc, aes(
      y = modeled_concentration, x = Serum_Concentration_mg_L,
      color = PFAS, shape = Model_Type,
      text = paste(
        "Chemical:", PFAS, "<br>",
        "Species:", Species, "<br>",
        "Sex:", Sex, "<br>",
        "Model:", Model_Type, "<br>",
        "Compartment: Plasma<br>",
        "Dose:", Dose_mg_per_kg, "mg/kg-d", "<br>",
        "Exposure:", Exposure_Duration_Days, "d<br>",
        "Interval:", Interval_Hours, "hrs<br>",
        "Measured:", signif(Serum_Concentration_mg_L,3), "mg/L<br>",
        "Modeled:", signif(modeled_concentration,3), "mg/L"
      )
    )) +
      geom_point(size = 2, alpha = 0.9, na.rm = TRUE) +
      scale_x_log10(limits = expanded_limits) +
      scale_y_log10(limits = expanded_limits) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray10") +
      geom_abline(slope = 1, intercept = 1/3,   linetype = "dashed", linewidth = 0.4, color = "gray50") +
      geom_abline(slope = 1, intercept = -1/3,  linetype = "dashed", linewidth = 0.4, color = "gray50") +
      labs(
        title = "Modeled vs. Measured Concentration",
        x = "Measured (mg/L)", y = "Modeled (mg/L)",
        color = "PFAS", shape = "Model"
      ) +
      theme_minimal(base_size = 15) +
      theme(legend.title = element_blank())
    
    session_store$scatterPlot <- ggplotly(p, tooltip = "text")
    session_store$scatterPlot
  })
  
  
  #download button for scatterplotly widget
  output$download_scatterPlot_widget <- downloadHandler(
    filename = function() {
      paste("scatterPlot-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$scatterPlot), file, selfcontained = TRUE)
    }
  )
  
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
    if (is.null(sim_results) || !nrow(sim_results)) return(NULL)
    
    num_keys <- c("Dose_mg_per_kg","Exposure_Duration_Days","Interval_Hours")
    
    
    sim_keyed <- sim_results %>%
      norm_keys() %>%
      mutate(
        Time_key = time_key_days(Time),
        across(all_of(num_keys), as.double)
      )
    
    exp_keyed <- experiment_data() %>%
      filter(!is.na(Time_Serum_Collected_hr)) %>%
      norm_keys() %>%
      mutate(
        Time_in_days = Time_Serum_Collected_hr / 24,
        Time_key     = time_key_days(Time_in_days),
        across(all_of(num_keys), as.double)
      ) %>%
      select(-any_of(c("Compartment","Concentration","Time")))
    
    data <- closest_join_multi(
      exp_keyed, sim_keyed,
      comps = c("Plasma","Blood","Liver","Kidneys","Gut","Rest"),
      tol_hours = Inf
    ) %>%
      dplyr::rename(
        "Predicted Concentration @ Hr Serum Collected (mg/L)" = Concentration,
        "Model" = Model_Type,
        "Dose (mg/kg)" = Dose_mg_per_kg,
        "Dosing Interval" = Interval_Hours,
        "Exposure (Days)" = Exposure_Duration_Days,
        "Hr Serum Collected" = Time_Serum_Collected_hr,
        "Measured Concentration (mg/L)" = Serum_Concentration_mg_L
      ) %>%
      dplyr::select(-Time_in_days, -Time_key, -Time) %>%
      dplyr::mutate(across(c(PFAS, Species, Sex, Model, Compartment), as.factor))
    

    
    dt <- datatable(
      data,
      rownames = FALSE,
      extensions = 'Buttons',
      filter = "top",
      options = list(
        pageLength = 25, autoWidth = TRUE, width = '100%', scrollX = TRUE,
        dom = 'Blrtip',
        buttons = list(
          I('colvis'), 'copy',
          list(
            extend = 'collection',
            buttons = list(
              list(extend = "csv",   filename = "page",
                   exportOptions = list(columns = ":visible", modifier = list(page = "current"))),
              list(extend = "excel", filename = "page", title = NULL,
                   exportOptions = list(columns = ":visible", modifier = list(page = "current")))
            ),
            text = 'Download current page'
          ),
          list(
            extend = 'collection',
            buttons = list(
              list(extend = "csv",   filename = "data",
                   exportOptions = list(columns = ":visible", modifier = list(page = "all"))),
              list(extend = 'excel', filename = "data", title = NULL,
                   exportOptions = list(columns = ":visible", modifier = list(page = "all")))
            ),
            text = 'Download all data'
          )
        ),
        # <-- keep lengthMenu INSIDE options
        lengthMenu = list(c(10, 30, 50, -1), c('10', '30', '50', 'All'))
      ),
      class = "display"   # <-- keep class as an arg to datatable()
    ) %>%
      formatStyle(
        target = 'row',
        backgroundColor = styleEqual(
          c("PBPK", "single-compartment", "two-compartment", "biphasic", "MassTransferPBPK"),
          c("#3A9AB2", "#6FB2C1", "#91BAB6", "#91BAD1", "#5A8356")
        ),
        columns = "Model"
      )
    
  
  # format numerics
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  for (col in numeric_cols) {
    dt <- dt %>% formatSignif(columns = col, digits = 3)
  }
  dt
  })


########################################## Calculate Summary Statistics ############  
  # Calculate summary statistics (C_max, C_TWA, AUC) for a specified time period
  summary_stats <- reactive({
    message("Calculating summary statistics (serum/plasma only)...")
    
    sim_results <- simulation_results()
    if (is.null(sim_results) || nrow(sim_results) == 0) return(NULL)
    
    # Keep *plasma/serum only* across all model types (MassTransferPBPK included)
    serum_only <- sim_results %>%
      dplyr::filter(Compartment == "Plasma")
    
    serum_only %>%
      # bring in exposure window per-row to truncate AUC correctly
      dplyr::left_join(
        experiment_data(),
        by = c("PFAS","Species","Sex","Model_Type",
               "Dose_mg_per_kg","Interval_Hours","Exposure_Duration_Days")
      ) %>%
      dplyr::group_by(Species, PFAS, Sex, Model_Type,
                      Dose_mg_per_kg, Interval_Hours, Exposure_Duration_Days) %>%
      dplyr::summarize(
        C_max = max(Concentration, na.rm = TRUE),
        # Time is in *days* everywhere in the app
        AUC_day = {
          ok <- is.finite(Time) & is.finite(Concentration) &
            Time >= 0 & Time <= Exposure_Duration_Days
          if (!any(ok)) NA_real_ else calculate_auc(Time[ok], Concentration[ok])
        },
        .groups = "drop"
      ) %>%
      # Report AUC in mg*hr/L (multiply day-integral by 24).
      # TWA should be mg/L, so use AUC_day / days.
      dplyr::mutate(
        AUC  = AUC_day * 24,                                # mg*hr/L
        C_TWA = dplyr::if_else(Exposure_Duration_Days > 0,
                               AUC_day / Exposure_Duration_Days,  # mg/L
                               NA_real_)
      ) %>%
      dplyr::select(-AUC_day) %>%
      dplyr::mutate(dplyr::across(c(C_max, AUC, C_TWA), ~ signif(.x, 3)))
  })
  

  ##### Summary Stats Table #####
  # Render the summary statistics table with export options
  output$summary_table <- renderDT({
    message("Rendering summary table...")
    
    summary <- summary_stats() %>%
      dplyr::rename(
        "Model" = Model_Type,
        "Dose (mg/kg)" = Dose_mg_per_kg,
        "Exposure (Days)" = Exposure_Duration_Days,
        "Maximum Serum Concentration (mg/L)" = C_max,
        "Area Under Curve Serum Concentration (mg·hr/L)" = AUC,
        "Time-Weighted Average Serum Concentration (mg/L)" = C_TWA
      )
    
    if (is.null(summary) || nrow(summary) == 0) return(NULL)
    
    # Use your global model_type_colors defined up top
    mt_cols_names  <- names(model_type_colors)
    mt_cols_values <- unname(model_type_colors)
    
    datatable(
      summary,
      rownames = FALSE,
      extensions = 'Buttons',
      filter = "top",
      options = list(
        pageLength = 25, autoWidth = TRUE, width = '100%', scrollX = TRUE,
        dom = 'Blrtip',
        buttons = list(
          I('colvis'), 'copy',
          list(extend = 'collection',
               buttons = list(
                 list(extend = "csv",   filename = "page",
                      exportOptions = list(columns = ":visible", modifier = list(page = "current"))),
                 list(extend = 'excel', filename = "page", title = NULL,
                      exportOptions = list(columns = ":visible", modifier = list(page = "current")))
               ),
               text = 'Download current page'),
          list(extend = 'collection',
               buttons = list(
                 list(extend = "csv",   filename = "data",
                      exportOptions = list(columns = ":visible", modifier = list(page = "all"))),
                 list(extend = 'excel', filename = "data", title = NULL,
                      exportOptions = list(columns = ":visible", modifier = list(page = "all")))
               ),
               text = 'Download all data')
        ),
        lengthMenu = list(c(10, 30, 50, -1), c('10', '30', '50', 'All'))
      ),
      class = "display"
    ) %>%
      formatStyle(
        target  = 'row',
        columns = "Model",
        backgroundColor = styleEqual(mt_cols_names, mt_cols_values)
      )
  })
  

#### Summary stats heatmap ####
  output$heatmap <- renderPlotly({
    params_heat <- summary_stats() %>%
      tidyr::pivot_longer(
        cols = c("C_max", "C_TWA", "AUC"),
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      dplyr::mutate(
        Metric = dplyr::recode(Metric,
                               C_max = "Maximum Serum Concentration",
                               AUC   = "Area Under the Curve (serum)",
                               C_TWA = "Time-Weighted Average Serum Concentration"
        ),
        Units = dplyr::case_when(
          Metric == "Area Under the Curve (serum)" ~ "mg·hr/L",
          TRUE                                     ~ "mg/L"
        )
      ) %>%
      dplyr::group_by(Metric) %>%
      dplyr::mutate(
        normalized_value = (Value - min(Value, na.rm = TRUE)) /
          (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(PFAS) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = Metric,
        y = paste(PFAS, Sex, Species, Model_Type, sep = " - "),
        group = interaction(Species, PFAS, Sex, Model_Type,
                            Exposure_Duration_Days, Interval_Hours, Dose_mg_per_kg),
        fill = normalized_value + 0.001,
        text = paste0(
          "Chemical: ", PFAS, "<br>",
          "Species: ", Species, "<br>",
          "Sex: ", Sex, "<br>",
          "Exposure: ", Exposure_Duration_Days, " d<br>",
          "Dose: ", Dose_mg_per_kg, " mg/kg<br>",
          "Model: ", Model_Type, "<br>",
          "Metric: ", Metric, "<br>",
          "Value: ", signif(Value, 3), " ", Units
        )
      )) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(
        trans = scales::log_trans(base = 10),
        low = "#56B1F7", high = "red4", space = "Lab", na.value = "grey50"
      ) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 12, angle = 90, hjust = 1),
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    ggplotly(params_heat, tooltip = "text")
  })
  
  
  
} # close server

shinyApp(ui = ui, server = server)

