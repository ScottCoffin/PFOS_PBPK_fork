# Load libraries ----------------------------------------------------------
library(mrgsolve)
library(magrittr)
library(dplyr)
library(ggplot2)

# Load the best-fit parameters --------------------------------------------
# Replace with the actual file path where the MCMC results or best-fit parameters are stored.
MCMC <- readRDS("Additional files/Results/Workplace/rat.MCMC.rds")
best_fit_params <- MCMC$bestpar

# Define the prediction function ------------------------------------------
predict_plasma_concentration <- function(dose_mg_per_kg, interval_hours = 24, exposure_duration_days = 7, 
                                         dose_fraction = 1, dose_frequency_per_day = 1) {
  # Get out of log domain
  pars <- lapply(best_fit_params, exp)
  
  # Body weight of the rat (kg)
  BW <- 0.3
  
  # Calculate the dose per administration (fractional dose)
  dose_mg <- dose_mg_per_kg * BW * dose_fraction
  
  # Calculate the interval between doses
  interval_between_doses <- 24 / dose_frequency_per_day
  
  # Total number of additional doses
  total_doses <- exposure_duration_days * dose_frequency_per_day
  
  # Create an exposure event
  ex <- ev(ID = 1, amt = dose_mg, ii = interval_between_doses, addl = total_doses - 1, cmt = "AST", replicate = FALSE)
  
  # Set up the time grid for simulation
  tgrid <- tgrid(0, interval_hours * (exposure_duration_days - 1) + 24 * 100, 1)
  
  # Get the model output
  output <- mod %>%
    param(pars) %>%
    Req(Plasma) %>%
    update(atol = 1E-8, maxsteps = 5000) %>%
    mrgsim_d(data = ex, tgrid = tgrid)
  
  # Return the plasma concentration over time
  return(data.frame(Time = output$time / 24, Concentration = output$Plasma))
}

# Example usage -----------------------------------------------------------
# Predict plasma concentration for a dose of 50 mg/kg, with a half-dose every 12 hours, over 30 days
prediction <- predict_plasma_concentration(dose_mg_per_kg = 50, 
                                           dose_fraction = 0.5,  # Half-dose
                                           dose_frequency_per_day = 2,  # Every 12 hours
                                           exposure_duration_days = 30)

# Plot the prediction
ggplot(prediction, aes(x = Time, y = Concentration)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Predicted Plasma Concentration", x = "Time (days)", y = "Concentration (ug/ml)")

