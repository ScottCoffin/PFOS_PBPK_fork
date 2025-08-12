# -----------------------------------------------
# PBPK Model for PFAS — Full Working Example
# -----------------------------------------------

# --- Load Libraries ---
library(readxl)
library(dplyr)
library(stringr)
library(deSolve)

# --- Helper: Clean and Parse P-values from "x ± y" or "no hit" ---
parse_p <- function(p) {
  if (is.na(p) || grepl("no hit", p, ignore.case = TRUE)) return(0)
  p_clean <- gsub("±.*", "", p)
  p_clean <- gsub("×10", "e", p_clean)
  return(as.numeric(p_clean))
}

# --- Physiological Parameters (Simplified Demo) ---
load_physiological_parameters <- function() {
  list(
    # Blood volume fractions (example values, update as needed)
    V_water_blood = 0.6,
    V_SA_blood = 0.1,
    V_Glob_blood = 0.1,
    V_SP_blood = 0.1,
    V_ML_blood = 0.1,
    V_sorb_blood = 0.4,  # SA + Glob + ML + SP
    
    # Liver
    V_water_liver = 0.7,
    V_FABP_liver = 0.05,
    V_SA_liver = 0.05,
    V_SP_liver = 0.05,
    V_ML_liver = 0.15,
    V_sorb_liver = 0.3,
    V_liver_tissue = 1,
    
    # Gut
    V_water_gut = 0.7,
    V_FABP_gut = 0.05,
    V_SA_gut = 0.05,
    V_SP_gut = 0.05,
    V_ML_gut = 0.15,
    V_sorb_gut = 0.3,
    V_gut_tissue = 1,
    
    # Kidneys
    V_water_kidneys = 0.7,
    V_FABP_kidneys = 0.05,
    V_SA_kidneys = 0.05,
    V_SP_kidneys = 0.05,
    V_ML_kidneys = 0.15,
    V_sorb_kidneys = 0.3,
    V_kidneys_tissue = 1,
    
    # Rest of Body
    V_water_rest = 0.7,
    V_FABP_rest = 0.05,
    V_SA_rest = 0.05,
    V_SP_rest = 0.05,
    V_ML_rest = 0.15,
    V_sorb_rest = 0.3,
    V_rest_tissue = 1,
    
    # Bile
    V_water_bile = 0.9,
    V_SA_bile = 0.05,
    V_ML_bile = 0.05,
    V_bile = 1,
    
    # Blood and gut lumen volumes
    V_blood = 1.0068,
    V_gut_lumen = 3.15
  )
}

# --- Partition Coefficient Calculation ---
calculate_partition_coeffs <- function(physio, compound) {
  with(compound, {
    with(physio, {
      K_blood <- (V_water_blood / V_blood) +
        (V_SA_blood / V_blood) * K_SA +
        (V_Glob_blood / V_blood) * K_Glob +
        (V_SP_blood / V_blood) * K_SP +
        (V_ML_blood / V_blood) * K_ML
      
      K_liver <- (V_water_liver / V_liver_tissue) +
        (V_FABP_liver / V_liver_tissue) * K_FABP +
        (V_SA_liver / V_liver_tissue) * K_SA +
        (V_SP_liver / V_liver_tissue) * K_SP +
        (V_ML_liver / V_liver_tissue) * K_ML
      
      K_bile <- (V_water_bile / V_bile) +
        (V_SA_bile / V_bile) * K_SA +
        (V_ML_bile / V_bile) * K_ML
      
      K_gut <- (V_water_gut / V_gut_tissue) +
        (V_FABP_gut / V_gut_tissue) * K_FABP +
        (V_SP_gut / V_gut_tissue) * K_SP +
        (V_ML_gut / V_gut_tissue) * K_ML +
        (V_SA_gut / V_gut_tissue) * K_SA
      
      K_kidneys <- (V_water_kidneys / V_kidneys_tissue) +
        (V_SA_kidneys / V_kidneys_tissue) * K_SA +
        (V_FABP_kidneys / V_kidneys_tissue) * K_FABP +
        (V_SP_kidneys / V_kidneys_tissue) * K_SP +
        (V_ML_kidneys / V_kidneys_tissue) * K_ML
      
      K_rest <- (V_water_rest / V_rest_tissue) +
        (V_FABP_rest / V_rest_tissue) * K_FABP +
        (V_SP_rest / V_rest_tissue) * K_SP +
        (V_ML_rest / V_rest_tissue) * K_ML +
        (V_SA_rest / V_rest_tissue) * K_SA
      
      list(
        K_blood = K_blood,
        K_liver = K_liver,
        K_bile = K_bile,
        K_gut = K_gut,
        K_kidneys = K_kidneys,
        K_rest = K_rest
      )
    })
  })
}

# --- Free/Bound Fraction Calculation ---
calculate_free_bound_fractions <- function(physio, K_values) {
  with(physio, {
    with(K_values, {
      f_free_blood    <- 1 / (1 + K_blood   * V_sorb_blood   / V_water_blood)
      f_free_liver    <- 1 / (1 + K_liver   * V_sorb_liver   / V_water_liver)
      f_free_gut      <- 1 / (1 + K_gut     * V_sorb_gut     / V_water_gut)
      f_free_kidneys  <- 1 / (1 + K_kidneys * V_sorb_kidneys / V_water_kidneys)
      f_free_rest     <- 1 / (1 + K_rest    * V_sorb_rest    / V_water_rest)
      f_free_filtrate <- f_free_blood
      
      list(
        f_free_blood = f_free_blood,
        f_bound_blood = 1 - f_free_blood,
        f_free_liver = f_free_liver,
        f_bound_liver = 1 - f_free_liver,
        f_free_gut = f_free_gut,
        f_bound_gut = 1 - f_free_gut,
        f_free_kidneys = f_free_kidneys,
        f_bound_kidneys = 1 - f_free_kidneys,
        f_free_rest = f_free_rest,
        f_bound_rest = 1 - f_free_rest,
        f_free_filtrate = f_free_filtrate
      )
    })
  })
}

# --- Master Initializer ---
initialize_pfas_model <- function(pfas_name = "PFOS",
                                  file_path = "es5c05473_si_001.xlsx",
                                  physio,
                                  bw = 30, #g
                                  total_IV_dose_mg_kg = 0,
                                  total_oral_dose_mg_kg = 1000) {
  
  S7 <- read_excel(file_path, sheet = "Table_S7") %>% filter(PFAS == pfas_name)
  S9 <- read_excel(file_path, sheet = "Table_S9") %>% filter(PFAS == pfas_name)
  S10 <- read_excel(file_path, sheet = "Table_S10") %>% filter(PFAS == pfas_name)
  
  # Extract logK values
  K_FABP <- 10^as.numeric(str_extract(S7$logK_FABP_W, "[0-9.]+"))
  K_SA   <- 10^as.numeric(str_extract(S7$logK_SA_w, "[0-9.]+"))
  K_Glob <- 10^as.numeric(str_extract(S7$logK_Glob_W, "[0-9.]+"))
  K_SP   <- 10^as.numeric(str_extract(S7$logK_SP_W, "[0-9.]+"))
  K_ML   <- 10^as.numeric(str_extract(S7$logK_PL_w, "[0-9.]+"))
  
  # Extract permeability and transport
  P_app      <- parse_p(S9$P_app)
  P_OATP1B2  <- parse_p(S9$P_OATP1B2)
  P_OATP2B1  <- parse_p(S9$P_OATP2B1)
  P_OAT1     <- parse_p(S10$P_OAT1)
  P_OAT2     <- parse_p(S10$P_OAT2)
  P_OAT3     <- parse_p(S10$P_OAT3)
  P_NTCP     <- parse_p(S10$P_NTCP)
  
  compound <- list(K_FABP = K_FABP, K_SA = K_SA, K_Glob = K_Glob,
                   K_SP = K_SP, K_ML = K_ML)
  
  K_values <- calculate_partition_coeffs(physio, compound)
  f_values <- calculate_free_bound_fractions(physio, K_values)
  
  total_IV_dose   <- total_IV_dose_mg_kg * bw / 1000
  total_oral_dose <- total_oral_dose_mg_kg * bw / 1000
  total_dose      <- total_IV_dose + total_oral_dose
  
  C_blood_t0       <- total_IV_dose / physio$V_blood
  C_free_blood_t0  <- C_blood_t0 * f_values$f_free_blood * physio$V_blood / physio$V_water_blood
  C_bound_blood_t0 <- C_blood_t0 * f_values$f_bound_blood * physio$V_blood / physio$V_sorb_blood
  C_gut_lumen_t0   <- total_oral_dose / physio$V_gut_lumen
  
  yini <- c(
    C_free_blood = C_free_blood_t0,
    C_bound_blood = C_bound_blood_t0,
    C_gut_lumen = C_gut_lumen_t0
  )
  
  return(list(
    yini = yini,
    K_values = K_values,
    f_values = f_values,
    P_values = list(
      P_app = P_app,
      P_OATP1B2 = P_OATP1B2,
      P_OATP2B1 = P_OATP2B1,
      P_OAT1 = P_OAT1,
      P_OAT2 = P_OAT2,
      P_OAT3 = P_OAT3,
      P_NTCP = P_NTCP
    )
  ))
}

# --- Run the Model for PFOS ---
physio <- load_physiological_parameters()

model_data <- initialize_pfas_model(
  pfas_name = "PFOS",
  file_path = "Additional files/Code/Fischer 2025/es5c05473_si_001.xlsx",
  physio = physio
)

print(model_data$yini)
print(model_data$K_values)

##### ODE ########
rigidode <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    # Example: Passive absorption from gut lumen
    J_absorption <- P_app * C_gut_lumen
    
    # Example: transfer from gut lumen to blood
    dC_gut_lumen <- -J_absorption
    dC_free_blood <- J_absorption / V_blood  # Assume instant mix for demo
    
    # Example: Bound/free equilibrium (you can add more detail)
    dC_bound_blood <- 0  # Placeholder
    
    # Add remaining compartments as needed
    
    # Return list of derivatives
    list(c(
      dC_free_blood,
      dC_bound_blood,
      dC_gut_lumen
    ))
  })
}


run_simulation <- function(yini, parms, times = seq(0, 65*86400, by = 3600)) {
  out <- ode(
    y = yini,
    times = times,
    func = rigidode,
    parms = parms,
    method = "lsoda"
  )
  as.data.frame(out)
}

# Build parameters list for ODEs
parms <- c(
  model_data$P_values,
  model_data$K_values,
  model_data$f_values,
  V_blood = physio$V_blood,
  V_gut_lumen = physio$V_gut_lumen
)

# Run simulation
sim_results <- run_simulation(
  yini = model_data$yini,
  parms = parms
)

# Plot gut lumen and blood concentrations
plot(sim_results$time / 86400, sim_results$C_gut_lumen,
     type = "l", col = "blue", ylab = "Concentration", xlab = "Days", ylim = c(0, max(sim_results$C_gut_lumen)))
lines(sim_results$time / 86400, sim_results$C_free_blood, col = "red")
legend("topright", legend = c("Gut Lumen", "Free Blood"), col = c("blue", "red"), lty = 1)

### Get final concentrations
get_final_tissue_concentrations <- function(sim_results, physio) {
  # Extract last row
  final <- sim_results[nrow(sim_results), ]
  
  # Total concentrations (mg/mL) in each compartment
  tissue_concs <- data.frame(
    Compartment = c("Blood", "Liver", "Kidneys", "Gut", "Rest"),
    Total_Concentration = c(
      (final$C_free_blood * physio$V_water_blood + final$C_bound_blood * physio$V_sorb_blood) / physio$V_blood,
      (final$C_free_liver * physio$V_water_liver + final$C_bound_liver * physio$V_sorb_liver) / physio$V_liver,
      (final$C_free_kidneys * physio$V_water_kidneys + final$C_bound_kidneys * physio$V_sorb_kidneys) / physio$V_kidneys,
      (final$C_free_gut * physio$V_water_gut + final$C_bound_gut * physio$V_sorb_gut) / physio$V_gut,
      (final$C_free_rest * physio$V_water_rest + final$C_bound_rest * physio$V_sorb_rest) / physio$V_rest
    )
  )
  
  return(tissue_concs)
}

# example usage
tissue_concs_final <- get_final_tissue_concentrations(sim_results, physio)
print(tissue_concs_final)

ggplot(tissue_concs_final, aes(x = Compartment, y = Total_Concentration)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Final Total Tissue Concentrations", y = "Concentration (mg/mL)") +
  theme_minimal()

