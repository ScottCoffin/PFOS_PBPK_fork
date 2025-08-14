Fischer_PBPK_mouse <- function(
    PFAS,
    dose_mg_per_kg,                 # repeated oral dose per event (mg/kg); set 0 to disable
    exposure_duration_days,
    interval_hours,                  # hours between repeated oral doses; ignored if dose_mg_per_kg == 0
    # dosing options
    single_oral_dose_mg_per_kg = 0, # single bolus oral dose at t = 0 (mg/kg)
    constant_oral_mg_per_kg_day = 0,# constant oral intake rate (mg/kg/day)
    constant_iv_mg_per_kg_day = 0,  # constant IV intake rate (mg/kg/day)
    pfas_param_path = "Additional files/Datasets/Fischer/pfas_parameters.csv",
    extend_hours_post = 96 #hrs to extend following final dose
) {
  # --- Dependencies (namespaced for Shiny safety) ---
  if (!requireNamespace("deSolve", quietly = TRUE)) stop("Package 'deSolve' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))   stop("Package 'tidyr' is required.")
  
  # --- Load PFAS parameters ---
  stopifnot(file.exists(pfas_param_path))
  params_df <- utils::read.csv(pfas_param_path, stringsAsFactors = FALSE)
  pfas_params <- params_df[params_df$PFAS == PFAS, , drop = FALSE]
  if (nrow(pfas_params) == 0) stop("PFAS '", PFAS, "' not found in parameter table: ", pfas_param_path)
  
  # Chemical descriptors (convert from exp10 columns)
  P_app     <- pfas_params$P_app
  K_FABP    <- 10^pfas_params$K_FABP_exp
  K_SA      <- 10^pfas_params$K_SA_exp
  K_Glob    <- 10^pfas_params$K_Glob_exp
  K_SP      <- 10^pfas_params$K_SP_exp
  K_ML      <- 10^pfas_params$K_ML_exp
  P_OATP1B2 <- pfas_params$P_OATP1B2
  P_OATP2B1 <- pfas_params$P_OATP2B1
  P_NTCP    <- pfas_params$P_NTCP
  P_OAT1    <- pfas_params$P_OAT1
  P_OAT2    <- pfas_params$P_OAT2
  P_OAT3    <- pfas_params$P_OAT3
  
  # --- Simulation horizon & dosing schedule ---
  t_end_sec <- exposure_duration_days * 86400 + extend_hours_post * 3600 #extend 96hr or specific time post
  times_sec <- seq(0, t_end_sec, by = 3600)  # 1-hour resolution
  
  # Body weight
  bw_g  <- 30 #body weight (g)
  bw_kg <- bw_g / 1000
  
  # --- Physiological constants ---
  V_body <- 30 #body volume (cm^3)
  Q_blood_liver      <- 0.021
  Q_blood_gut        <- 0.025
  Q_blood_kidneys    <- 0.0217
  Q_blood_rest       <- 0.0717
  Q_blood_liver_port <- 0.009
  
  # Liver
  V_liver_tissue    <- 1.65
  V_blood_liver_hep <- 0.357
  V_bile            <- 0.03
  V_blood_liver     <- 0.153
  V_liver           <- V_liver_tissue + V_blood_liver + V_blood_liver_hep
  A_liver_hep       <- 700
  A_liver           <- 300
  
  V_SA_blood_liver    <- 0.0147 * V_blood_liver
  V_Glob_blood_liver  <- 0.0133 * V_blood_liver
  V_SP_blood_liver    <- 0.1657 * V_blood_liver
  V_ML_blood_liver    <- 0.00392 * V_blood_liver
  V_water_blood_liver <- V_blood_liver - V_SA_blood_liver - V_Glob_blood_liver - V_SP_blood_liver - V_ML_blood_liver
  V_sorb_blood_liver  <- V_blood_liver - V_water_blood_liver
  
  V_SA_blood_liver_hep    <- 0.0147 * V_blood_liver_hep
  V_Glob_blood_liver_hep  <- 0.0133 * V_blood_liver_hep
  V_SP_blood_liver_hep    <- 0.1657 * V_blood_liver_hep
  V_ML_blood_liver_hep    <- 0.00392 * V_blood_liver_hep
  V_water_blood_liver_hep <- V_blood_liver_hep - V_SA_blood_liver_hep - V_Glob_blood_liver_hep - V_SP_blood_liver_hep - V_ML_blood_liver_hep
  V_sorb_blood_liver_hep  <- V_blood_liver_hep - V_water_blood_liver_hep
  
  V_FABP_liver  <- 0.0025  * V_liver_tissue
  V_SA_liver    <- 0.0012  * V_liver_tissue
  V_SP_liver    <- 0.2075  * V_liver_tissue
  V_ML_liver    <- 0.02046 * V_liver_tissue
  V_water_liver <- V_liver_tissue - V_FABP_liver - V_SA_liver - V_SP_liver - V_ML_liver
  V_sorb_liver  <- V_liver_tissue - V_water_liver
  
  # Bile
  V_SA_bile     <- 0.01 * V_bile
  V_ML_bile     <- 0.03 * V_bile
  V_water_bile  <- V_bile - V_SA_bile - V_ML_bile
  V_sorb_bile   <- V_SA_bile + V_ML_bile
  Q_bile        <- 1e-7
  
  # Gut
  V_gut_tissue <- 1.71
  V_blood_gut  <- 0.03
  V_gut_lumen  <- 3.15
  A_gut        <- 1200
  A_gut_lumen  <- 4800
  V_gut        <- V_gut_tissue + V_blood_gut
  
  V_SA_blood_gut    <- 0.0147 * V_blood_gut
  V_Glob_blood_gut  <- 0.0133 * V_blood_gut
  V_SP_blood_gut    <- 0.1657 * V_blood_gut
  V_ML_blood_gut    <- 0.00392 * V_blood_gut
  V_water_blood_gut <- V_blood_gut - V_SA_blood_gut - V_Glob_blood_gut - V_SP_blood_gut - V_ML_blood_gut
  V_sorb_blood_gut  <- V_blood_gut - V_water_blood_gut
  
  V_FABP_gut   <- 0.0002
  V_SA_gut     <- 0.0007  * V_gut_tissue
  V_SP_gut     <- 0.0563  * V_gut_tissue
  V_ML_gut     <- 0.0126  * V_gut_tissue
  V_water_gut  <- V_gut_tissue - V_FABP_gut - V_SA_gut - V_SP_gut - V_ML_gut
  V_sorb_gut   <- V_gut_tissue - V_water_gut
  
  # Kidneys
  V_kidneys_tissue <- 0.51
  V_blood_kidneys  <- 0.12
  V_kidneys_lumen  <- 0.0765
  VF_water_filtrate <- 0.81
  A_kidneys        <- 800
  V_kidneys        <- V_kidneys_tissue + V_blood_kidneys
  
  V_SA_blood_kidneys    <- 0.0147 * V_blood_kidneys
  V_Glob_blood_kidneys  <- 0.0133 * V_blood_kidneys
  V_SP_blood_kidneys    <- 0.1657 * V_blood_kidneys
  V_ML_blood_kidneys    <- 0.00392 * V_blood_kidneys
  V_water_blood_kidneys <- V_blood_kidneys - V_SA_blood_kidneys - V_Glob_blood_kidneys - V_SP_blood_kidneys - V_ML_blood_kidneys
  V_sorb_blood_kidneys  <- V_blood_kidneys - V_water_blood_kidneys
  
  V_SA_kidneys    <- 0.0014 * V_kidneys_tissue
  V_FABP_kidneys  <- 0.0031 * V_kidneys_tissue
  V_SP_kidneys    <- 0.1767 * V_kidneys_tissue
  V_ML_kidneys    <- 0.0206 * V_kidneys_tissue
  V_water_kidneys <- V_kidneys_tissue - V_SA_kidneys - V_FABP_kidneys - V_SP_kidneys - V_ML_kidneys
  V_sorb_kidneys  <- V_kidneys_tissue - V_water_kidneys
  
  # Rest
  V_rest_tissue <- 24.13
  V_blood_rest  <- 0.33
  A_rest        <- 400
  V_rest        <- V_rest_tissue + V_blood_rest
  
  V_SA_blood_rest    <- 0.0147 * V_blood_rest
  V_Glob_blood_rest  <- 0.0133 * V_blood_rest
  V_SP_blood_rest    <- 0.1657 * V_blood_rest
  V_ML_blood_rest    <- 0.00392 * V_blood_rest
  V_water_blood_rest <- V_blood_rest - V_SA_blood_rest - V_Glob_blood_rest - V_SP_blood_rest - V_ML_blood_rest
  V_sorb_blood_rest  <- V_blood_rest - V_water_blood_rest
  
  V_SA_rest    <- 0.00086
  V_FABP_rest  <- 0
  V_SP_rest    <- 0.1128 * V_rest_tissue
  V_ML_rest    <- 0.0114 * V_rest_tissue
  V_water_rest <- V_rest_tissue - V_FABP_rest - V_SP_rest - V_ML_rest
  V_sorb_rest  <- V_rest_tissue - V_water_rest
  
  # Central blood
  V_blood      <- 1.0068
  V_SA_blood   <- 0.0147 * V_blood
  V_Glob_blood <- 0.0133 * V_blood
  V_SP_blood   <- 0.1657 * V_blood
  V_ML_blood   <- 0.00392 * V_blood
  V_water_blood <- V_blood - V_SA_blood - V_Glob_blood - V_SP_blood - V_ML_blood
  V_sorb_blood  <- V_blood - V_water_blood
  
  # Volume fractions in plasma (unitless, L_constituent / L_plasma)
  VF_water_plasma <- 0.92800
  VF_SA_plasma    <- 0.02940
  VF_Glob_plasma  <- 0.02000
  VF_SP_plasma    <- 0.02210
  VF_ML_plasma    <- 0.00170
  # Note: VF_plasma (plasma fraction of whole blood) ~ 0.60 but cancels in the final formula
  
  # Excretion
  Q_feces   <- 3.48/86400
  Q_urine   <- 2.26/86400
  A_GF      <- 35
  A_tubular <- 11.6
  
  # Partition coefficients
  K_blood   <- V_water_blood / V_blood + V_SA_blood / V_blood * K_SA + V_Glob_blood / V_blood * K_Glob + V_SP_blood / V_blood * K_SP + V_ML_blood / V_blood * K_ML
  K_liver   <- V_water_liver / V_liver_tissue + V_FABP_liver / V_liver_tissue * K_FABP + V_SA_liver / V_liver_tissue * K_SA + V_SP_liver / V_liver_tissue * K_SP + V_ML_liver / V_liver_tissue * K_ML
  K_bile    <- V_water_bile / V_bile + V_SA_bile / V_bile * K_SA + V_ML_bile / V_bile * K_ML
  K_bile_liver <- K_bile / K_liver
  K_gut     <- V_water_gut / V_gut_tissue + V_FABP_gut / V_gut_tissue * K_FABP + V_SP_gut / V_gut_tissue * K_SP + V_ML_gut / V_gut_tissue * K_ML + V_SA_gut / V_gut_tissue * K_SA
  K_kidneys <- V_water_kidneys / V_kidneys_tissue + V_SA_kidneys / V_kidneys_tissue * K_SA + V_FABP_kidneys / V_kidneys_tissue * K_FABP + V_SP_kidneys / V_kidneys_tissue * K_SP + V_ML_kidneys / V_kidneys_tissue * K_ML
  K_rest    <- V_water_rest / V_rest_tissue + V_FABP_rest / V_rest_tissue * K_FABP + V_SP_rest / V_rest_tissue * K_SP + V_ML_rest / V_rest_tissue * K_ML + V_SA_rest / V_rest_tissue * K_SA
  K_plasma  <- VF_water_plasma + VF_SA_plasma * K_SA + VF_Glob_plasma * K_Glob + VF_SP_plasma * K_SP + VF_ML_plasma * K_ML
  
  # Free fractions
  f_free_blood    <- 1/(1 + K_blood    * V_sorb_blood    / V_water_blood)
  f_free_liver    <- 1/(1 + K_liver    * V_sorb_liver    / V_water_liver)
  f_free_gut      <- 1/(1 + K_gut      * V_sorb_gut      / V_water_gut)
  f_free_kidneys  <- 1/(1 + K_kidneys  * V_sorb_kidneys  / V_water_kidneys)
  f_free_rest     <- 1/(1 + K_rest     * V_sorb_rest     / V_water_rest)
  f_free_filtrate <- f_free_blood
  
  # Desorption/absorption rates
  k_des_blood   <- 0.118
  k_des_liver   <- 0.118
  k_des_bile    <- 0.118
  k_ab_bile     <- (k_des_bile * V_liver_tissue) / (V_bile * K_bile_liver)
  k_des_gut     <- 0.118
  k_des_kidneys <- 0.118
  k_des_rest    <- 0.118
  
  # --- Initial conditions (IV bolus = 0; oral via events/constant input below) ---
  total_IV_dose <- 0
  C_blood_t0    <- total_IV_dose / V_blood
  C_free_blood_t0  <- C_blood_t0 * f_free_blood * V_blood / V_water_blood
  C_bound_blood_t0 <- C_blood_t0 * (1 - f_free_blood) * V_blood / V_sorb_blood
  C_gut_lumen_t0 <- 0
  
  yini <- c(
    C_free_blood = C_free_blood_t0,
    C_bound_blood = C_bound_blood_t0,
    C_free_blood_liver = C_free_blood_t0,
    C_bound_blood_liver = C_bound_blood_t0,
    C_free_blood_liver_hep = 0,
    C_bound_blood_liver_hep = 0,
    C_free_blood_gut = C_free_blood_t0,
    C_bound_blood_gut = C_bound_blood_t0,
    C_free_blood_kidneys = C_free_blood_t0,
    C_bound_blood_kidneys = C_bound_blood_t0,
    C_free_blood_rest = C_free_blood_t0,
    C_bound_blood_rest = C_bound_blood_t0,
    C_free_liver = 0, C_bound_liver = 0,
    C_bile = 0,
    C_free_gut = 0, C_bound_gut = 0,
    C_free_kidneys = 0, C_bound_kidneys = 0,
    C_filtrate = 0,
    C_free_rest = 0, C_bound_rest = 0,
    C_gut_lumen = C_gut_lumen_t0,
    J_blood_to_filtrate = 0, J_liver_to_bile = 0,
    Excreted_feces = 0, Excreted_urine = 0,
    Absorbed_total = 0
  )
  
  # --- Dosing controls ---
  # Continuous (mg/s) inputs derived from mg/kg/day parameters
  amount_in_blood <- constant_iv_mg_per_kg_day * bw_kg / 86400    # goes to central blood
  amount_in_gut   <- constant_oral_mg_per_kg_day * bw_kg / 86400  # goes to gut lumen
  
  # Event-based oral dosing (repeated &/or single bolus)
  event_rows <- list()
  
  # Repeated oral dose events, if requested
  if (!isTRUE(all.equal(dose_mg_per_kg, 0))) {
    if (interval_hours <= 0) stop("interval_hours must be > 0 when dose_mg_per_kg > 0")
    dose_mg_each <- dose_mg_per_kg * bw_kg
    dose_times_h <- seq(0, exposure_duration_days * 24, by = interval_hours)
    dose_times_s <- unique(pmax(0, pmin(exposure_duration_days * 24, dose_times_h))) * 3600
    event_rows[[length(event_rows) + 1]] <- data.frame(
      var = "C_gut_lumen",
      time = dose_times_s,
      value = dose_mg_each / V_gut_lumen,
      method = "add"
    )
  }
  
  # Single oral bolus at t = 0, if requested
  if (!isTRUE(all.equal(single_oral_dose_mg_per_kg, 0))) {
    bolus_mg <- single_oral_dose_mg_per_kg * bw_kg
    event_rows[[length(event_rows) + 1]] <- data.frame(
      var = "C_gut_lumen",
      time = 0,
      value = bolus_mg / V_gut_lumen,
      method = "add"
    )
  }
  
  dose_events <- if (length(event_rows)) do.call(rbind, event_rows) else NULL
  events_list <- if (!is.null(dose_events) && nrow(dose_events) > 0) list(data = dose_events) else NULL
  
  # --- ODE system ---
  rigidode <- function(t, y, parms) {
    with(as.list(y), {
      # Incoming fluxes (mg/s) from continuous inputs
      J_in_blood_free  <- amount_in_blood * f_free_blood
      J_in_blood_bound <- amount_in_blood * (1 - f_free_blood)
      J_in_gut_lumen   <- amount_in_gut
      
      # Blood – central
      J_bound_free_blood <- k_des_blood * V_sorb_blood * (C_bound_blood - C_free_blood * K_blood)
      
      # Liver blood
      J_blood_liver_free_in  <- Q_blood_liver * V_water_blood_liver / V_blood_liver * (C_free_blood - C_free_blood_liver)
      J_blood_liver_bound_in <- Q_blood_liver * V_sorb_blood_liver  / V_blood_liver * (C_bound_blood - C_bound_blood_liver)
      J_bound_free_blood_liver <- k_des_blood * V_sorb_blood_liver * (C_bound_blood_liver - C_free_blood_liver * K_blood)
      
      # Hepatic portal (from gut)
      J_bound_free_blood_liver_hep <- k_des_blood * V_sorb_blood_liver_hep * (C_bound_blood_liver_hep - C_free_blood_liver_hep * K_blood)
      J_hep_vein_to_liver_Papp  <- P_app * A_liver_hep * (C_free_blood_liver_hep - C_free_liver)
      J_hep_vein_to_liver_OATPs <- (P_OATP1B2 + P_OATP2B1 + P_NTCP) * A_liver_hep * C_free_blood_liver_hep
      
      # Gut blood
      J_blood_gut_free_in  <- Q_blood_gut * V_water_blood_gut / V_blood_gut * (C_free_blood - C_free_blood_gut)
      J_blood_gut_bound_in <- Q_blood_gut * V_sorb_blood_gut  / V_blood_gut * (C_bound_blood - C_bound_blood_gut)
      J_bound_free_blood_gut <- k_des_blood * V_sorb_blood_gut * (C_bound_blood_gut - C_free_blood_gut * K_blood)
      
      # Kidneys blood
      J_blood_kidneys_free_in  <- Q_blood_kidneys * V_water_blood_kidneys / V_blood_kidneys * (C_free_blood - C_free_blood_kidneys)
      J_blood_kidneys_bound_in <- Q_blood_kidneys * V_sorb_blood_kidneys  / V_blood_kidneys * (C_bound_blood - C_bound_blood_kidneys)
      J_bound_free_blood_kidneys <- k_des_blood * V_sorb_blood_kidneys * (C_bound_blood_kidneys - C_free_blood_kidneys * K_blood)
      
      # Rest blood
      J_blood_rest_free_in  <- Q_blood_rest * V_water_blood_rest / V_blood_rest * (C_free_blood - C_free_blood_rest)
      J_blood_rest_bound_in <- Q_blood_rest * V_sorb_blood_rest  / V_blood_rest * (C_bound_blood - C_bound_blood_rest)
      J_bound_free_blood_rest <- k_des_blood * V_sorb_blood_rest * (C_bound_blood_rest - C_free_blood_rest * K_blood)
      
      # Liver tissue
      J_blood_liver_Papp  <- P_app * A_liver * (C_free_blood_liver - C_free_liver)
      J_blood_liver_OATPs <- (P_OATP1B2 + P_OATP2B1 + P_NTCP) * A_liver * C_free_blood_liver
      J_bound_free_liver  <- k_des_liver * V_sorb_liver * (C_bound_liver - C_free_liver * K_liver)
      J_liver_to_bile     <- k_ab_bile * V_water_bile * (C_free_liver - C_bile / K_bile)
      J_bile_to_gut       <- Q_bile * C_bile
      
      # Gut tissue / portal
      J_gut_to_hep_blood_Papp <- P_app * A_gut * (C_free_gut - C_free_blood_liver_hep)
      J_gut_to_hep_blood_AT   <- (P_OATP2B1) * A_gut * C_free_gut
      J_bound_free_gut        <- k_des_gut * V_sorb_gut * (C_bound_gut - C_free_gut * K_gut)
      
      # Gut lumen
      J_gut_lumen_to_gut <- P_app * A_gut_lumen * (C_gut_lumen - C_free_gut)
      J_gut_to_feces     <- Q_feces * C_gut_lumen
      
      # Kidneys tissue & filtrate
      J_blood_kidneys_Papp <- P_app * A_kidneys * (C_free_blood_kidneys - C_free_kidneys)
      J_blood_kidneys_OAT  <- (P_OAT1 + P_OAT2 + P_OAT3) * A_kidneys * C_free_blood_kidneys
      J_kidneys_filtrate   <- P_app * A_tubular * (C_filtrate * f_free_filtrate * 1/VF_water_filtrate - C_free_kidneys)
      J_bound_free_kidneys <- k_des_kidneys * V_sorb_kidneys * (C_bound_kidneys - C_free_kidneys * K_kidneys)
      J_GlomFil            <- P_app * A_GF * C_free_blood_kidneys
      J_filtrate_to_urine  <- Q_urine * C_filtrate
      
      # Rest tissue
      J_blood_rest_Papp <- P_app * A_rest * (C_free_blood_rest - C_free_rest)
      J_bound_free_rest <- k_des_rest * V_sorb_rest * (C_bound_rest - C_free_rest * K_rest)
      
      # Mass balances
      dC_free_blood  <- (J_in_blood_free + J_bound_free_blood - J_blood_liver_free_in - J_blood_gut_free_in - J_blood_kidneys_free_in - J_blood_rest_free_in) / V_water_blood
      dC_bound_blood <- (J_in_blood_bound - J_bound_free_blood - J_blood_liver_bound_in - J_blood_gut_bound_in - J_blood_kidneys_bound_in - J_blood_rest_bound_in) / V_sorb_blood
      
      dC_free_blood_liver  <- (J_blood_liver_free_in + J_bound_free_blood_liver - J_blood_liver_Papp - J_blood_liver_OATPs) / V_water_blood_liver
      dC_bound_blood_liver <- (J_blood_liver_bound_in - J_bound_free_blood_liver) / V_sorb_blood_liver
      
      dC_free_blood_liver_hep  <- (J_gut_to_hep_blood_AT + J_gut_to_hep_blood_Papp + J_bound_free_blood_liver_hep - J_hep_vein_to_liver_Papp - J_hep_vein_to_liver_OATPs) / V_water_blood_liver_hep
      dC_bound_blood_liver_hep <- (- J_bound_free_blood_liver_hep) / V_sorb_blood_liver_hep
      
      dC_free_liver  <- (J_blood_liver_Papp + J_blood_liver_OATPs + J_bound_free_liver - J_liver_to_bile + J_hep_vein_to_liver_Papp + J_hep_vein_to_liver_OATPs) / V_water_liver
      dC_bound_liver <- (- J_bound_free_liver) / V_sorb_liver
      
      dC_bile <- (J_liver_to_bile - J_bile_to_gut) / V_bile
      
      dC_free_blood_gut  <- (J_blood_gut_free_in + J_bound_free_blood_gut) / V_water_blood_gut
      dC_bound_blood_gut <- (J_blood_gut_bound_in - J_bound_free_blood_gut) / V_sorb_blood_gut
      
      dC_free_gut  <- (J_bound_free_gut - J_gut_to_hep_blood_Papp - J_gut_to_hep_blood_AT + J_gut_lumen_to_gut) / V_water_gut
      dC_bound_gut <- (- J_bound_free_gut) / V_sorb_gut
      
      dC_gut_lumen <- (J_in_gut_lumen + J_bile_to_gut - J_gut_lumen_to_gut - J_gut_to_feces) / V_gut_lumen
      
      dC_free_blood_kidneys  <- (J_blood_kidneys_free_in + J_bound_free_blood_kidneys - J_blood_kidneys_Papp - J_blood_kidneys_OAT - J_GlomFil) / V_water_blood_kidneys
      dC_bound_blood_kidneys <- (J_blood_kidneys_bound_in - J_bound_free_blood_kidneys) / V_sorb_blood_kidneys
      
      dC_free_kidneys  <- (J_blood_kidneys_Papp + J_bound_free_kidneys + J_blood_kidneys_OAT + J_kidneys_filtrate) / V_water_kidneys
      dC_bound_kidneys <- (- J_bound_free_kidneys) / V_sorb_kidneys
      
      dC_filtrate <- (- J_kidneys_filtrate + J_GlomFil - J_filtrate_to_urine) / V_kidneys_lumen
      
      dC_free_blood_rest  <- (J_blood_rest_free_in + J_bound_free_blood_rest - J_blood_rest_Papp) / V_water_blood_rest
      dC_bound_blood_rest <- (J_blood_rest_bound_in - J_bound_free_blood_rest) / V_sorb_blood_rest
      
      dC_free_rest  <- (J_blood_rest_Papp + J_bound_free_rest) / V_water_rest
      dC_bound_rest <- (- J_bound_free_rest) / V_sorb_rest
      
      dExcreted_feces <- J_gut_to_feces
      dExcreted_urine <- J_filtrate_to_urine
      dAbsorbed_total <- J_gut_to_hep_blood_Papp + J_gut_to_hep_blood_AT
      
      list(c(
        dC_free_blood, dC_bound_blood,
        dC_free_blood_liver, dC_bound_blood_liver,
        dC_free_blood_liver_hep, dC_bound_blood_liver_hep,
        dC_free_blood_gut, dC_bound_blood_gut,
        dC_free_blood_kidneys, dC_bound_blood_kidneys,
        dC_free_blood_rest, dC_bound_blood_rest,
        dC_free_liver, dC_bound_liver,
        dC_bile,
        dC_free_gut, dC_bound_gut,
        dC_free_kidneys, dC_bound_kidneys,
        dC_filtrate,
        dC_free_rest, dC_bound_rest,
        dC_gut_lumen,
        J_blood_to_filtrate = 0,  # placeholder
        J_liver_to_bile,
        dExcreted_feces, dExcreted_urine, dAbsorbed_total
      ))
    })
  }
  
  # --- Solve ---
  out <- deSolve::ode(
    y = yini,
    times = times_sec,
    func = rigidode,
    parms = NULL,
    method = "lsoda",
    events = events_list
  )
  res <- as.data.frame(out, check.names = FALSE)
  
  # --- Compute total concentration per major compartment ---
  total_conc <- function(C_free, C_bound, V_w, V_s, V_tot) (C_free * V_w + C_bound * V_s) / V_tot
  
  res$time_h <- res$time / 3600
  # Blood (central only)
  C_blood <- total_conc(res$C_free_blood, res$C_bound_blood, V_water_blood, V_sorb_blood, V_blood)
  C_plasma <- C_blood * (K_plasma / K_blood)  # convert blood to plasma concentration via composition-based conversion
  # Liver (tissue only)
  C_liver <- total_conc(res$C_free_liver, res$C_bound_liver, V_water_liver, V_sorb_liver, V_liver)
  # Kidneys (tissue)
  C_kidneys <- total_conc(res$C_free_kidneys, res$C_bound_kidneys, V_water_kidneys, V_sorb_kidneys, V_kidneys)
  # Gut (tissue)
  C_gut <- total_conc(res$C_free_gut, res$C_bound_gut, V_water_gut, V_sorb_gut, V_gut)
  # Rest (tissue)
  C_rest <- total_conc(res$C_free_rest, res$C_bound_rest, V_water_rest, V_sorb_rest, V_rest)
  
  out_df <- dplyr::tibble(
    time_h = res$time_h,
    C_blood = C_blood,
    C_plasma = C_plasma,
    C_liver = C_liver,
    C_kidneys = C_kidneys,
    C_gut = C_gut,
    C_rest = C_rest
  )
  
  out_df
}

########## EXAMPLE USAGE ##############
# 1) Single oral bolus, no repeats
#Fischer_PBPK_mouse("PFBS", dose_mg_per_kg = 0, exposure_duration_days = 30, interval_hours = 24,
#                  single_oral_dose_mg_per_kg = 10)
# 
# # 2) Constant oral intake (5 mg/kg/day), no bolus/repeats
# constant <- Fischer_PBPK_mouse("PFBS", dose_mg_per_kg = 5, exposure_duration_days = 30, interval_hours = 24,
#                                constant_oral_mg_per_kg_day = 5)
# # plot to test
# constant_plot <- constant %>% 
#   pivot_longer(cols = contains("C_"), names_to = "Compartment", values_to = "Concentration (mg/L)") %>% 
#   ggplot(aes(x = time_h, y = `Concentration (mg/L)`, color = Compartment)) +
#   geom_line() +
#   scale_y_log10() +
#   theme_minimal(base_size = 15)
# 
# plotly::ggplotly(constant_plot)
 
# # 3) Constant IV intake (1 mg/kg/day) + daily oral (2 mg/kg per event)
# Fischer_PBPK_mouse("PFBS", dose_mg_per_kg = 2, exposure_duration_days = 14, interval_hours = 24,
#                    constant_iv_mg_per_kg_day = 1)

# # 4) Repeat oral dose (daily gavage) (5 mg/kg/day)
# repeat_dose <- Fischer_PBPK_mouse("PFBS", dose_mg_per_kg = 5, exposure_duration_days = 30, interval_hours = 24,
#                                constant_oral_mg_per_kg_day = 0)
# # plot to test
# repeat_plot <- repeat_dose %>% 
#   pivot_longer(cols = contains("C_"), names_to = "Compartment", values_to = "Concentration (mg/L)") %>% 
#   ggplot(aes(x = time_h, y = `Concentration (mg/L)`, color = Compartment)) +
#   geom_line() +
#   scale_y_log10() +
#   theme_minimal(base_size = 15)
# 
# plotly::ggplotly(repeat_plot)
