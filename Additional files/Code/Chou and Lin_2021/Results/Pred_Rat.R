## Define the prediction function based on study design of literature of Thibodeaux et al., 2003
## Exposure scenario: Dosing druing pregnnacy from GD 2 - GD20 
pred.A <- function(Gpars, DOSE) { ## Gpars: input parameters, Dose: input dose, Dose regimen: 1, 2, 3, 5, 10 mg/kg/day
    
    ## Get out of log domain
    Gpars <- exp(Gpars)                   ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.225                  ## Body weight during gestation (measrument data if available); Default value adopted from Loccisano et al. (2012); 0.225-0.307 kg
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 20                     ## Total dosing/Dose times; Repeat oral dose from GD2 - GD20
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    ## To create a data set of 1 subject receiving DOSEoral every 24 hours for 1 total doses
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 24*2,          ## Dossed strat time (GD2)
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*2, 1) ## set up the output time; dtout = 1 hours 
    
    ## Simulation of exposure scenaior (Repeated oral dose to 1/2/3/5/10 mg/kg)
    Gout <- 
        Gmod_R %>% ## Gestational PBPK model
        param (Gpars) %>% ## Update the parameter list with Gpars
        update(atol = 1E-6, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
    
    ## Extract the "Time", "CPlas", "CL" and "CFL" from Gout; 
    Goutdf = cbind.data.frame(Time   = Gout$time, 
                              CPlas  = Gout$Plasma, # CPlas: PFOS concentration in maternal plasma
                              CL     = Gout$Liver,  # CL: PFOS concentration in maternal liver
                              CFL    = Gout$Liver_Fet) # CFL: PFOS concentration in fetal liver
    
    Goutdf <- Goutdf %>% filter (Time > 0)
    return (Goutdf) # Return Goutdf
}


## Exposure scenario B. Dosing from GD0 - PND20; the dams were sacrificed on GD 20; PND20 based on study design of Chang et al. (2009) 
pred.B <- function(Lpars, DOSE) {
    
    ## Get out of log domain
    Gpars <- lapply(GFit_R$par, exp)        ## Return a list of exp (parametrs for gestational model) from log scale
    Lpars <- lapply(Lpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.225                  ## Body weight during gestation (measrument data if available); Default value adopted from Loccisano et al. (2012); 0.225-0.307 kg
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from GD0 - GD20
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    ## To create a data set of 1 subject receiving DOSEoral every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## set up the output time; dtout = 1 hours 
    
    ## Simulation of exposure scenaior during gestation (oral repeat dose to 0.1, 0.3, 1 mg/kg)
    Gout <- 
        Gmod_R %>%
        param (Gpars) %>%
        update(atol = 1E-6, maxsteps = 50000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Goutdf = cbind.data.frame(Time      = Gout$time, 
                              CPlas     = Gout$Plasma,     # Maternal plasma
                              CL        = Gout$Liver,      # Maternal liver
                              CPlas_pup = Gout$Plasma_Fet, # Fetal plasma
                              CL_pup    = Gout$Liver_Fet)  # Fetal liver
    
    Init <- Gout %>% filter (time == 21*24) %>% select(-c("AUC_CPlas","AUC_CPlas_Fet", "Plasma", "Liver", "Plasma_Fet", "Liver_Fet"))
    
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.242                  ## Body weight; Default value adopted from Loccisano et al. (2012); 0.242-0.292 kg
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from PND0 - PND21
    LDOSE        = DOSE                   ## Repeat oral dose (0.1, 0.3, 1 mg/kg);  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Lex.oral <- ev (ID   = 1,            ## One individual
                    amt  = LDOSEoral,    ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = LTDOSE - 1,   ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,         ## Infusion time;  
                    replicate = FALSE)   ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 1) ## Simulation time from PND0 to PND21 with dosing of LTDOSE - 1 + 2 days (No dosing)
    
    ## The concentration at the end of gestation used as initial concentration of lactation model 
    Lout <- Lmod_R %>% 
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux, APlas_free_pup = Init$APlas_Fet_free, AL_pup = Init$AL_Fet, ARest_pup = Init$ARest_Fet) %>%
        param (Lpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    Loutdf = cbind.data.frame(Time   = Lout$time + 22*24,   # Simulation time + GD22 (22*24 hours)
                              CPlas  = Lout$Plasma,         # Maternal plasma
                              CL     = Lout$Liver,          # Matnerla liver
                              CPlas_pup = Lout$Plasma_pup,  # Neonatal palsma
                              CL_pup    = Lout$Liver_pup)   # Neonatal liver
    
    return (Loutdf)
}


