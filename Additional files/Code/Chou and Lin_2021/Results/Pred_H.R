## Define the prediction function based on study design of literature
## Exposure scenario:  Dosing druing pregnnacy from pre-pregnant to gestational stage

Pred.preG <- function (DOSE) {
    
    ## Defined the exposure scenario for age specific data
    ex <- tibble(ID   = rep(1, 365*30 + 1), # individual ID 
                 time = seq(from = 0, to = 24*365*30, by = 24)) %>% # time from 0 to 30 years
        mutate (DAY   = time/24) %>%  # DAY           
        mutate (YEAR  = DAY/365) %>%  # AGE
        mutate (BW    = if_else(YEAR <= 18, # if age small than or equal to 18 years old use the equation; otherwise bodyweight equal to about 54 kg
                                true  = (-2.561*YEAR^4 + 85.576*YEAR^3 - 855.95*YEAR^2 + 5360.6*YEAR + 4428.5)/1000,
                                ifelse(YEAR >= 30, 67, 54)))
    
    tsamp <- tgrid (start = 0, end = 365*30, delta = 0.1) # simulation time from 0 to 30 years old, and the time interval is 0.1
    
    ## Exposure scenario: A dynamic exposure model for PFOS daily intakes 
    ## a constant exposure (per kg of bodyweight) from birth until the year 2000 (Assumed 20 years old), then 66% decrease for PFOS  (Paul et al., 2009) 
    PDOSEoral_1 = DOSE ## assumed the daily oral dose before 2000
    ex_1 <- mutate (ex, amt = case_when (YEAR < 20 ~ PDOSEoral_1*BW,
                                         YEAR >= 20 ~ PDOSEoral_1*BW), 
                    cmt  = "AST", ii = 24, evid = 1, time = DAY)
    
    
    out <- PreGmod_H %>%
        update(atol = 1E-3, maxsteps = 50000) %>%          
        mrgsim_d(data = ex_1, tgrid = tsamp)%>%
        filter (time > 0)
    
    return (out)
}

## Simulation the initial concentration of gestational model based on different exposure scenario
## The assumed epxosure dose were estiamted from previous literatures
## Init_A1: Japan population;   Dose: 1.2 ng/kg/day estiamted from Table 5 in Loccisano et al. (2013);
## Init_A2: Danish population;  Dose: 3.7 ng/kg/day estiamted from Table 5 in Loccisano et al. (2013);  
## Init_A3: U.S. population     Dose: 0.1- 2 ng/kg/day estimated from Loccisano et al. (2014)
## Init_A4: Chinese population  Dose: 1.1 ng/kg/day estimaed from Zhang et al. (2010)
## Init_A5: Sweden population   Dose: 0.860 - 1.44 ng/kg/day estimated from Robin Vestergren et al. (2012)
## Init_A6: Germany population  Dose: 1.35 ng/kg/day estimated from Table 5 in Loccisano et al. (2013);

DOSE_A1 <- 1.2e-6
DOSE_A2 <- 3.7e-6
DOSE_A3 <- 1.5e-6
DOSE_A4 <- 1.1e-6
DOSE_A5 <- 0.86e-6
DOSE_A6 <- 1.35e-6

Init_A1 <- Pred.preG (DOSE = DOSE_A1) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A2 <- Pred.preG (DOSE = DOSE_A2) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A3 <- Pred.preG (DOSE = DOSE_A3) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A4 <- Pred.preG (DOSE = DOSE_A4) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A5 <- Pred.preG (DOSE = DOSE_A5) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A6 <- Pred.preG (DOSE = DOSE_A6) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))

## Exposure sceinario 
pred.G <- function(Gpars, DOSE, Init) {
    
    ## Get out of log domain
    Gpars <- lapply(Gpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 60                     ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; Repeat oral dose from GA0 - GA40 (weeks)
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time from GD0 to GD 40
    
    ## Simulation of exposure scenaior (Repeated oral dose to 1/2/3/5/10 mg/kg)
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>%  ## Input the intial concentration
        param (Gpars) %>%                 ## Update the parameter list with Gpars
        update(atol = 1E-6,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    Goutdf = cbind.data.frame(Time     = Gout$time, 
                              CPlas    = Gout$Plasma*1000, 
                              CL       = Gout$Liver*1000,
                              CFPlas   = Gout$CordB*1000,
                              CPla     = Gout$Placenta*1000,
                              CFL      = Gout$Fliver*1000)
    
    return (Goutdf) # Return Goutdf
}

## Prediction function from pre-pregnant to gestational exposure
pred.G2 <- function(DOSE) {
    
    Init <- Pred.preG (DOSE) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ## Get out of log domain
    Gpars <- lapply(GFit_H$par, exp)        ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 67                     ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; Repeat oral dose from GD2 - GD21
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSE every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 0,             ## Dossed strat time (GD0)
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD0 to GD 40 (weeks) + 1 day
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>% ## Input the intial concentrations 
        param (Gpars) %>% ## Update the parameter list with Gpars
        update(atol = 1E-3,  maxsteps = 50000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
    
    Init_G <- Gout %>% filter (time == 40*24*7) ## Extract the concentration at GD40
    return (Init_G)
}   

## Simulation the initial concentration of gestational model based on different exposure scenario
## The assumed epxosure dose were estiamted from previous literatures
## DOSE_B1: Sweden population (Karrman et al. (2007);                Dose: 2.3 ng/kg/day estimated from Table 5 in Loccisano et al. (2013);   
## DOSE_B2: U.S. population exposure (Von Ehrenstein et al., 2009);  Dose: 2.25 ng/kg/day estimated from Table 5 in Loccisano et al. (2013); 
## DOSE_B3: Gemerny population exposure (Fromme et al., 2010);       Dose: 0.38 ng/kg/day estimated from Table 5 in Loccisano et al. (2013); 
## DOSE_B4: Sourth Korea population exposure (Lee et al., 2018);     Dose: 0.35 ng/kg/day estimated from Table 5 in Loccisano et al. (2013); 

DOSE_B1 <- 2.3E-6
DOSE_B2 <- 2.25E-6
DOSE_B3 <- 0.38E-6
DOSE_B4 <- 0.35E-6

Init_G1 <- pred.G2 (DOSE_B1)
Init_G2 <- pred.G2 (DOSE_B2)
Init_G3 <- pred.G2 (DOSE_B3)
Init_G4 <- pred.G2 (DOSE_B4)

## Prediction function for lactational
pred.L <- function(Lpars, DOSE, Init_G) {
    
    ## Get out of log domain
    Lpars <- lapply(Lpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 67                     ## Rat body weight on PND0 
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 41*7                   ## Total dosing/Dose times; Repeat oral dose from PND0 - PND41 (weeks)
    LDOSE        = DOSE                   ## Repeat oral dose  (mg/kg/day);  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose (mg)
    
    
    # To create a data set of 1 subject receiving DOSE every 24 hours
    Lex.oral <- ev (ID   = 1,             ## One individual
                    amt  = LDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = LTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 24) ## Simulation time from PND0 to PND74 with dosing of LTDOSE - 1+ 2 days (No dosing)
    
    
    Lout <- Lmod_H %>% # Lactational model
        init(APlas_free =  Init_G$APlas_free, APTC =  Init_G$APTC, AFil =  Init_G$AFil, AKb =  Init_G$AKb, ARest =  Init_G$ARest,
             AL =  Init_G$AL, AM =  Init_G$AM, AF =  Init_G$AF, A_baso =  Init_G$A_baso, A_apical = Init_G$A_apical, Adif =  Init_G$Adif,
             Aefflux =  Init_G$Aefflux, APlas_free_neo =  Init_G$APlas_Fet_free, 
             AL_neo =  Init_G$AL_Fet, ARest_neo =  Init_G$ARest_Fet) %>% # Input the intial concentration
        param (Lpars) %>% # Update the parameter list with pars
        update(atol = 1E-6, maxsteps = 500000) %>% # Atol: absolute tolerance parameter; hmax: The maximum step size; maxstep: maximum number of steps the solver will take when advancing from one time to the next         
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    ## Extract the concentration 
    Loutdf = cbind.data.frame(Time   = Lout$time + 40*7*24, 
                              CPlas  = Lout$Plasma*1000, 
                              CPlas_pup = Lout$CPneo*1000,
                              CMilk = Lout$Milk*1000)
    
    outdf  <- rbind.data.frame (Loutdf)
    return (Loutdf)
    
}