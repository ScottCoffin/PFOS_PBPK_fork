##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(ggExtra)     ## R-package for function "grid.arrange" 
library(EnvStats)

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")
source (file = "HMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_R <- mcode ("PreG_RatPBPK.code", PreG_RatPBPK.code)
Gmod_R    <- mcode ("GRatPBPK.code", GRatPBPK.code)
Lmod_R    <- mcode ("LRatPBPK.code", LRatPBPK.code)
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code) 
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Loading datasets
RDat <- read.csv(file = "Data_R.csv") 
HDat <- read.csv(file = "Data_H.csv")

## Loading Fit results from rds files
GFit_R <- readRDS(file = "GFit_R.rds")
LFit_R <- readRDS(file = "LFit_R.rds")
GFit_H <- readRDS(file = "GFit_H.rds")
LFit_H <- readRDS(file = "LFit_H.rds")

## Input mrgsolve-based PBPK Model
source (file = "Pred_Rat.R")
source (file = "Pred_H.R")

## Descrtion of selected animal studies
#==============================================================================================================================
# A1. Lau et al. 2003       : Pregnant SD rats were dosing daily by gavage from gestation day (GD) 2 to GD 21                                                         #
# A2. Luebker et al. 2005a  : Female SD rats were dosing daily beginning 42 days prior to mating and through the mating period (14 days ) to GD20    
# A3. Luebker et al. 2005b  : Female SD rats were dosing daily beginning 42 days prior to mating and through the mating period (14 days ) to PND 20                                                       
# A4. Butenhoff et al. 2009 : Female SD rats dosing daily beginning gestation day (GD) 0 through postnatal day (PND) 20.                                        
#========================================================================================================================= 

## Define prediction function based on the exposure scenario from literatures 
pred.rat <- function(Gpars_R, Lpars_R, DOSE) {
    
    ## Exposure scenario for Lactational exposure; through GD(-53) - GD21 (end of gestation)
    PreG_BW          = 0.25                         ## Rat body weight; 
    tinterval        = 24                           ## Time interval; 
    PreG_TDOSE       = 42 + 14                      ## Total dosing/Dose times; Repeat oral dose from beginning 42 days prior to cohabitation; (14 days matting)
    PreG_DOSE        = DOSE                         ## Repeat oral dose of GDOSE mg/kg/day;  
    PreG_DOSEoral    = PreG_DOSE*PreG_BW            ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    PreG_ex.oral_R <- ev (ID   = 1,                 ## One individual
                          amt  = PreG_DOSEoral,     ## Amount of dose 
                          ii   = tinterval,         ## Time interval
                          addl = PreG_TDOSE - 1,    ## Addtional doseing 
                          cmt  = "AST",             ## The dosing comaprtment: AST Stomach  
                          replicate = FALSE)        ## No replicate
    
    PreG_tsamp_R       = tgrid(0, tinterval*(PreG_TDOSE - 1), 1) ## Simulation time from 6 weeks prior to mating and during mating (maximum of 14 days)
    
    PreG_out_R <- PreGmod_R %>% 
                  update(atol = 1E-8, maxsteps = 500000) %>%          
                  mrgsim_d (data = PreG_ex.oral_R, tgrid = PreG_tsamp_R)
    
    PreG_Init <- PreG_out_R %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver","Kidney"))
    
    ## Get out of log domain
    Gpars <- exp(Gpars_R)           ## Return a list of exp (parametrs for gestational model) from log scale

    
    ## Exposure scenario for gestational exposure
    GBW          = 0.3                    ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 21                     ## Total dosing/Dose times; Repeat oral dose from GD0 - GD20
    GDOSE        = DOSE                   ## Repeat oral dose 
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    

    # To create a data set of 1 subject receiving DOSEoral every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp_A1  = tgrid(2, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD2 to GD 21 
    Gtsamp_A2  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD0 to GD 21 
    Gtsamp_A3  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD0 to GD 21 
    Gtsamp_A4  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD0 to GD 21 
    
    ## Simulation of exposure scenaior 
    Gout_A1 <- Gmod_R %>% param (Gpars) %>%
               update(atol = 1E-5, maxsteps = 500000) %>%          
               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A1)
    
    Gout_A2 <- Gmod_R %>% init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
                               AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
                               Aefflux = PreG_Init$Aefflux) %>%
                               param (Gpars) %>% update(atol = 1E-5, maxsteps = 500000) %>%          
                               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A2)
    
    Gout_A3 <- Gmod_R %>% init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
                               AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
                               Aefflux = PreG_Init$Aefflux) %>%
                               param (Gpars) %>% update(atol = 1E-5, maxsteps = 500000) %>%          
                               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A3)
    
    Gout_A4 <- Gmod_R %>% param (Gpars) %>%
               update(atol = 1E-5, maxsteps = 500000) %>%          
               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A4)
    

    Goutdf_A1 <- as.data.frame(Gout_A1) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Goutdf_A2 <- as.data.frame(Gout_A2) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Ginit_A3  <- as.data.frame(Gout_A3) %>% filter (row_number()== n()) 
    Ginit_A4  <- as.data.frame(Gout_A4) %>% filter (row_number()== n()) 
    
    ############### Lactational model ##############################
    Lpars <- exp(Lpars_R)          ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.35                  ## Rat body weight on PND0 (GD21); 
    tinterval    = 24                    ## Time interval; 
    LTDOSE       = 21                    ## Total dosing/Dose times; Repeat oral dose from PND0 - PND20
    LDOSE        = DOSE                  ## Repeat oral dose   
    LDOSEoral    = LDOSE*LBW             ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Lex.oral <- ev (ID   = 1,            ## One individual
                    amt  = LDOSEoral,    ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = LTDOSE - 1,   ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*1, 1) ## Simulation time from PND0 to PND20 plus one day without dosing
    
    
    Lout_A3 <- Lmod_R %>% 
               init(APlas_free = Ginit_A3$APlas_free, APTC = Ginit_A3$APTC, AFil = Ginit_A3$AFil, AKb = Ginit_A3$AKb, ARest = Ginit_A3$ARest,
               AL = Ginit_A3$AL, AM = Ginit_A3$AM, AF = Ginit_A3$AF, A_baso = Ginit_A3$A_baso, A_apical = Ginit_A3$A_apical, Adif = Ginit_A3$Adif,
               Aefflux = Ginit_A3$Aefflux, APlas_free_pup = Ginit_A3$APlas_Fet) %>%
               param (Lpars) %>%
               update(atol = 1E-3, maxsteps = 500000) %>%          
               mrgsim_d (data = Lex.oral, tgrid = Ltsamp)
    
    Lout_A4 <- Lmod_R %>% 
               init(APlas_free = Ginit_A4$APlas_free, APTC = Ginit_A4$APTC, AFil = Ginit_A4$AFil, AKb = Ginit_A4$AKb, ARest = Ginit_A4$ARest,
               AL = Ginit_A4$AL, AM = Ginit_A4$AM, AF = Ginit_A4$AF, A_baso = Ginit_A4$A_baso, A_apical = Ginit_A4$A_apical, Adif = Ginit_A4$Adif,
               Aefflux = Ginit_A4$Aefflux, APlas_free_pup = Ginit_A4$APlas_Fet) %>%
               param (Lpars) %>%
               update(atol = 1E-3, maxsteps = 500000) %>%          
               mrgsim_d (data = Lex.oral, tgrid = Ltsamp)
    
    
    Loutdf_A3 <- as.data.frame(Lout_A3) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Loutdf_A4 <- as.data.frame(Lout_A4) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Loutdf_A3$Time <- Loutdf_A3$Time + (24*20)
    Loutdf_A4$Time <- Loutdf_A4$Time + (24*20)
    
    return (list("A1" = Goutdf_A1, "A2" = Goutdf_A2, "A3" = Loutdf_A3, "A4" = Loutdf_A4))
    
}


## Define human model
pred.human <- function(Gpars_H, Lpars_H, DOSE) {
    
    ## Exposure scenario for Lactational exposure; 
    ## Defined the exposure scenario for age specific data
    ex <- tibble(ID   = rep(1, 365*30 + 1), # individual ID 
                 time = seq(from = 0, to = 24*365*30, by = 24)) %>% # time from 0 to 30 years
        mutate (DAY   = time/24, YEAR  = DAY/365, 
        BW    = if_else(YEAR <= 18, # if age small than or equal to 18 years old use the equation; otherwise bodyweight equal to about 54 kg
                                true  = (-2.561*YEAR^4 + 85.576*YEAR^3 - 855.95*YEAR^2 + 5360.6*YEAR + 4428.5)/1000,
                                ifelse(YEAR >= 30, 67, 54)))
    
    tsamp <- tgrid (start = 0, end = 365*30, delta = 0.1) # simulation time from 0 to 30 years old, and the time interval is 0.1
    
    ## Exposure scenarior:  
    # ## Exposure scenarior: Daily dose to 4.5 ng/kg/day for 0-25 years (1976 - 2001); daily dose to 4.5*0.34 ng/kg/day for 6 years (2001 - 2007) matrix: plasma (mg/L)
    PDOSEoral_1 = DOSE #0.0000045  ## mg/kg-d, daily oral dose before 2000
    ex_1 <- mutate (ex, amt = case_when (YEAR < 20 ~ PDOSEoral_1*BW,
                                         YEAR >= 20 ~ PDOSEoral_1*BW*1), # adj.factor: default is 0.34
                    cmt  = "AST", ii = 24, evid = 1, time = DAY)
    
    
    
    out <- PreGmod_H %>%
           update(atol = 1E-8, maxsteps = 500000) %>%          
           mrgsim_d(data = ex_1, tgrid = tsamp)%>%
           filter (time > 0)
    
    PreG_Init <- out %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ############## Gestatonal model
    ## Get out of log domain
    Gpars <- exp(Gpars_H)                 ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 67                     ## Body weight before gestation (measrument data if available); 
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; Repeat oral dose from GA40 (Weeks)
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSEoralevery 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 0,             ## Dossed strat time (GD0)
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp_A1  = tgrid(0, tinterval*(GTDOSE - 1) + 24*7, 1) ## Simulation time from GA0 to GA 40 
    Gtsamp_A2  = tgrid(0, tinterval*(GTDOSE - 1) + 24*7, 1) ## Simulation time from GA0 to GA 40 

    ## Simulation of exposure scenaior A 
    
    
    Gout_A1 <- Gmod_H %>% init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
                               AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
                               Aefflux = PreG_Init$Aefflux) %>%
                               param (Gpars) %>% update(atol = 1E-3, maxsteps = 500000) %>%          
                               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A2)
    
    Gout_A2 <- Gmod_H %>% init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
                               AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
                               Aefflux = PreG_Init$Aefflux) %>%
                               param (Gpars) %>% update(atol = 1E-3, maxsteps = 500000) %>%          
                               mrgsim_d (data = Gex.oral, tgrid = Gtsamp_A2)
    
    
    Goutdf_A1 <- as.data.frame(Gout_A1) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Ginit_A2  <- as.data.frame(Gout_A2) %>% filter (row_number()== n()) 

    ############### Lactational model ##############################
    Lpars <- exp(Lpars_H)          ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for Lactational exposure; 
    LBW          = 67                    ## body weight 
    tinterval    = 24                    ## Time interval; 
    LTDOSE       = 7*41                  ## Total dosing/Dose times; 
    LDOSE        = DOSE                  ## Repeat oral dose of GDOSE mg/kg/day;  
    LDOSEoral    = LDOSE*LBW             ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving LDOSEoral every 24 hours
    Lex.oral <- ev (ID   = 1,            ## One individual
                    amt  = LDOSEoral,    ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = LTDOSE - 1,   ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)   ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*7, 1) ## Simulation time 
    
    
    Lout_A2 <- Lmod_H %>% 
               init(APlas_free = Ginit_A2$APlas_free, APTC = Ginit_A2$APTC, AFil = Ginit_A2$AFil, AKb = Ginit_A2$AKb, ARest = Ginit_A2$ARest,
               AL = Ginit_A2$AL, AM = Ginit_A2$AM, AF = Ginit_A2$AF, A_baso = Ginit_A2$A_baso, A_apical = Ginit_A2$A_apical, Adif = Ginit_A2$Adif,
               Aefflux = Ginit_A2$Aefflux, APlas_free_neo = Ginit_A2$APlas_Fet) %>%
               param (Lpars) %>%
               update(atol = 1E-3, maxsteps = 500000) %>%          
               mrgsim_d (data = Lex.oral, tgrid = Ltsamp)
    
    
    
    
    Loutdf_A2 <- as.data.frame(Lout_A2) %>% select (Time = time, CPlas = Plasma, AUC_CPlas = AUC_CPlas)
    Loutdf_A2$Time <- Loutdf_A2$Time + (24*7*40)

    return (list("G" = Goutdf_A1, "L" = Loutdf_A2))
}

###########
## Defined the Monte Carlo simulation function for rat to estiamte the probablistic AUC based on NOAEL
# A1 Lau et al. 2003       : (GD2- GD21)    : NOAEL = 1.0 mg/kg-day
# A2 Luebker et al. 2005b  : (PREG - G20)   : NOAEL = 0.4 mg/kg-day
# A3 Luebker et al. 2005a  : (PREG - PND20) : NOAEL = 0.1 mg/kg-day
# A4 Butenhoff et al. 2009 : (GD0 - PND20)  : NOAEL = 0.3 mg/kg-day

## Create the list object
HED_A1_A <- list()
HED_A2_A <- list()
HED_A3_A <- list()
HED_A4_A <- list()
HED_A1_B <- list()
HED_A2_B <- list()
HED_A3_B <- list()
HED_A4_B <- list()
ASC_A1   <- list()
ASC_A2   <- list()
ASC_A3   <- list()
ASC_A4   <- list()

N = 50 # For testing purpose, we use 50 iterations. 
## Use N = 1000 to reproduce results in Table 5. The results using 50 vs. 1000 are similar.

for (i in 1:N){
#Defined the parameters with Monte carlo analysis (please refer to Tables S9-10)
    Gpars_R  <- log(c(
        BW        = rnormTrunc (1, min = 0.15, max = 0.57, mean = 0.36, sd = 0.108),
        VLC       = rnormTrunc  (1, min = 0.01, max = 0.06, mean = 0.04, sd = 0.01),
        KbileC    = rlnormTrunc (1, min = exp(-5.38), max = exp(-4.63), meanlog = -5.00, sdlog = 0.19),
        PL        = rlnormTrunc (1, min = exp(0.88), max = exp(1.39), meanlog = 1.13, sdlog = 0.13),
        PRest     = rlnormTrunc (1, min = exp(-1.79), max = exp(-1.28), meanlog = -1.53, sdlog = 0.13),
        Ktrans1C  = rlnormTrunc (1, min = exp(-0.18), max = exp(0.58), meanlog = 0.20, sdlog = 0.19),
        Ktrans2C  = rlnormTrunc (1, min = exp(-0.42), max = exp(0.34), meanlog = -0.04, sdlog = 0.19),
        Ktrans3C  = 0.23,
        Free      = 0.019,
        Vmax_baso_invitro = 221,
        Km_baso   = 19.9,
        Free_Fet  = 0.022,
        PL_Fet    = 2.55
    ))    
    
    Lpars_R  <- log(c(
        Free        = rlnormTrunc (1, min = exp(-6.01), max = exp(-5.32), meanlog = -5.62, sdlog = 0.198),
        Free_pup    = rlnormTrunc (1, min = exp(-4.22), max = exp(-3.44), meanlog = -3.84, sdlog = 0.198),
        GFRC        = rnormTrunc  (1, min = 24.95, max = 57.13, mean = 41, sd = 8.208),
        KbileC      = rlnormTrunc (1, min = exp(-6.57), max = exp(-5.41), meanlog = -5.99, sdlog = 0.294),
        Km_apical_p = rlnormTrunc (1, min = exp(5.009), max = exp(6.16), meanlog = 5.58, sdlog = 0.29),
        KMilk0      = rlnormTrunc (1, min = exp(-1.68), max = exp(-0.90), meanlog = -1.29, sdlog = 0.13), 
        KurineC_pup = rlnormTrunc (1, min = exp(-0.148), max = exp(1.002), meanlog = 0.45, sdlog = 0.19),
        PL          = rlnormTrunc (1, min = exp(0.59), max = exp(1.37), meanlog = 0.98, sdlog = 0.198),
        PL_pup      = rlnormTrunc (1, min = exp(0.53), max = exp(1.30), meanlog = 0.92, sdlog = 0.198),
        PMilkM      = rlnormTrunc (1, min = exp(0.22), max = exp(1.17), meanlog = 0.60, sdlog = 0.19),
        PRest_pup   = rlnormTrunc (1, min = exp(-1.92), max = exp(-1.14), meanlog = -1.53, sdlog = 0.13),
        RAFapi      = rlnormTrunc (1, min = exp(0.80), max = exp(1.955), meanlog = 1.38, sdlog = 0.29),
        VFilC       = rnormTrunc  (1, min = 0.00035, max = 0.00133, mean = 0.00084, sd = 0.00025),
        Vmax_apical_invitro_p = rlnormTrunc (1, min = exp(5.35), max = exp(6.507), meanlog = 5.93, sdlog = 0.29),
        Vmax_apical_invitro = 4141,
        PRest       = 0.03,
        SA_BW       = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_VLC      = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_KMilkC   = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_VL_pup   = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2)
    ))
    
    Gpars_H  <- log(c(
        BW        = rnormTrunc (1, min = 27.89, max = 107.5, mean = 67.7, sd = 20.3),
        SA_Htc    = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_QC     = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_QK_P   = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        Free      = rlnormTrunc (1, min = exp(-6.53), max = exp(-5.38), meanlog = -5.96, sdlog = 0.29),  
        Free_Fet  = rlnormTrunc (1, min = exp(-6.19), max = exp(-5.04), meanlog = -5.61, sdlog = 0.29),
        KeffluxC  = rlnormTrunc (1, min = exp(-4.82), max = exp(-3.67), meanlog = -4.24, sdlog = 0.29), 
        Ktrans1C  = rlnormTrunc (1, min = exp(-0.85), max = exp(0.29), meanlog = -0.28, sdlog = 0.29), 
        Ktrans2C  = rlnormTrunc (1, min = exp(-0.51), max = exp(0.65), meanlog = 0.07, sdlog = 0.29), 
        QCC       = rnormTrunc  (1, min = 6.76, max = 26.04, mean = 16.4, sd = 4.92), 
        Km_apical = 248,
        PPla      = 0.13,
        PL_Fet    = 0.58,
        PRest_Fet = 2.3
    ))    
    
    Lpars_H <- log(c(
        SA_BW     = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        SA_BW_neo = rnormTrunc (1, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
        Free      = rlnormTrunc (1, min = exp(-3.68), max = exp(-2.94), meanlog = -3.32, sdlog = 0.07), 
        Free_neo  = rlnormTrunc (1, min = exp(-4.71), max = exp(-3.87), meanlog = -4.29, sdlog = 0.13),
        KeffluxC  = rlnormTrunc (1, min = exp(-2.52), max = exp(-1.36), meanlog = -1.94, sdlog = 0.19),
        KeffluxC_neo = rlnormTrunc (1, min = exp(-2.52), max = exp(-1.36), meanlog = -1.94, sdlog = 0.19),
        PAMilkC    = rlnormTrunc (1, min = exp(-6.49), max = exp(-5.35), meanlog = -3.6, sdlog = 0.29),
        PRest     = rlnormTrunc (1, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
        PRest_neo = rlnormTrunc (1, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
        QCC       = rnormTrunc (1, min = 6.76, max = 26.04, mean = 16.4, sd = 3.28),
        QKC       = rnormTrunc (1, min = 0.07, max = 0.28, mean = 0.175, sd = 0.035),
        RAFbaso   = rlnormTrunc (1, min = exp(-0.618), max = exp(-0.532), meanlog = -0.043, sdlog = 0.29),
        VKC       = rnormTrunc (1, min = 0.002, max = 0.01, mean = 0.004, sd = 0.001),
        RAFapi   = 0.525
    ))    
    
    

    outdf_A1_R  <- pred.rat(Gpars_R, Lpars_R, 1)[[1]]   %>% filter (Time == 20*24) %>% select(AUC_CPlas)
    outdf_A2_R  <- pred.rat(Gpars_R, Lpars_R, 0.4)[[2]] %>% filter (Time == 20*24) %>% select(AUC_CPlas)
    outdf_A3_R  <- pred.rat(Gpars_R, Lpars_R, 0.1)[[3]] %>% filter (Time == 40*24) %>% select(AUC_CPlas)
    outdf_A4_R  <- pred.rat(Gpars_R, Lpars_R, 0.3)[[4]] %>% filter (Time == 40*24) %>% select(AUC_CPlas)
    
    outdf_A1_H  <- pred.human(Gpars_H, Lpars_H, 1)[[1]]   %>% filter (Time == 40*24*7) %>% select(AUC_CPlas)
    outdf_A2_H  <- pred.human(Gpars_H, Lpars_H, 0.4)[[1]] %>% filter (Time == 40*24*7) %>% select(AUC_CPlas)
    outdf_A3_H  <- pred.human(Gpars_H, Lpars_H, 0.1)[[2]] %>% filter (Time == 70*24*7) %>% select(AUC_CPlas)
    outdf_A4_H  <- pred.human(Gpars_H, Lpars_H, 0.3)[[2]] %>% filter (Time == 70*24*7) %>% select(AUC_CPlas)
    
       ###########
    ## Defined the Monte Carlo simulation function for rat to estiamte the probablistic AUC based on NOAEL
    
    ASC_A1 [[i]]   <- outdf_A1_R/(19*24) # A1 Lau et al. 2003       : (GD2- GD21)    : NOAEL = 1.0 mg/kg-day
    ASC_A2 [[i]]   <- outdf_A2_R/(78*24) # A2 Luebker et al. 2005b  : (PREG - G20)   : NOAEL = 0.4 mg/kg-day
    ASC_A3 [[i]]   <- outdf_A3_R/(98*24) # A3 Luebker et al. 2005a  : (PREG - PND20) : NOAEL = 0.1 mg/kg-day
    ASC_A4 [[i]]   <- outdf_A4_R/(41*24) # A4 Butenhoff et al. 2009 : (GD0 - PND20)  : NOAEL = 0.3 mg/kg-day
    
    ## EPA method
    HED_A1_A [[i]]   <- ((ASC_A1[[i]])*0.000081) 
    HED_A2_A [[i]]   <- ((ASC_A2[[i]])*0.000081) 
    HED_A3_A [[i]]   <- ((ASC_A3[[i]])*0.000081) 
    HED_A4_A [[i]]   <- ((ASC_A4[[i]])*0.000081) 
    
    ## AUC method
    HED_A1_B [[i]]   <- 1*  ((outdf_A1_R) / (outdf_A1_H )) 
    HED_A2_B [[i]]   <- 0.4*((outdf_A2_R) / (outdf_A2_H )) 
    HED_A3_B [[i]]   <- 0.1*((outdf_A3_R) / (outdf_A3_H )) 
    HED_A4_B [[i]]   <- 0.3*((outdf_A4_R) / (outdf_A4_H )) 
    
    
    cat(paste ("iteration = ", i, "\n"))
}


## Estimate the median (95% CI) of ASC (mg/L)
# A1 Lau et al. 2003       : (GD2- GD21)    : NOAEL = 1.0 mg/kg-day
# A2 Luebker et al. 2005b  : (PREG - G20)   : NOAEL = 0.4 mg/kg-day
# A3 Luebker et al. 2005a  : (PREG - PND20) : NOAEL = 0.1 mg/kg-day
# A4 Butenhoff et al. 2009 : (GD0 - PND20)  : NOAEL = 0.3 mg/kg-day
Range_ASC_A1 <- do.call (rbind, ASC_A1)%>%
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95  = quantile (AUC_CPlas , probs = 0.975),
        ci_05  = quantile (AUC_CPlas , probs = 0.025))

Range_ASC_A2 <- do.call (rbind, ASC_A2)%>%  
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95  = quantile (AUC_CPlas , probs = 0.975),
        ci_05  = quantile (AUC_CPlas , probs = 0.025))

Range_ASC_A3 <- do.call (rbind, ASC_A3)%>% 
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95  = quantile (AUC_CPlas , probs = 0.975),
        ci_05  = quantile (AUC_CPlas , probs = 0.025))

Range_ASC_A4 <- do.call (rbind, ASC_A4)%>%  
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95  = quantile (AUC_CPlas , probs = 0.975),
        ci_05  = quantile (AUC_CPlas , probs = 0.025))



## Method 1:ASC 
Range_HED_A1_A <- do.call (rbind, HED_A1_A)%>%  
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95 = quantile (AUC_CPlas , probs = 0.975),
        ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A2_A <- do.call (rbind, HED_A2_A)%>%  
        summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95 = quantile (AUC_CPlas , probs = 0.975),
        ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A3_A <- do.call (rbind, HED_A3_A)%>%   
    summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95 = quantile (AUC_CPlas , probs = 0.975),
        ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A4_A <- do.call (rbind, HED_A4_A)%>%  
    summarize (
        median = quantile (AUC_CPlas , probs = 0.5), 
        ci_95 = quantile (AUC_CPlas , probs = 0.975),
        ci_05 = quantile (AUC_CPlas , probs = 0.025))


### Method 2:AUC ratios

Range_HED_A1_B <- do.call (rbind, HED_A1_B)%>%  
                summarize (
                median = quantile (AUC_CPlas , probs = 0.5), 
                ci_95 = quantile (AUC_CPlas , probs = 0.975),
                ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A2_B <- do.call (rbind, HED_A2_B)%>%  
                summarize (
                median = quantile (AUC_CPlas , probs = 0.5), 
                ci_95 = quantile (AUC_CPlas , probs = 0.975),
                ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A3_B <- do.call (rbind, HED_A3_B)%>%  
                summarize (
                median = quantile (AUC_CPlas , probs = 0.5), 
                ci_95 = quantile (AUC_CPlas , probs = 0.975),
                ci_05 = quantile (AUC_CPlas , probs = 0.025))

Range_HED_A4_B <- do.call (rbind, HED_A4_B)%>%  
                summarize (
                median = quantile (AUC_CPlas , probs = 0.5), 
                ci_95 = quantile (AUC_CPlas , probs = 0.975),
                ci_05 = quantile (AUC_CPlas , probs = 0.025))





