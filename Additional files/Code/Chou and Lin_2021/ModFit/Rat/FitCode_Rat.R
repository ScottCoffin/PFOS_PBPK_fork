##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for MCMC simulation and model fitting
library(minpack.lm)  ## R-package for model fitting

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_R <- mcode ("PreG_RatPBPK.code", PreG_RatPBPK.code) # bulid the pre-pregnant model
Gmod_R    <- mcode ("GRatPBPK.code", GRatPBPK.code) # bulid the gestational model
Lmod_R    <- mcode ("LRatPBPK.code", LRatPBPK.code) # bulid the lactational model

## Loading datasets for model calibration
Data_0    <- read.csv(file = "Data_R.csv")

#==================================================================================================
# Model calibration for gestational PBPK model based on the data of Thibodeaux et al., 2003                                                                                        #
# Dose regimen: 1, 2, 3, 5, 10 mg/kg/day, exposure from GD2 - GD20
# Abreviation: maternal plasma (MP), maternal liver (ML), fetal liver (FL)
# A1. : Pregnant SD rat oral dose to 1 mg/kg,  matrix: MP, Sampling time: GD7, 14 and 21                                                     
# A2. : Pregnant SD rat oral dose to 2 mg/kg,  matrix: MP, Sampling time: GD7, 14 and 21                                   
# A3. : Pregnant SD rat oral dose to 3 mg/kg,  matrix: MP, Sampling time: GD7, 14 and 21                                   
# A4. : Pregnant SD rat oral dose to 5 mg/kg,  matrix: MP, Sampling time: GD7, 14 and 21                                   
# A5. : Pregnant SD rat oral dose to 10 mg/kg, matrix: MP, Sampling time: GD7, 14 and 21                                   
# A6. : Pregnant SD rat oral dose to 1 mg/kg,  matrix: ML, Sampling time: GD21                                                      
# A7. : Pregnant SD rat oral dose to 2 mg/kg,  matrix: ML, Sampling time: GD21                                   
# A8. : Pregnant SD rat oral dose to 3 mg/kg,  matrix: ML, Sampling time: GD21                                    
# A9. : Pregnant SD rat oral dose to 5 mg/kg,  matrix: ML, Sampling time: GD21                                    
# A10 : Pregnant SD rat oral dose to 10 mg/kg, matrix: ML, Sampling time: GD21                                   
# A11.: Pregnant SD rat oral dose to 1 mg/kg,  matrix: FL, Sampling time: GD21                                                      
# A12.: Pregnant SD rat oral dose to 2 mg/kg,  matrix: FL, Sampling time: GD21                                   
# A13.: Pregnant SD rat oral dose to 3 mg/kg,  matrix: FL, Sampling time: GD21                                    
# A14.: Pregnant SD rat oral dose to 5 mg/kg,  matrix: FL, Sampling time: GD21                                    
# A15 : Pregnant SD rat oral dose to 10 mg/kg, matrix: FL, Sampling time: GD21
#===================================================================================================

## Read these datasets and later used in model calibration
OBS.A1  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 1 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A2  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 2 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A3  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 3 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A4  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 5 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A5  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 10& Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A6  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 1 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A7  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 2 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A8  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 3 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A9  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 5 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A10 <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 10 & Time!= 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A11 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 1 & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS.A12 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 2 & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS.A13 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 3 & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS.A14 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 5 & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS.A15 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 10 & Time != 0)%>% select(Time = "Time", CFL   = "Conc")

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
    
    Goutdf <- Goutdf %>% filter (Time > 0) # filter the value at time = 0
    return (Goutdf) # Return Goutdf
}

    
## Create a cost fuction and later used in model calibration 
## Estimate the model residual with experimental data by modCost function (from FME package)
Gcost<-function (pars){
    outdf.A1 <- pred.A (Gpars = pars,  DOSE = 1)
    outdf.A2 <- pred.A (Gpars = pars,  DOSE = 2)
    outdf.A3 <- pred.A (Gpars = pars,  DOSE = 3)
    outdf.A4 <- pred.A (Gpars = pars,  DOSE = 5)
    outdf.A5 <- pred.A (Gpars = pars,  DOSE = 10)
    
    cost<- modCost  (model = outdf.A1, obs = OBS.A1, x ="Time")
    cost<- modCost  (model = outdf.A2, obs = OBS.A2, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS.A3, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A4, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS.A5, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A1, obs = OBS.A6, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS.A7, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS.A8, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A9, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS.A10, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A1, obs = OBS.A11, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS.A12, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS.A13, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A14, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS.A15, x ="Time", cost = cost)
    
    return(cost)
}


## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
Gtheta.int <- log(c(
    Vmax_baso_invitro              = 393.45,                      
    Km_baso                        = 27.2,                        
    Vmax_apical_invitro            = 1808,                        
    Km_apical                      = 278,                        
    RAFbaso                        = 1.90,                        
    RAFapi                         = 4.15,                       
    Kdif                           = 5.1e-4,                       
    KeffluxC                       = 2.09,                        
    KbileC                         = 0.0026,                       
    KurineC                        = 1.60,                        
    Free                           = 0.022,                        
    PL                             = 3.66,                        
    PK                             = 0.8,                         
    PM                             = 0.16,
    PF                             = 0.13,
    PRest                          = 0.22,                         
    PPla                           = 0.41,
    K0C                            = 1,                           
    Kabsc                          = 2.12,                        
    KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 0.46,
    Ktrans2C                       = 1,
    Ktrans3C                       = 0.008,
    Ktrans4C                       = 0.001,
    Free_Fet                       = 0.022,
    PL_Fet                         = 3.66,
    PRest_Fet                      = 0.22
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = Gcost, parms = Gtheta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(summary(SnsPlasma))

## Selected senstive parameters
Gtheta <- Gtheta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
Gtheta 

## Selected parameters
Gtheta <- log(c(
    Vmax_baso_invitro              = 393.45,                      
    #Km_baso                        = 27.2,                        
    #Vmax_apical_invitro            = 1808,                        
    Km_apical                      = 278,                        
    #RAFbaso                        = 1.90,                        
    #RAFapi                         = 4.15,                       
    #Kdif                           = 5.1e-4,                       
    #KeffluxC                       = 2.09,                        
    KbileC                         = 0.0026,                       
    #KurineC                        = 1.60,                        
    Free                           = 0.022,                        
    PL                             = 3.66,                        
    #PK                             = 0.8,                         
    #PM                             = 0.16,
    #PF                             = 0.13,
    #PRest                          = 0.22,                         
    #PPla                           = 0.41,
    #K0C                            = 1,                           
    #Kabsc                          = 2.12,                        
    #KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 0.46,
    #Ktrans2C                       = 1,
    Ktrans3C                       = 0.008,
    #Ktrans4C                       = 0.001,
    #Free_Fet                       = 0.022,
    PL_Fet                         = 3.66,
    PRest_Fet                      = 0.22
))

## PBPK model calibration 
## Least squres fit using levenberg-marquart (method "Marq") algorithm
GFit <- modFit(f = Gcost, p = Gtheta, method ="Marq",
               control = nls.lm.control(nprint = 1))

## Summary of fitting results
summary(GFit) 
exp(GFit$par)   ## Get the arithmetic value out of the log domain
Gcost(GFit$par) ## Check the results 

#==================================================================================================================== 
# Model calibration for lactational PBPK model based on the data of Chang et al., 2009          
# Exposure scenario: Dose regimen: 0.1, 0.3, 1 mg/kg/day; exposure from GD0 (day positive for mating) through PND20.                                                                                 #
# Abreviation: D: Dam, P: Pup, MP: maternal plasma, ML: maternal liver, NP: neonatal plasma, NL: neonatal liver
# Note: Oringinal data included the samples at GD20, PND4 and PND72. This study used the samples at GD20 and PND4 for model calibration and not include the sample at PND72 due to the model limitation
# B1. : Pregnant SD rat oral dose to 1 mg/kg,  matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                                            
# B2. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                           
# B3. : Pregnant SD rat oral dose to 0.1 mg/kg,    matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                            
# B4. : Pregnant SD rat oral dose to 1 mg/kg,  matrix: NP, Sampling time: GD20, PND4, PND21 and PND72                                             
# B5. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: NP, Sampling time: GD20, PND4, PND21 and PND72                           
# B6. : Pregnant SD rat oral dose to 0.1 mg/kg,    matrix: NP, Sampling time: GD20, PND4, PND21 and PND72  
# B7. : Pregnant SD rat oral dose to 1 mg/kg,  matrix: NL, Sampling time: GD20, PND4, PND21 and PND72                                             
# B8. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: NL, Sampling time: GD20, PND4, PND21 and PND72                           
# B9. : Pregnant SD rat oral dose to 0.1 mg/kg,    matrix: NL, Sampling time: GD20, PND4, PND21 and PND72  
#======================================================================================================================                                                   

## Read the data and later used in model calibration
## Filter the dataset and get the samples during lactation (Time > 480)
OBS.B1  <- Data_0 %>% filter(Study == 2 & Sample == "MP" & Dose == 1 )  %>% select(Time = "Time", CPlas = "Conc")%>%filter(Time > 480) 
OBS.B2  <- Data_0 %>% filter(Study == 2 & Sample == "MP" & Dose == 0.3) %>% select(Time = "Time", CPlas = "Conc")%>%filter(Time > 480)
OBS.B3  <- Data_0 %>% filter(Study == 2 & Sample == "MP" & Dose == 0.1) %>% select(Time = "Time", CPlas = "Conc")%>%filter(Time > 480)
OBS.B4  <- Data_0 %>% filter(Study == 2 & Sample == "NP" & Dose == 1)   %>% select(Time = "Time", CPlas_pup = "Conc")%>%filter(Time > 480)
OBS.B5  <- Data_0 %>% filter(Study == 2 & Sample == "NP" & Dose == 0.3) %>% select(Time = "Time", CPlas_pup = "Conc")%>%filter(Time > 480)
OBS.B6  <- Data_0 %>% filter(Study == 2 & Sample == "NP" & Dose == 0.1) %>% select(Time = "Time", CPlas_pup = "Conc")%>%filter(Time > 480)
OBS.B7  <- Data_0 %>% filter(Study == 2 & Sample == "NL" & Dose == 1)   %>% select(Time = "Time", CL_pup = "Conc")%>%filter(Time > 480)
OBS.B8  <- Data_0 %>% filter(Study == 2 & Sample == "NL" & Dose == 0.3) %>% select(Time = "Time", CL_pup = "Conc")%>%filter(Time > 480)
OBS.B9  <- Data_0 %>% filter(Study == 2 & Sample == "NL" & Dose == 0.1) %>% select(Time = "Time", CL_pup = "Conc")%>%filter(Time > 480)


## Exposure scenario. Dosing from GD0 - PND20; the dams were sacrificed on GD 20 and PND20 based on study design of Chang et al. (2009) 
pred.B <- function(Lpars, DOSE) {
    
    ## Get out of log domain
    Gpars <- lapply(GFit$par, exp)        ## Return a list of exp (parametrs for gestational model) from log scale
    Lpars <- lapply(Lpars, exp)           ## Return a list of exp (parametrs for lactational model) from log scale
    
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
    
    Goutdf = cbind.data.frame(Time      = Gout$time,       # output time (hours)
                              CPlas     = Gout$Plasma,     # Maternal plasma
                              CL        = Gout$Liver,      # Maternal liver
                              CPlas_pup = Gout$Plasma_Fet, # Fetal plasma
                              CL_pup    = Gout$Liver_Fet)  # Fetal liver
    
    # Extract the predicted value at GD21 and then used as initial values of lactational model
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

## Cost fuction (FME) 
## Estimate the model residual by modCost function
Lcost<-function (pars){
    outdf.A <- pred.B (Lpars = pars,  DOSE = 1)
    outdf.B <- pred.B (Lpars = pars,  DOSE = 0.3)
    outdf.C <- pred.B (Lpars = pars,  DOSE = 0.1)
    
    cost<- modCost  (model = outdf.A, obs = OBS.B1, x ="Time")
    cost<- modCost  (model = outdf.B, obs = OBS.B2, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C, obs = OBS.B3, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A, obs = OBS.B4, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B, obs = OBS.B5, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C, obs = OBS.B6, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A, obs = OBS.B7, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B, obs = OBS.B8, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C, obs = OBS.B9, x ="Time", cost = cost)
    
    return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
Ltheta.int <- log(c(
    Vmax_baso_invitro              = 393.45,
    Km_baso                        = 27.2,
    Vmax_apical_invitro            = 1808,
    Km_apical                      = 278,
    RAFbaso                        = 1.90,
    RAFapi                         = 4.15,
    KeffluxC                       = 2.09,
    KbileC                         = 0.0026,
    KurineC                        = 1.6,
    Free                           = 0.022,
    PL                             = 3.66,
    PK                             = 0.8,
    PM                             = 0.16,
    PF                             = 0.13,
    PRest                          = 0.26,
    PMilkM                         = 1.9,
    PMilkP                         = 0.11,
    PAMilkC                        = 0.5,
    K0C                            = 1,
    KabsC                          = 2.12,
    KunabsC                        = 7.05e-5,
    Kdif                           = 5.1e-4,
    KMilk0                         = 0.0268,
    Free_pup                       = 0.022,
    PL_pup                         = 3.66,
    PK_pup                         = 0.8,
    PRest_pup                      = 0.22,
    KabsC_pup                      = 2.12,
    KbileC_pup                     = 0.35,
    KurineC_pup                    = 1.6,
    KeffluxC_p                     = 2.09,
    Kdif_pup_0                     = 0.001,
    Vmax_baso_invitro_p            = 393.45,
    Km_baso_p                      = 27.2,
    Vmax_apical_invitro_p          = 1808,
    Km_apical_p                    = 278
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = Lcost, parms = Ltheta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(summary(SnsPlasma))


## Selected senstive parameters; 
Ltheta <- Ltheta.int[abs(Sen$Min) > 1.2*mean(abs(Sen$Min))]
Ltheta 

## Selected senstive parmaeters and relavant parameters
Ltheta <- log(c(
    #Vmax_baso_invitro              = 393.45,
    #Km_baso                        = 27.2,
    Vmax_apical_invitro            = 1808,
    #Km_apical                      = 278,
    #RAFbaso                        = 1.90,
    #RAFapi                         = 4.15,
    #KeffluxC                       = 2.09,
    #KbileC                         = 0.0026,
    #KurineC                        = 1.6,
    Free                           = 0.022,
    PL                             = 3.66,
    #PK                             = 0.8,
    #PM                             = 0.16,
    #PF                             = 0.13,
    PRest                          = 0.26,
    #PMilkM                         = 1.9,
    #PMilkP                         = 0.11,
    #PAMilkC                        = 0.5,
    #K0C                            = 1,
    #KabsC                          = 2.12,
    #KunabsC                        = 7.05e-5,
    #Kdif                           = 5.1e-4,
    KMilk0                         = 0.0268,
    #Free_pup                       = 0.022,
    PL_pup                         = 3.66
    #PK_pup                         = 0.8,
    #PRest_pup                      = 0.22
    #KabsC_pup                      = 2.12,
    #KbileC_pup                     = 0.0026,
    #KurineC_pup                    = 1.6
    #KeffluxC_p                     = 2.09,
    #Kdif_pup_0                     = 0.001,
    #Vmax_baso_invitro_p            = 393.45,
    #Km_baso_p                      = 27.2,
    #Vmax_apical_invitro_p          = 1808
    #Km_apical_p                    = 278
))

## PBPK model fitting 
## Least squres fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters

LFit <- modFit(f = Lcost, p = Ltheta, method ="Marq",
               control = nls.lm.control(nprint = 1))


summary(LFit)    ## Summary of fit 
exp(LFit$par)    ## Get the arithmetic value out of the log domain
                                       
## Check the model residuces
Lcost(LFit$par)

## Save the fitting results to RDS files
saveRDS(GFit, file = "GFit_R.rds")
saveRDS(LFit, file = "LFit_R.rds")








