##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for ggplot
library(FME)         ## R-package for MCMC simulation and model fitting
library(minpack.lm)  ## R-package for model fitting

## Input mrgsolve-based PBPK Model
source (file = "HMod.R")  ## PBPK model for human


## Build mrgsolve-based PBPK Model
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code)
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Loading datasets
Data_0    <- read.csv(file = "Data_H.csv")
Data_Cal  <- Data_0 %>% filter(Cal == 1) 

#=============================================================================================
# Model calibration for gestational PBPK model
# Abbrevation: MP: Maternal plasma, CB: cord blood, Pla: placenta; Fli: Fetal liver
# A1. (Inoue et al., 2004) : Japan population;   Matrix: MP, CB        
# A2. (Fei et al., 2007)   : Danish population;  Matrix: MP, CB           
# A3. (Kato et al., 2014)  : U.S. population;    Matrix: MP, CB          
# A4. (Pan et al., 2017)   : Chinese population; Matrix: MP, CB         
# A5. (Mamsen et al., 2019): Sweden population;  Matrix: MP, Pla
# A6. (Midasch et al., 2007) : German population; Matrix: MP, CB
#=============================================================================================

## Read the data and later used in model calibration 
OBS.A1_MP     <- Data_Cal %>% filter(Study == 1 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A1_CB     <- Data_Cal %>% filter(Study == 1 & Matrix == "CB") %>% select(Time = "Time",  CFPlas  = "Conc")
OBS.A2_MP     <- Data_Cal %>% filter(Study == 2 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A2_CB     <- Data_Cal %>% filter(Study == 2 & Matrix == "CB") %>% select(Time = "Time",  CFPlas  = "Conc")
OBS.A3_MP     <- Data_Cal %>% filter(Study == 3 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A3_CB     <- Data_Cal %>% filter(Study == 3 & Matrix == "CB") %>% select(Time = "Time",  CFPlas  = "Conc")
OBS.A4_MP     <- Data_Cal %>% filter(Study == 4 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A4_CB     <- Data_Cal %>% filter(Study == 4 & Matrix == "CB") %>% select(Time = "Time",  CFPlas  = "Conc")
OBS.A5_MP     <- Data_Cal %>% filter(Study == 5 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A5_Pla    <- Data_Cal %>% filter(Study == 5 & Matrix == "Pla") %>% select(Time = "Time", CPla    = "Conc")
OBS.A5_Fli    <- Data_Cal %>% filter(Study == 5 & Matrix == "Fli") %>% select(Time = "Time", CFL     = "Conc")
OBS.A6_MP     <- Data_Cal %>% filter(Study == 6 & Matrix == "MP") %>% select(Time = "Time",  CPlas   = "Conc")
OBS.A6_CB     <- Data_Cal %>% filter(Study == 6 & Matrix == "CB") %>% select(Time = "Time",  CFPlas  = "Conc")

## Define the prediction function 
## Exposure scenario:  Dosing druing pregnnacy from pre-pregnant to gestational stage
pred.preG <- function (DOSE) {
    
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
    ## a constant exposure (per kg of bodyweight) from birth until the year 2000 (Assumed 20 years old)
    PDOSEoral_1 = DOSE ## assumed the daily oral dose before 2000
    ex_1 <- mutate (ex, 
                    amt = case_when (YEAR < 20 ~ PDOSEoral_1*BW,
                                         YEAR >= 20 ~ PDOSEoral_1*BW), 
                    cmt  = "AST", 
                    ii = 24, 
                    evid = 1, 
                    time = DAY)
    
    
    out <- PreGmod_H %>%
        update(atol = 1E-3, maxsteps = 500000) %>%          
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

Init_A1 <- pred.preG (DOSE = DOSE_A1) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A2 <- pred.preG (DOSE = DOSE_A2) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A3 <- pred.preG (DOSE = DOSE_A3) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A4 <- pred.preG (DOSE = DOSE_A4) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A5 <- pred.preG (DOSE = DOSE_A5) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
Init_A6 <- pred.preG (DOSE = DOSE_A6) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))

## Exposure sceinario 
pred.G <- function(Gpars, DOSE, Init) {
    
    ## Get out of log domain
    Gpars <- lapply(Gpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 60                     ## Body weight 
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; 
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time 
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>%  ## Input the intial concentration
        param (Gpars) %>%                 ## Update the parameter list with Gpars
        update(atol = 1E-3,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
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

    
## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)
Gcost<-function (pars){
    outdf.A1 <- pred.G (Gpars = pars,  DOSE = DOSE_A1, Init = Init_A1)
    outdf.A2 <- pred.G (Gpars = pars,  DOSE = DOSE_A2, Init = Init_A2)
    outdf.A3 <- pred.G (Gpars = pars,  DOSE = DOSE_A3, Init = Init_A3)
    outdf.A4 <- pred.G (Gpars = pars,  DOSE = DOSE_A4, Init = Init_A4)
    outdf.A5 <- pred.G (Gpars = pars,  DOSE = DOSE_A5, Init = Init_A5)
    outdf.A6 <- pred.G (Gpars = pars,  DOSE = DOSE_A6, Init = Init_A6)
    
    cost<- modCost  (model = outdf.A1, obs = OBS.A1_MP, x ="Time")
    cost<- modCost  (model = outdf.A1, obs = OBS.A1_CB, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS.A2_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS.A2_CB, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS.A3_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS.A3_CB, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A4_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A4_CB, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A5_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A5_Pla, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS.A5_Fli, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A6, obs = OBS.A6_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A6, obs = OBS.A6_CB, x ="Time", cost = cost)
    
    return(cost)
}

## Local sensitivity analysis
## initial parmaeters
Gtheta.int <- log(c(
    Vmax_baso_invitro              = 479,                      
    Km_baso                        = 20.1,                        
    Vmax_apical_invitro            = 51803,                        
    Km_apical                      = 64.4,                        
    RAFapi                         = 0.001,                       
    RAFbaso                        = 1,                        
    KeffluxC                       = 0.15,                        
    KbileC                         = 1.3e-4,                       
    KurineC                        = 0.096,                        
    Free                           = 0.014,                        
    PL                             = 2.03,                        
    PK                             = 1.26,                         
    PM                             = 0.13,
    PF                             = 0.16,
    PRest                          = 0.20,                         
    PPla                           = 0.41,
    K0C                            = 1,                           
    Kabsc                          = 2.12,                        
    Kdif                           = 0.001,                       
    KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 0.46,
    Ktrans2C                       = 1.01,
    Ktrans3C                       = 0.006,
    Ktrans4C                       = 0.001,
    Free_Fet                       = 0.014,
    PL_Fet                         = 2.03,
    PRest_Fet                      = 0.20
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = Gcost, parms = Gtheta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(summary(SnsPlasma))

## Selected senstive parameters;
## Selecte the parameters that senstive scores higher than 1.2x average socres of all parameters
Gtheta <- Gtheta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))] 
Gtheta

## 
Gtheta <- log(c(
    #Vmax_baso_invitro              = 479,                      
    #Km_baso                        = 20.1,                        
    #Vmax_apical_invitro            = 51803,                        
    Km_apical                      = 64.4,                        
    #RAFapi                         = 0.001,                       
    #RAFbaso                        = 1,                        
    KeffluxC                       = 0.15,                        
    #KbileC                         = 1.3e-4,                       
    #KurineC                        = 0.096,                        
    Free                           = 0.014,                        
    #PL                             = 2.03,                        
    #PK                             = 1.26,                         
    #PM                             = 0.13,
    #PF                             = 0.16,
    #PRest                          = 0.20,                         
    PPla                           = 0.41,
    #K0C                            = 1,                           
    #Kabsc                          = 2.12,                        
    #Kdif                           = 0.001,                       
    #KunabsC                        = 7.05e-5,                     
    Ktrans1C                       = 0.46,
    Ktrans2C                       = 1.01,
    #Ktrans3C                       = 0.006,
    #Ktrans4C                       = 0.001,
    Free_Fet                       = 0.014,
    PL_Fet                         = 2.03,
    PRest_Fet                      = 0.20
))


## PBPK model fitting 
## Least squres fit using levenberg-marquart (method "Marq") algorithm
GFit <- modFit(f = Gcost, p = Gtheta, method ="Marq",
               control = nls.lm.control(nprint = 1))

summary(GFit)                                  ## Summary of fit 
exp(GFit$par)                                  ## Get the arithmetic value out of the log domain
Gcost(GFit$par)  # Check


#==================================================================================================
# Model calibration for lactational PBPK model based on human biomonitoring data
# Abrreviation: MP: Maternal plasma, NP: Neonatal plasma
# B1. (Karrman et al. (2007)       : Sweden population;  Matrix: MP, Milk                                 
# B2. (Von Ehrenstein et al., 2009): U.S. population;    Matrix: MP                       
# B3. (Fromme et al., 2010)        : Gemerny population; Matrix: MP, NP     
# B4. (Lee et al., 2018)           : Sourth Korean;      Matrix: Milk

## Read the data and later used in model calibration and evluation
OBS.B1_MP   <- Data_Cal %>% filter(Study == 16 & Matrix == "MP")   %>% select(Time = "Time", CPlas = "Conc")
OBS.B1_Milk <- Data_Cal %>% filter(Study == 16 & Matrix == "Milk") %>% select(Time = "Time", CMilk = "Conc")
OBS.B2_MP   <- Data_Cal %>% filter(Study == 17 & Matrix == "MP")   %>% select(Time = "Time", CPlas = "Conc")
OBS.B3_MP   <- Data_Cal %>% filter(Study == 18 & Matrix == "MP")   %>% select(Time = "Time", CPlas = "Conc")
OBS.B3_NP   <- Data_Cal %>% filter(Study == 18 & Matrix == "NP")   %>% select(Time = "Time", CPlas_pup = "Conc")
OBS.B4_Milk <- Data_Cal %>% filter(Study == 19 & Matrix == "Milk") %>% select(Time = "Time", CMilk = "Conc")

## Prediction function from pre-pregnant to gestational exposure
pred.G2 <- function(DOSE) {
    
    Init <- pred.preG (DOSE) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ## Get out of log domain
    Gpars <- lapply(GFit$par, exp)        ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 67                     ## Body weight before gestation 
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; 
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSE every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 0,             ## Dossed strat time 
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time 
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>% ## Input the intial concentrations 
        param (Gpars) %>% ## Update the parameter list with Gpars
        update(atol = 1E-3,  maxsteps = 50000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
    
    Init_G <- Gout %>% filter (time == 40*24*7) ## Extract the concentration at GA40
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
    
    ## Exposure scenario for Lactational exposure; 
    LBW          = 67                     ## Rat body weight 
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 41*7                   ## Total dosing/Dose times; 
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
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 24) ## Simulation time 
    
    
    Lout <- Lmod_H %>% # Lactational model
        init(APlas_free =  Init_G$APlas_free, APTC =  Init_G$APTC, AFil =  Init_G$AFil, AKb =  Init_G$AKb, ARest =  Init_G$ARest,
             AL =  Init_G$AL, AM =  Init_G$AM, AF =  Init_G$AF, A_baso =  Init_G$A_baso, A_apical = Init_G$A_apical, Adif =  Init_G$Adif,
             Aefflux =  Init_G$Aefflux, APlas_free_neo =  Init_G$APlas_Fet_free, 
             AL_neo =  Init_G$AL_Fet, ARest_neo =  Init_G$ARest_Fet) %>% # Input the intial concentration
        param (Lpars) %>% # Update the parameter list with pars
        update(atol = 1E-3, maxsteps = 500000) %>% # Atol: absolute tolerance parameter; hmax: The maximum step size; maxstep: maximum number of steps the solver will take when advancing from one time to the next         
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    ## Extract the concentration 
    Loutdf = cbind.data.frame(Time   = Lout$time + 40*7*24, 
                              CPlas  = Lout$Plasma*1000, 
                              CPlas_pup = Lout$CPneo*1000,
                              CMilk = Lout$Milk*1000)
    
    outdf  <- rbind.data.frame (Loutdf)
    return (Loutdf)
    
}

## Cost fuction (from FME package) 
## Estimate the model residual by modCost function
Lcost<-function (pars){
    outdf.B1 <- pred.L (Lpars = pars,  DOSE = DOSE_B1, Init_G1)
    outdf.B2 <- pred.L (Lpars = pars,  DOSE = DOSE_B2, Init_G2)
    outdf.B3 <- pred.L (Lpars = pars,  DOSE = DOSE_B3, Init_G3)
    outdf.B4 <- pred.L (Lpars = pars,  DOSE = DOSE_B4, Init_G4)
    
    cost<- modCost  (model = outdf.B1, obs = OBS.B1_MP, x ="Time")
    cost<- modCost  (model = outdf.B1, obs = OBS.B1_Milk, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B2, obs = OBS.B2_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B3, obs = OBS.B3_MP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B3, obs = OBS.B3_NP, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B4, obs = OBS.B4_Milk, x ="Time", cost = cost)
    
    return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
Ltheta.int <- log(c(
    Vmax_baso_invitro              = 479,
    Km_baso                        = 20.1,
    Vmax_apical_invitro            = 51803,
    Km_apical                      = 64.4,
    RAFbaso                        = 1,
    RAFapi                         = 0.001,
    KeffluxC                       = 0.150,
    KbileC                         = 1.3e-4,
    KurineC                        = 0.096,
    Free                           = 0.014,
    PL                             = 2.03,
    PK                             = 1.26,
    PM                             = 0.16,
    PF                             = 0.13,
    PRest                          = 0.20,
    PMilkM                         = 1.9,
    PMilkP                         = 0.01,
    PAMilkC                        = 0.5,
    K0C                            = 1,
    KabsC                          = 2.12,
    KunabsC                        = 7.05e-5,
    Kdif                           = 0.001,
    Free_neo                       = 0.014,
    PL_neo                         = 2.03,
    PK_neo                         = 1.26,
    PRest_neo                      = 0.20,
    Vmax_baso_invitro_neo          = 479,
    Km_baso_neo                    = 20.1,
    Vmax_apical_invitro_neo        = 51803,
    Km_apical_neo                  = 64.4,
    KeffluxC_neo                   = 0.150,
    Kdif_neo                       = 0.001,
    KabsC_neo                      = 31.3,
    KbileC_neo                     = 1.3e-4,
    KurineC_neo                    = 0.001
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = Lcost, parms = Ltheta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(summary(SnsPlasma))

## Selected senstive parameters
## Selecte the parameters that senstive scores higher than 2x average socres of all parameters
Ltheta <- Ltheta.int[abs(Sen$Min) > 2*mean(abs(Sen$Min))]
Ltheta 

Ltheta <- log(c(
    #Vmax_baso_invitro              = 479,
    #Km_baso                        = 20.1,
    #Vmax_apical_invitro            = 51803,
    #Km_apical                      = 64.4,
    #RAFbaso                        = 1,
    RAFapi                         = 0.001,
    #KeffluxC                       = 0.150,
    #KbileC                         = 1.3e-4,
    #KurineC                        = 0.096,
    Free                           = 0.014,
    #PL                             = 2.03,
    #PK                             = 1.26,
    #PM                             = 0.16,
    #PF                             = 0.13,
    #PRest                          = 0.20,
    #PMilkM                         = 1.9,
    #PMilkP                         = 0.01,
    PAMilkC                        = 0.5
    #K0C                            = 1,
    #KabsC                          = 2.12,
    #KunabsC                        = 7.05e-5,
    #Kdif                           = 0.001,
    #Free_neo                       = 0.014
    #PL_neo                         = 2.03,
    #PK_neo                         = 1.26,
    #PRest_neo                      = 0.20,
    #Vmax_baso_invitro_neo          = 479,
    #Km_baso_neo                    = 20.1,
    #Vmax_apical_invitro_neo        = 51803,
    #Km_apical_neo                  = 64.4,
    #KeffluxC_neo                   = 0.150,
    #Kdif_neo                       = 0.001,
    #KabsC_neo                      = 31.3,
    #KbileC_neo                     = 1.3e-4,
))

## PBPK model fitting 
## Least squres fit using levenberg-marquart (method "Marq") algorithm
## training gestational parameters
LFit <- modFit(f = Lcost, p = Ltheta, method ="Marq",
               control = nls.lm.control(nprint = 1))

summary(LFit)                                  ## Summary of fit 
exp(LFit$par)                                  ## Get the arithmetic value out of the log domain

## Check results
Lcost(LFit$par)

## Save the fitting results to RDS files
saveRDS(GFit, file = "GFit_H.rds")
saveRDS(LFit, file = "LFit_H.rds")



