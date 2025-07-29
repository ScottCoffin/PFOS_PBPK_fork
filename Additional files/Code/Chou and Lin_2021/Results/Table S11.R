
##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for ggplot
library(FME)         ## R-package for MCMC simulation and model fitting
library(minpack.lm)  ## R-package for model fitting

## Input mrgsolve-based PBPK Model
source (file = "HMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code)
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Prediction function for gestaiton 
pred.G <- function(DOSE, Route) {
    
    route = Route
    # ## Get out of log domain
    # Gpars <- lapply(Gpars, exp)           ## Return a list of exp (parametrs for gestational model) from log scale
     
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
                    addl = 1 - 1,         ## Addtional doseing 
                    cmt  = route,         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, 24*7*30, 24) 
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        update(atol = 1E-15,  rtol= 1e-8, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    Goutdf = cbind.data.frame (Time    = Gout$time,
                               A_Plas  = (Gout$APlas_free/GDOSEoral)*100,
                               A_ST    = (Gout$AST/GDOSEoral)*100,
                               A_SI    = (Gout$ASI/GDOSEoral)*100,
                               A_Feces = (Gout$Afeces/GDOSEoral)*100,
                               A_Urine = (Gout$Aurine/GDOSEoral)*100,
                               CPlas   = Gout$Plasma*1000)
    
    return (Goutdf) # Return Goutdf
}

## 
Oral_G <- pred.G (1e-6, Route = "AST")
IV_G   <- pred.G (1e-6, Route = "APlas_free")


Lose_CPlas_G <- (1-Oral_G$CPlas[length(Oral_G$CPlas)]/IV_G$CPlas[length(IV_G$CPlas)])*100
Lose_CFPlas_G<-(Oral_G$CFPlas[length(Oral_G$CFPlas)]/IV_G$CFPlas[length(IV_G$CFPlas)]-1)*100

Lose_CPlas_G
Lose_CFPlas_G

## Prediction function for lactational
pred.L <- function(DOSE, Route) {
    
    route = Route
    
    ## Exposure scenario for Lactational exposure; 
    LBW          = 67                     ## Body weight  
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 180                    ## Total dosing/Dose times; Repeat oral dose from PND0 - PND41 (weeks)
    LDOSE        = DOSE                   ## Repeat oral dose  (mg/kg/day);  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose (mg)
    
    
    # To create a data set of 1 subject receiving DOSE every 24 hours
    Lex.oral <- ev (ID   = 1,             ## One individual
                    amt  = LDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = LTDOSE - 1,    ## Addtional doseing 
                    cmt  = route,         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 24) ## Simulation time from PND0 to PND74 with dosing of LTDOSE - 1+ 2 days (No dosing)
    
    
    Lout <- Lmod_H %>% # Lactational model
        update(atol = 1E-15,  rtol= 1e-8, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    Loutdf = cbind.data.frame (Time    = Lout$time,
                               A_Plas  = (Lout$APlas_free/LDOSEoral)*100,
                               A_ST    = (Lout$AST/LDOSEoral)*100,
                               A_SI    = (Lout$ASI/LDOSEoral)*100,
                               A_Feces = (Lout$Afeces/LDOSEoral)*100,
                               A_Urine = (Lout$Aurine/LDOSEoral)*100,
                               CPlas   = Lout$Plasma*1000)
    
    return (Loutdf)
    
}

## 
Oral_L <- pred.L (1e-6, Route = "AST")
IV_L   <- pred.L (1e-6, Route = "APlas_free")


Lose_CPlas_L <- (Oral_L$CPlas[length(Oral_L$CPlas)]/IV_L$CPlas[length(IV_L$CPlas)] - 1)*100
Lose_CPlas_pup_L <-(1 - Oral_L$CPlas_pup[length(Oral_L$CPlas_pup)]/IV_L$CPlas_pup[length(IV_L$CPlas_pup)])*100

Lose_CPlas_L
Lose_CPlas_pup_L



