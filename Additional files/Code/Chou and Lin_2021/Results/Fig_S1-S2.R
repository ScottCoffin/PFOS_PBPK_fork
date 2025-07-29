##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(gridExtra)     ## R-package for function "grid.arrange" 
library(EnvStats)    ## R-package for defined the truncated distribution

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_R <- mcode ("PreG_RatPBPK.code", PreG_RatPBPK.code)
Gmod_R    <- mcode ("GRatPBPK.code", GRatPBPK.code)
Lmod_R    <- mcode ("LRatPBPK.code", LRatPBPK.code)

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

## Loading datasets
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
OBS.A1  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 1 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD= SD)
OBS.A2  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 2 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = SD)
OBS.A3  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 3 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = SD)
OBS.A4  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 5 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = SD)
OBS.A5  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 10& Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = SD)
OBS.A6  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 1 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = SD)
OBS.A7  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 2 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = SD)
OBS.A8  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 3 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = SD)
OBS.A9  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 5 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = SD)
OBS.A10 <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 10 & Time!= 0) %>% select(Time = "Time", CL    = "Conc", SD = SD)
OBS.A11 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 1 & Time != 0) %>% select(Time = "Time", CFL   = "Conc", SD = SD)
OBS.A12 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 2 & Time != 0) %>% select(Time = "Time", CFL   = "Conc", SD = SD)
OBS.A13 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 3 & Time != 0) %>% select(Time = "Time", CFL   = "Conc", SD = SD)
OBS.A14 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 5 & Time != 0) %>% select(Time = "Time", CFL   = "Conc", SD = SD)
OBS.A15 <- Data_0 %>% filter(Study == 1 & Sample == "FL" & Dose == 10 & Time != 0)%>% select(Time = "Time", CFL   = "Conc", SD = SD)

###############################################################
Pop_Rat_G_A <- function (pars, DOSE, N) {
    
    pars <- exp(pars)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.225                  ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 21                     ## Total dosing/Dose times; Repeat oral dose from GD2 - GD21
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    ## Repeat oral dose of GDOSE mg/kg/day;  
    ## Amount of oral dose
    idata <- 
        tibble(ID = 1:N) %>% 
        mutate( BW        = rnormTrunc (N, min = 0.15, max = 0.57, mean = 0.36, sd = 0.108),
                VLC       = rnormTrunc  (N, min = 0.01, max = 0.06, mean = 0.04, sd = 0.01),
                KbileC    = rlnormTrunc (N, min = exp(-5.38), max = exp(-4.63), meanlog = -5.00, sdlog = 0.19),
                PL        = rlnormTrunc (N, min = exp(0.88), max = exp(1.39), meanlog = 1.13, sdlog = 0.13),
                PRest     = rlnormTrunc (N, min = exp(-1.79), max = exp(-1.28), meanlog = -1.53, sdlog = 0.13),
                Ktrans1C  = rlnormTrunc (N, min = exp(-0.18), max = exp(0.58), meanlog = 0.20, sdlog = 0.19),
                Ktrans2C  = rlnormTrunc (N, min = exp(-0.42), max = exp(0.34), meanlog = -0.04, sdlog = 0.19),
                Ktrans3C  = 0.23,
                Free      = 0.019,
                Vmax_baso_invitro = 221,
                Km_baso   = 19.9,
                Free_Fet  = 0.022,
                PL_Fet    = 2.55,
                GDOSEoral = GDOSE*BW
        )   
    
    
    Gex.oral_1 <- ev (ID   = 1:N,             ## One individual
                      time = 24*2,          ## Dossed strat time (GD2)
                      amt  = idata$GDOSEoral,     ## Amount of dose 
                      ii   = tinterval,     ## Time interval
                      addl = GTDOSE - 1,    ## Addtional doseing 
                      cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                      replicate = FALSE)    ## No replicate
    
    Gex.oral_2 <- ev (ID   = 1,             ## One individual
                      time = 24*2,          ## Dossed strat time (GD2)
                      amt  = GDOSEoral,     ## Amount of dose 
                      ii   = tinterval,     ## Time interval
                      addl = GTDOSE - 1,    ## Addtional doseing 
                      cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                      replicate = FALSE)    ## No replicate
    
    
    ex_1 <- Gex.oral_1
    ex_2 <- Gex.oral_2
    
    ## set up time grid
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD0 to GD 22 
    
    # Combine data and run the simulation
    Gout_1 <- Gmod_R %>%data_set(ex_1) %>%
              idata_set(idata) %>% 
              update(atol = 1e-6, maxstep = 50000) %>%
              mrgsim(obsonly=TRUE, tgrid = Gtsamp)   
    
    Goutdf_1 = cbind.data.frame ( ID     = Gout_1$ID,
                                  Time   = Gout_1$time/24, 
                                  CPlas  = Gout_1$Plasma, 
                                  CL     = Gout_1$Liver,
                                  CFL    = Gout_1$Liver_Fet)
    Gout_2 <- 
        Gmod_R %>%
        param (pars) %>%
        update(atol = 1E-6, maxsteps = 50000) %>%          
        mrgsim_d (data = ex_2, tgrid = Gtsamp)
    
    Goutdf_2 = cbind.data.frame(Time    = Gout_2$time/24, 
                              CPlas     = Gout_2$Plasma, 
                              CL        = Gout_2$Liver,
                              CPlas_pup = Gout_2$Plasma_Fet,
                              CFL    = Gout_2$Liver_Fet)
    
    return (list(Goutdf_1, Goutdf_2)) # Return Goutdf
}

####################################
R_Gpars <- GFit_R$par
N = 1000

PlotDat_A1   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_1", Tissue = "Plasma")
PlotDat_A2   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_2", Tissue = "Plasma")
PlotDat_A3   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_3", Tissue = "Plasma")
PlotDat_A4   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_4", Tissue = "Plasma")
PlotDat_A5   <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_5", Tissue = "Plasma")
PlotDat_A6   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_1", Tissue = "MLiver")
PlotDat_A7   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_2", Tissue = "MLiver")
PlotDat_A8   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_3", Tissue = "MLiver")
PlotDat_A9   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_4", Tissue = "MLiver")
PlotDat_A10   <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_5", Tissue = "MLiver")
PlotDat_A11   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_1", Tissue = "FLiver")
PlotDat_A12   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_2", Tissue = "FLiver")
PlotDat_A13   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_3", Tissue = "FLiver")
PlotDat_A14   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_4", Tissue = "FLiver")
PlotDat_A15   <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_5", Tissue = "FLiver")


#########################################################
PlotDat_A1_m   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_1", Tissue = "Plasma")
PlotDat_A2_m   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_2", Tissue = "Plasma")
PlotDat_A3_m   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_3", Tissue = "Plasma")
PlotDat_A4_m   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_4", Tissue = "Plasma")
PlotDat_A5_m   <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "A_5", Tissue = "Plasma")
PlotDat_A6_m   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_1", Tissue = "MLiver")
PlotDat_A7_m   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_2", Tissue = "MLiver")
PlotDat_A8_m   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_3", Tissue = "MLiver")
PlotDat_A9_m   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_4", Tissue = "MLiver")
PlotDat_A10_m  <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL)  %>% mutate(DOSE_Group = "A_5", Tissue = "MLiver")
PlotDat_A11_m   <- Pop_Rat_G_A (R_Gpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_1", Tissue = "FLiver")
PlotDat_A12_m   <- Pop_Rat_G_A (R_Gpars, 2, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_2", Tissue = "FLiver")
PlotDat_A13_m   <- Pop_Rat_G_A (R_Gpars, 3, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_3", Tissue = "FLiver")
PlotDat_A14_m   <- Pop_Rat_G_A (R_Gpars, 5, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_4", Tissue = "FLiver")
PlotDat_A15_m  <- Pop_Rat_G_A (R_Gpars, 10,N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CFL)  %>% mutate(DOSE_Group = "A_5", Tissue = "FLiver")

########################################################
OBS_A1 <- OBS.A1 %>% mutate(DOSE_Group = "A_1", Tissue = "Plasma") %>% rename(Conc = CPlas)
OBS_A2 <- OBS.A2 %>% mutate(DOSE_Group = "A_2", Tissue = "Plasma") %>% rename(Conc = CPlas)
OBS_A3 <- OBS.A3 %>% mutate(DOSE_Group = "A_3", Tissue = "Plasma") %>% rename(Conc = CPlas)
OBS_A4 <- OBS.A4 %>% mutate(DOSE_Group = "A_4", Tissue = "Plasma") %>% rename(Conc = CPlas)
OBS_A5 <- OBS.A5 %>% mutate(DOSE_Group = "A_5", Tissue = "Plasma") %>% rename(Conc = CPlas) 
OBS_A6 <- OBS.A6 %>% mutate(DOSE_Group = "A_1", Tissue = "Mliver") %>% rename(Conc = CL)
OBS_A7 <- OBS.A7 %>% mutate(DOSE_Group = "A_2", Tissue = "Mliver") %>% rename(Conc = CL)
OBS_A8 <- OBS.A8 %>% mutate(DOSE_Group = "A_3", Tissue = "Mliver") %>% rename(Conc = CL)
OBS_A9 <- OBS.A9 %>% mutate(DOSE_Group = "A_4", Tissue = "Mliver") %>% rename(Conc = CL)
OBS_A10 <- OBS.A10 %>% mutate(DOSE_Group = "A_5", Tissue = "Mliver") %>% rename(Conc = CL) 
OBS_A11 <- OBS.A11 %>% mutate(DOSE_Group = "A_1", Tissue = "Fliver") %>% rename(Conc = CFL)
OBS_A12 <- OBS.A12 %>% mutate(DOSE_Group = "A_2", Tissue = "Fliver") %>% rename(Conc = CFL)
OBS_A13 <- OBS.A13 %>% mutate(DOSE_Group = "A_3", Tissue = "Fliver") %>% rename(Conc = CFL)
OBS_A14 <- OBS.A14 %>% mutate(DOSE_Group = "A_4", Tissue = "Fliver") %>% rename(Conc = CFL)
OBS_A15 <- OBS.A15 %>% mutate(DOSE_Group = "A_5", Tissue = "Fliver") %>% rename(Conc = CFL) 



PlotDat_Plas  <- rbind.data.frame (PlotDat_A1, PlotDat_A2, PlotDat_A3, PlotDat_A4, PlotDat_A5)
PlotDat_Liv   <- rbind.data.frame (PlotDat_A6, PlotDat_A7, PlotDat_A8, PlotDat_A9, PlotDat_A10)
PlotDat_FLiv   <- rbind.data.frame (PlotDat_A11, PlotDat_A12, PlotDat_A13, PlotDat_A14, PlotDat_A15)

PlotDat_MPlas <- rbind.data.frame (PlotDat_A1_m, PlotDat_A2_m, PlotDat_A3_m, PlotDat_A4_m, PlotDat_A5_m)
PlotDat_MLiv  <- rbind.data.frame (PlotDat_A6_m, PlotDat_A7_m, PlotDat_A8_m, PlotDat_A9_m, PlotDat_A10_m)
PlotDat_MFLiv  <- rbind.data.frame (PlotDat_A11_m, PlotDat_A12_m, PlotDat_A13_m, PlotDat_A14_m, PlotDat_A15_m)

OBS_Plas <- rbind.data.frame(OBS_A1,OBS_A2,OBS_A3,OBS_A4,OBS_A5)
OBS_Liv  <- rbind.data.frame(OBS_A6,OBS_A7,OBS_A8,OBS_A9,OBS_A10)
OBS_FLiv  <- rbind.data.frame(OBS_A11,OBS_A12,OBS_A13,OBS_A14,OBS_A15)

labe_A <- c("Dose: 1 mg/kg/day", "Dose: 2 mg/kg/day", "Dose: 3 mg/kg/day", "Dose: 5 mg/kg/day", "Dose: 10 mg/kg/day")
names(labe_A)  <- c("A_1", "A_2", "A_3", "A_4", "A_5")

plot.Plas <- 
    ggplot() + 
    geom_line  (data = PlotDat_Plas, aes(x = Time, y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MPlas, aes(x = Time, y = Conc, group = factor(DOSE_Group)), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_Plas, aes(x = Time/24, y = Conc), colour = "lightpink4", size = 3) +
    facet_wrap(~DOSE_Group, ncol = 5, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  

plot.Liv <- 
    ggplot() + 
    geom_line  (data = PlotDat_Liv, aes(x = Time, y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MLiv, aes(x = Time, y = Conc), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_Liv, aes(x = Time/24, y = Conc), colour = "lightpink4", size = 3) +
    facet_wrap(~DOSE_Group, ncol = 5, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  

plot.FLiv <- 
    ggplot() + 
    geom_line  (data = PlotDat_FLiv, aes(x = Time, y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MFLiv, aes(x = Time, y = Conc), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_FLiv, aes(x = Time/24, y = Conc), colour = "lightpink4", size = 3) +
    scale_x_continuous(expand = c (0,0), limits = c(17, 21.5)) +
    facet_wrap(~DOSE_Group, ncol = 5, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  
   

PlotTheme <- theme_bw() + theme (axis.text = element_text(size = rel(1.2), face = "bold",colour = "black"),
             strip.text = element_text(size = rel(1.2),  face = "bold", colour = "grey10"))  
             
 
plot.Plas <- plot.Plas + PlotTheme + xlab("") + ylab("") 
plot.Liv  <- plot.Liv + PlotTheme  + xlab("") + ylab("")
plot.FLiv <- plot.FLiv + PlotTheme  + xlab("") + ylab("")

grid.arrange (plot.Plas, plot.Liv, plot.FLiv, nrow = 3)


ggsave ("Fig.S1.tiff",scale = 1,
       plot = grid.arrange (plot.Plas, plot.Liv, plot.FLiv, nrow = 3),
       #path = "C:/Users/weichunc/Desktop",
       width = 28, height = 28, units = "cm", dpi = 320)

#==================================================================================================================== 
# Model calibration for lactational PBPK model based on the data of Chang et al., 2009          
# Exposure scenario: Dose regimen: 0.1, 0.3, 1 mg/kg/day; exposure from GD0 (day positive for mating) through PND20.                                                                                 #
# Abreviation: D: Dam, P: Pup, MP: maternal plasma, ML: maternal liver, NP: neonatal plasma, NL: neonatal liver
# Note: The samples collected at GD20 were used for model evaluation; PND72 can't be used due to the mode limitation
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


## Exposure scenario B. Dosing from GD0 - PND21; the dams were sacrificed on GD 20; PND20 based on study design of Chang et al. (2009) 
Pop_Rat_G_B <- function(Gpars, Lpars, DOSE, N) {
    
    ## Get out of log domain
    Gpars <- exp(Gpars)        ## Return a list of exp (parametrs for gestational model) from log scale
    Lpars <- exp(Lpars)           ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.225                  ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from GD0 - GD19
    GDOSE        = DOSE                   ## Repeat oral dose of 1 mg/kg/day 
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Simulation time from GD2 to GD 22 with dosing of GTDOSE - 1
    
    ## Simulation of exposure scenaior A and B (oral single dose to 15 mg/kg)
    Gout <- 
        Gmod_R %>%
        param (Gpars) %>%
        update(atol = 1E-8, maxsteps = 50000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Init <- Gout %>% filter (time == 21*24) %>% select(-c("Plasma", "Liver", "Plasma_Fet"))
    
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.242                 ## Rat body weight on PND0 (GD21); from Thibodeaux et al., 2003
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from PND0 - PND21
    LDOSE        = DOSE                   ## Repeat oral dose of GDOSE mg/kg/day;  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose
    
    ## Amount of oral dose
    idata <- 
        tibble(ID = 1:N) %>% 
        mutate( Free        = rlnormTrunc (N, min = exp(-6.01), max = exp(-5.32), meanlog = -5.62, sdlog = 0.198),
                Free_pup    = rlnormTrunc (N, min = exp(-4.22), max = exp(-3.44), meanlog = -3.84, sdlog = 0.198),
                GFRC        = rnormTrunc  (N, min = 24.95, max = 57.13, mean = 41, sd = 8.208),
                KbileC      = rlnormTrunc (N, min = exp(-6.57), max = exp(-5.41), meanlog = -5.99, sdlog = 0.294),
                Km_apical_p = rlnormTrunc (N, min = exp(5.009), max = exp(6.16), meanlog = 5.58, sdlog = 0.29),
                KMilk0      = rlnormTrunc (N, min = exp(-1.68), max = exp(-0.90), meanlog = -1.29, sdlog = 0.13), 
                KurineC_pup = rlnormTrunc (N, min = exp(-0.148), max = exp(1.002), meanlog = 0.45, sdlog = 0.19),
                PL          = rlnormTrunc (N, min = exp(0.59), max = exp(1.37), meanlog = 0.98, sdlog = 0.198),
                PL_pup      = rlnormTrunc (N, min = exp(0.53), max = exp(1.30), meanlog = 0.92, sdlog = 0.198),
                PMilkM      = rlnormTrunc (N, min = exp(0.22), max = exp(1.17), meanlog = 0.60, sdlog = 0.19),
                PRest_pup   = rlnormTrunc (N, min = exp(-1.92), max = exp(-1.14), meanlog = -1.53, sdlog = 0.13),
                RAFapi      = rlnormTrunc (N, min = exp(0.80), max = exp(1.955), meanlog = 1.38, sdlog = 0.29),
                VFilC       = rnormTrunc  (N, min = 0.00035, max = 0.00133, mean = 0.00084, sd = 0.00025),
                Vmax_apical_invitro_p = rlnormTrunc (N, min = exp(5.35), max = exp(6.507), meanlog = 5.93, sdlog = 0.29),
                Vmax_apical_invitro = 4141,
                PRest       = 0.03,
                SA_BW       = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
                SA_VLC      = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
                SA_KMilkC   = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
                SA_VL_pup   = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
                DOSEoral = LDOSEoral)
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Lex.oral_1 <- ev (ID   = 1:N,            ## One individual
                      amt  = idata$DOSEoral,    ## Amount of dose 
                      ii   = tinterval,    ## Time interval
                      addl = LTDOSE - 1,   ## Addtional doseing 
                      cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                      replicate = FALSE)   ## No replicate
    
    Lex.oral_2 <- ev (ID   = 1,             ## One individual
                      amt  = LDOSEoral,     ## Amount of dose 
                      ii   = tinterval,     ## Time interval
                      addl = LTDOSE - 1,    ## Addtional doseing 
                      cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                      replicate = FALSE)    ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 1) ## Simulation time from PND0 to PND74 with dosing of LTDOSE - 1+ 2 days (No dosing)
    
    
    Lout_1 <- Lmod_R %>% idata_set(idata) %>% 
              init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
              AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
              Aefflux = Init$Aefflux, APlas_free_pup = Init$APlas_Fet_free, AL_pup = Init$AL_Fet, ARest_pup = Init$ARest_Fet) %>%
              update(atol = 1E-8, maxsteps = 500000) %>%          
              mrgsim (data = Lex.oral_1, tgrid = Ltsamp) 
    
    Lout_2 <- 
        Lmod_R %>%
        param (Lpars) %>% 
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux, APlas_free_pup = Init$APlas_Fet_free, AL_pup = Init$AL_Fet, ARest_pup = Init$ARest_Fet) %>%
        update(atol = 1E-8, maxsteps = 500000) %>%          
        mrgsim_d (data = Lex.oral_2, tgrid = Ltsamp)
    
    
    
    Loutdf_1 = cbind.data.frame(ID  = Lout_1$ID,
                                Time   = Lout_1$time/24 + 22, 
                                CPlas  = Lout_1$Plasma, 
                                CPlas_pup = Lout_1$Plasma_pup,
                                CL_pup    = Lout_1$Liver_pup)
    
    Loutdf_2 = cbind.data.frame(ID  = Lout_2$ID,
                                Time   = Lout_2$time/24 + 22, 
                                CPlas  = Lout_2$Plasma, 
                                CPlas_pup = Lout_2$Plasma_pup,
                                CL_pup    = Lout_2$Liver_pup)
    
    
    return (list(Loutdf_1, Loutdf_2))
}

####################################
R_Lpars <- LFit_R$par
N = 1000

PlotDat_B1   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B2   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B3   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_1")
PlotDat_B4   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B5   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B6   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_1")
PlotDat_B7   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B8   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B9   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[1]]  %>% filter(Time > 0) %>% select (ID = ID, Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_1")

#########################################################
PlotDat_B1_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B2_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B3_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas)  %>% mutate(DOSE_Group = "B_1")
PlotDat_B4_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B5_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[2]] %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B6_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CPlas_pup)  %>% mutate(DOSE_Group = "B_1")
PlotDat_B7_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_3")
PlotDat_B8_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.3, N = N)[[2]] %>% filter(Time > 0) %>% select (Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_2")
PlotDat_B9_m   <- Pop_Rat_G_B (R_Gpars, R_Lpars, 0.1, N = N)[[2]]  %>% filter(Time > 0) %>% select (Time = Time, Conc = CL_pup)  %>% mutate(DOSE_Group = "B_1")

########################################################
OBS_B1 <- OBS.B1 %>% mutate(DOSE_Group = "B_3") %>% rename(Conc = CPlas) 
OBS_B2 <- OBS.B2 %>% mutate(DOSE_Group = "B_2") %>% rename(Conc = CPlas)
OBS_B3 <- OBS.B3 %>% mutate(DOSE_Group = "B_1") %>% rename(Conc = CPlas)
OBS_B4 <- OBS.B4 %>% mutate(DOSE_Group = "B_3") %>% rename(Conc = CPlas_pup)%>% filter(Time <60*24)
OBS_B5 <- OBS.B5 %>% mutate(DOSE_Group = "B_2") %>% rename(Conc = CPlas_pup)%>% filter(Time <60*24)
OBS_B6 <- OBS.B6 %>% mutate(DOSE_Group = "B_1") %>% rename(Conc = CPlas_pup)%>% filter(Time <60*24)
OBS_B7 <- OBS.B7 %>% mutate(DOSE_Group = "B_3") %>% rename(Conc = CL_pup)%>% filter(Time <60*24)
OBS_B8 <- OBS.B8 %>% mutate(DOSE_Group = "B_2") %>% rename(Conc = CL_pup)%>% filter(Time <60*24)
OBS_B9 <- OBS.B9 %>% mutate(DOSE_Group = "B_1") %>% rename(Conc = CL_pup)%>% filter(Time <60*24)

PlotDat_Plas_L  <- rbind.data.frame (PlotDat_B1, PlotDat_B2, PlotDat_B3)
PlotDat_FPlas_L <- rbind.data.frame (PlotDat_B4, PlotDat_B5, PlotDat_B6)
PlotDat_FLiv_L  <- rbind.data.frame (PlotDat_B7, PlotDat_B8, PlotDat_B9)

PlotDat_MPlas_L <- rbind.data.frame (PlotDat_B1_m, PlotDat_B2_m, PlotDat_B3_m)
PlotDat_MFPlas_L  <- rbind.data.frame (PlotDat_B4_m, PlotDat_B5_m, PlotDat_B6_m)
PlotDat_MFLiv_L  <- rbind.data.frame (PlotDat_B7_m, PlotDat_B8_m, PlotDat_B9_m)

OBS_Plas_L <- rbind.data.frame(OBS_B1,OBS_B2,OBS_B3)
OBS_FPlas_L  <- rbind.data.frame(OBS_B4,OBS_B5,OBS_B6)
OBS_FLiv_L  <- rbind.data.frame(OBS_B7,OBS_B8,OBS_B9)

labe_A <- c("Dose: 0.1 mg/kg/day", "Dose: 0.3 mg/kg/day", "Dose: 1 mg/kg/day")
names(labe_A)  <- c("B_1", "B_2", "B_3")


PlotDat_Plas_L <- PlotDat_Plas_L #%>% filter (Conc < 200 & Conc >0)

plot.Plas_L <- 
    ggplot() + 
    geom_line  (data = PlotDat_Plas_L, aes(x = Time - 21 , y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MPlas_L, aes(x = Time - 21, y = Conc, group = factor(DOSE_Group)), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_Plas_L, aes(x = (Time/24)-21, y = Conc), colour = "lightpink4", size = 3) +
    facet_wrap(~DOSE_Group, ncol = 3, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  

plot.FPlas_L <- 
    ggplot() + 
    geom_line  (data = PlotDat_FPlas_L, aes(x = Time - 21, y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MFPlas_L, aes(x = Time - 21, y = Conc, group = factor(DOSE_Group)), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_FPlas_L, aes(x = (Time/24)-21, y = Conc), colour = "lightpink4", size = 3) +
    facet_wrap(~DOSE_Group, ncol = 3, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  

plot.FLiv_L <- 
    ggplot() + 
    geom_line  (data = PlotDat_FLiv_L, aes(x = Time - 21, y = Conc, group = factor(ID)), colour = "lightpink3", alpha = 0.2) +
    geom_line  (data = PlotDat_MFLiv_L, aes(x = Time - 21, y = Conc, group = factor(DOSE_Group)), colour = "lightpink4", size = 0.8) +
    geom_point (data = OBS_FLiv_L, aes(x = (Time/24)-21, y = Conc), colour = "lightpink4", size = 3) +
    facet_wrap(~DOSE_Group, ncol = 3, scales = "fixed",labeller = labeller(DOSE_Group = labe_A))  


PlotTheme <- theme_bw() + theme (axis.text = element_text(size = rel(1.2), face = "bold",colour = "black"),
                                 strip.text = element_text(size = rel(1.2),  face = "bold", colour = "grey10"))  


plot.Plas_L <- plot.Plas_L + PlotTheme + xlab("") + ylab("") 
plot.FPlas_L  <- plot.FPlas_L + PlotTheme  + xlab("") + ylab("")
plot.FLiv_L <- plot.FLiv_L + PlotTheme  + xlab("") + ylab("")


  
#grid.arrange (plot.Plas_L, plot.FPlas_L, plot.FLiv_L, ncol = 1)


ggsave ("Fig.S2.tiff",scale = 1,
        plot = grid.arrange (plot.Plas_L, plot.FPlas_L, plot.FLiv_L, ncol = 1),
        width = 28, height = 28, units = "cm", dpi = 320)


