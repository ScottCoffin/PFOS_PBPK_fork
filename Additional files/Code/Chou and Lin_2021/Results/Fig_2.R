#===========================================================================================================================
# The code for the Fig. 2: Comparisons of model calibration between predictions (y-axis) with experimental data (x-axis)
# - Author : Wei-Chun Chou
# - Date   : Feb, 2020
# - Fig.2A : Global evaluation of goodness of model fit 
# - Fig.2B : Predicted-to-observed ratio versus model prediction plot. 
#===========================================================================================================================

##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(gridExtra)   ## R-package for function "grid.arrange"   

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")
source (file = "HMod.R")

## Compile mrgsolve-based PBPK Model
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

## Input defined R files
source (file = "Pred_Rat.R")
source (file = "Pred_H.R")

## Read the data sets
#===========================================================================================================================
# Model calibration data A 
# Reference: Thibodeaux et al., 2003; 
# Dose regimen: 1, 2, 3, 5, 10 mg/kg/day, exposure from GD2 - GD20
# Abreviation: maternal plasma (MP), maternal liver (ML), fetal liver (FL)
# 
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
#
# Model calibration data B for rat:  Chang et al., 2009; 
# Dose regimen: 0.1, 0.3, 1 mg/kg/day; exposure from GD0 (day positive for mating) through PND20.                                                                                 #
# Abreviation: D: Dam, P: Pup, MP: maternal plasma, ML: maternal liver, NP: neonatal plasma, NL: neonatal liver
# Note: The samples collected at GD20 were used for model evaluation; PND72 can't be used due to the mode limitation
# 
# B1. : Pregnant SD rat oral dose to 0.1 mg/kg,  matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                                            
# B2. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                           
# B3. : Pregnant SD rat oral dose to 1 mg/kg,    matrix: MP, Sampling time: GD20, PND4, PND21 and PND72                            
# B4. : Pregnant SD rat oral dose to 0.1 mg/kg,  matrix: NP, Sampling time: GD20, PND4, PND21 and PND72                                             
# B5. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: NP, Sampling time: GD20, PND4, PND21 and PND72                           
# B6. : Pregnant SD rat oral dose to 1 mg/kg,    matrix: NP, Sampling time: GD20, PND4, PND21 and PND72  
# B7. : Pregnant SD rat oral dose to 0.1 mg/kg,  matrix: NL, Sampling time: GD20, PND4, PND21 and PND72                                             
# B8. : Pregnant SD rat oral dose to 0.3 mg/kg,  matrix: NL, Sampling time: GD20, PND4, PND21 and PND72                           
# B9. : Pregnant SD rat oral dose to 1 mg/kg,    matrix: NL, Sampling time: GD20, PND4, PND21 and PND72  
#==============================================================================================================================================================                                                   

## Read the data and later used in model calibration and evluation
OBS_R_A1  <- RDat %>% filter(Study == 1 & Sample == "MP" & Dose == 1   & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS_R_A2  <- RDat %>% filter(Study == 1 & Sample == "MP" & Dose == 2   & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS_R_A3  <- RDat %>% filter(Study == 1 & Sample == "MP" & Dose == 3   & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS_R_A4  <- RDat %>% filter(Study == 1 & Sample == "MP" & Dose == 5   & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS_R_A5  <- RDat %>% filter(Study == 1 & Sample == "MP" & Dose == 10  & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS_R_A6  <- RDat %>% filter(Study == 1 & Sample == "ML" & Dose == 1   & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS_R_A7  <- RDat %>% filter(Study == 1 & Sample == "ML" & Dose == 2   & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS_R_A8  <- RDat %>% filter(Study == 1 & Sample == "ML" & Dose == 3   & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS_R_A9  <- RDat %>% filter(Study == 1 & Sample == "ML" & Dose == 5   & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS_R_A10 <- RDat %>% filter(Study == 1 & Sample == "ML" & Dose == 10  & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS_R_A11 <- RDat %>% filter(Study == 1 & Sample == "FL" & Dose == 1   & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS_R_A12 <- RDat %>% filter(Study == 1 & Sample == "FL" & Dose == 2   & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS_R_A13 <- RDat %>% filter(Study == 1 & Sample == "FL" & Dose == 3   & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS_R_A14 <- RDat %>% filter(Study == 1 & Sample == "FL" & Dose == 5   & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS_R_A15 <- RDat %>% filter(Study == 1 & Sample == "FL" & Dose == 10  & Time != 0) %>% select(Time = "Time", CFL   = "Conc")
OBS_R_B1  <- RDat %>% filter(Study == 2 & Sample == "MP" & Dose == 1   & Time > 480)%>% select(Time = "Time", CPlas = "Conc")
OBS_R_B2  <- RDat %>% filter(Study == 2 & Sample == "MP" & Dose == 0.3 & Time > 480)%>% select(Time = "Time", CPlas = "Conc")
OBS_R_B3  <- RDat %>% filter(Study == 2 & Sample == "MP" & Dose == 0.1 & Time > 480)%>% select(Time = "Time", CPlas = "Conc")
OBS_R_B4  <- RDat %>% filter(Study == 2 & Sample == "NP" & Dose == 1   & Time > 480)%>% select(Time = "Time", CPlas_pup = "Conc")
OBS_R_B5  <- RDat %>% filter(Study == 2 & Sample == "NP" & Dose == 0.3 & Time > 480)%>% select(Time = "Time", CPlas_pup = "Conc")
OBS_R_B6  <- RDat %>% filter(Study == 2 & Sample == "NP" & Dose == 0.1 & Time > 480)%>% select(Time = "Time", CPlas_pup = "Conc")
OBS_R_B7  <- RDat %>% filter(Study == 2 & Sample == "NL" & Dose == 1   & Time > 480)%>% select(Time = "Time", CL_pup = "Conc")
OBS_R_B8  <- RDat %>% filter(Study == 2 & Sample == "NL" & Dose == 0.3 & Time > 480)%>% select(Time = "Time", CL_pup = "Conc")
OBS_R_B9  <- RDat %>% filter(Study == 2 & Sample == "NL" & Dose == 0.1 & Time > 480)%>% select(Time = "Time", CL_pup = "Conc")

#=====================================================================================================================================
# Model calibration data A for human
# Abbrevation: MP: Maternal plasma, CB: cord blood, Pla: placenta; Fli: Fetal liver
#
# A1. (Inoue et al., 2004)         : Japan population;   Matrix: MP, CB        
# A2. (Fei et al., 2007)           : Danish population;  Matrix: MP, CB           
# A3. (Kato et al., 2014)          : U.S. population;    Matrix: MP, CB          
# A4. (Pan et al., 2017)           : Chinese population; Matrix: MP, CB         
# A5. (Mamsen et al., 2019)        : Sweden population;  Matrix: MP, Pla                
# B1. (Karrman et al. (2007)       : Sweden population;  Matrix: MP, Milk                                 
# B2. (Von Ehrenstein et al., 2009): U.S. population;    Matrix: MP                       
# B3. (Fromme et al., 2010)        : Gemerny population; Matrix: MP, NP     
#========================================================================================================================================

## Read the data and later used in model calibration and evluation
OBS_H_A1_MP     <- HDat %>% filter(Study == 1  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A1_CB     <- HDat %>% filter(Study == 1  & Matrix == "CB")   %>% select(Time = "Time",  CFPlas    = "Conc")
OBS_H_A2_MP     <- HDat %>% filter(Study == 2  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A2_CB     <- HDat %>% filter(Study == 2  & Matrix == "CB")   %>% select(Time = "Time",  CFPlas    = "Conc")
OBS_H_A3_MP     <- HDat %>% filter(Study == 3  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A3_CB     <- HDat %>% filter(Study == 3  & Matrix == "CB")   %>% select(Time = "Time",  CFPlas    = "Conc")
OBS_H_A4_MP     <- HDat %>% filter(Study == 4  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A4_CB     <- HDat %>% filter(Study == 4  & Matrix == "CB")   %>% select(Time = "Time",  CFPlas    = "Conc")
OBS_H_A5_MP     <- HDat %>% filter(Study == 5  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A5_Pla    <- HDat %>% filter(Study == 5  & Matrix == "Pla")  %>% select(Time = "Time",  CPla      = "Conc")
OBS_H_A5_Fli    <- HDat %>% filter(Study == 5  & Matrix == "Fli")  %>% select(Time = "Time",  CFL       = "Conc")
OBS_H_A6_MP     <- HDat %>% filter(Study == 6  & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_A6_CB     <- HDat %>% filter(Study == 6  & Matrix == "CB")   %>% select(Time = "Time",  CFPlas    = "Conc")
OBS_H_B1_MP     <- HDat %>% filter(Study == 16 & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_B1_Milk   <- HDat %>% filter(Study == 16 & Matrix == "Milk") %>% select(Time = "Time",  CMilk     = "Conc")
OBS_H_B2_MP     <- HDat %>% filter(Study == 17 & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_B3_MP     <- HDat %>% filter(Study == 18 & Matrix == "MP")   %>% select(Time = "Time",  CPlas     = "Conc")
OBS_H_B3_NP     <- HDat %>% filter(Study == 18 & Matrix == "NP")   %>% select(Time = "Time",  CPlas_pup = "Conc")
OBS_H_B4_Milk   <- HDat %>% filter(Study == 19 & Matrix == "Milk") %>% select(Time = "Time",  CMilk     = "Conc")

## Cost function for rat
R_Cost<-function (Gpars_R, Lpars_R){ # DOSE is based on the exposure scenario of literature (Thibodeaux et al., 2003 and Change et al., 2009)
    outdf.A1 <- pred.A (Gpars = Gpars_R,  DOSE = 1)
    outdf.A2 <- pred.A (Gpars = Gpars_R,  DOSE = 2)
    outdf.A3 <- pred.A (Gpars = Gpars_R,  DOSE = 3)
    outdf.A4 <- pred.A (Gpars = Gpars_R,  DOSE = 5)
    outdf.A5 <- pred.A (Gpars = Gpars_R,  DOSE = 10)
    outdf.A  <- pred.B (Lpars = Lpars_R,  DOSE = 1)
    outdf.B  <- pred.B (Lpars = Lpars_R,  DOSE = 0.3)
    outdf.C  <- pred.B (Lpars = Lpars_R,  DOSE = 0.1)
    
    cost<- modCost  (model = outdf.A1, obs = OBS_R_A1, x ="Time")
    cost<- modCost  (model = outdf.A2, obs = OBS_R_A2, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS_R_A3, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_R_A4, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS_R_A5, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A1, obs = OBS_R_A6, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS_R_A7, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS_R_A8, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_R_A9, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS_R_A10,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A1, obs = OBS_R_A11,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS_R_A12,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS_R_A13,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_R_A14,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A5, obs = OBS_R_A15,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A,  obs = OBS_R_B1, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B,  obs = OBS_R_B2, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C,  obs = OBS_R_B3, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A,  obs = OBS_R_B4, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B,  obs = OBS_R_B5, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C,  obs = OBS_R_B6, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A,  obs = OBS_R_B7, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B,  obs = OBS_R_B8, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.C,  obs = OBS_R_B9, x ="Time", cost = cost)
    return(cost)
}

## Cost function for human
H_Cost<-function (Gpars_H, Lpars_H){ ## DOSE_A1 - B4 were defined in the prediction function "Pred_H.R"
    outdf.A1 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A1, Init = Init_A1)
    outdf.A2 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A2, Init = Init_A2)
    outdf.A3 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A3, Init = Init_A3)
    outdf.A4 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A4, Init = Init_A4)
    outdf.A5 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A5, Init = Init_A5)
    outdf.A6 <- pred.G (Gpars = Gpars_H,  DOSE = DOSE_A6, Init = Init_A6)
    outdf.B1 <- pred.L (Lpars = Lpars_H,  DOSE = DOSE_B1, Init_G = Init_G1)
    outdf.B2 <- pred.L (Lpars = Lpars_H,  DOSE = DOSE_B2, Init_G = Init_G2)
    outdf.B3 <- pred.L (Lpars = Lpars_H,  DOSE = DOSE_B3, Init_G = Init_G3)
    outdf.B4 <- pred.L (Lpars = Lpars_H,  DOSE = DOSE_B4, Init_G = Init_G4)
    
    cost<- modCost  (model = outdf.A1, obs = OBS_H_A1_MP,  x ="Time")
    cost<- modCost  (model = outdf.A1, obs = OBS_H_A1_CB,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS_H_A2_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A2, obs = OBS_H_A2_CB,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS_H_A3_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A3, obs = OBS_H_A3_CB,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_H_A4_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_H_A4_CB,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_H_A5_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_H_A5_Pla, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A4, obs = OBS_H_A5_Fli, x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A6, obs = OBS_H_A6_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.A6, obs = OBS_H_A6_CB,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B1, obs = OBS_H_B1_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B1, obs = OBS_H_B1_Milk,x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B2, obs = OBS_H_B2_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B3, obs = OBS_H_B3_MP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B3, obs = OBS_H_B3_NP,  x ="Time", cost = cost)
    cost<- modCost  (model = outdf.B4, obs = OBS_H_B4_Milk,x ="Time", cost = cost)
    
    
    return(cost)
}

## Input the optimized parameters value into the cost function
RatCost   <- R_Cost(GFit_R$par, LFit_R$par) 
HumanCost <- H_Cost(GFit_H$par, LFit_H$par)


############## Plot code for Figure 2a ################
GPDat <- cbind.data.frame (OBS = RatCost$residuals$obs, 
                           PRE = RatCost$residuals$mod)

HPDat <- cbind.data.frame (OBS = HumanCost$residuals$obs, 
                           PRE = HumanCost$residuals$mod)

GPDat %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Rat")
HPDat %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Human")
PlotDat <- rbind(GPDat, HPDat)


## Estimating the R squared and goodness of fit using linear regression model
fit <- lm(Log.OBS ~ Log.PRE, data = PlotDat)
summary(fit)

## Make the plot data using the fitting results 
PlotDat %<>% mutate(res = residuals(fit), 
                    prediction = predict(fit), 
                    OPR = PRE/OBS) ## OPR: the ratio of prediction value and observed data

p <- 
    ggplot(PlotDat, aes(Log.OBS, Log.PRE)) + ## using log-sacle axis
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,3), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-2,3),labels = scales::math_format(10^.x))

## Set up your theme and font
windowsFonts("Times" = windowsFont("Times New Roman"))

p1 <- p + 
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, size=2),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(), 
        axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position         ='none') +
    labs (x = "",  y = "")

p1

################## Plot code for Figure 2b ####################


p2 <-
    ggplot(PlotDat, aes(Log.PRE, log(OPR,10))) +
    geom_hline(yintercept = log10(2),linetype = 3, color   = "black", size =1) +
    geom_hline(yintercept = log10(0.5),linetype = 3, color   = "black", size =1) +
    geom_point(color   = "steelblue4", 
               aes(shape= as.factor(Species)),size = 3) +
    geom_smooth(color="steelblue4",method="loess", se = FALSE, size = 1) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,3), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-2,3),labels = scales::math_format(10^.x))

p2 <- p2 +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, size=2),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(), 
        axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position='none') +
    labs (x = "", y = "")



grid.arrange(p1, p2, nrow = 1) # Combine two plot 


## Save the plots as tiff files
ggsave("Fig.2.tiff",scale = 1,
       plot = grid.arrange(p1, p2, nrow = 1),
       path = "C:/Users/weichunc/Desktop",
       width = 25, height = 12, units = "cm", dpi=320)


dev.off()




