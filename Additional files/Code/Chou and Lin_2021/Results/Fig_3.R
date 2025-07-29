#===========================================================================================================================
# Figure 3: Example of illustration of the general features/trends and model evaluation for PFOS levels in maternal plasma 
# - Author : Wei-Chun Chou
# - Date   : March, 2020
# - Fig.3A : The PFOS time-crouse profiles in maternal palsma during gestation  
# - Fig.3B : The PFOS time-crouse profiles in maternal palsma during lactation
# - Fig.3C : The PFOS time-crouse profiles in pup palsma during lactation 
#===========================================================================================================================

##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(ggExtra)     ## R-package for function "grid.arrange" 

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

# ================================================================================================== 
# Model evaluation for lactational PBPK model based on the data of Luebker et al., 2005a,b         
# Exposure scenario: 42 days prior to matting, through GD0 to GD20 or PND4
# Abreviation: D: Dam, P: Pup, MP: maternal plasma, NP: neonatal plasma, 
#                               
# C1. : Pregnant SD rat oral single dose to 0.1 mg/kg,  matrix: plasma                                            
# C2. : Pregnant SD rat oral single dose to 0.4 mg/kg,  matrix: plasma                          
# C3. : Pregnant SD rat oral single dose to 1.6 mg/kg,  matrix: plasma                          
# C4. : Pregnant SD rat oral single dose to 3.2 mg/kg,  matrix: plasma                                                                           #
#=====================================================================================================

## Read the data and later used in model calibration and evluation
OBS.C1_MP <- RDat %>% filter(Study >2) %>% filter( Sample == "MP" & Dose == 0.1 )%>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.C2_MP <- RDat %>% filter(Study >2) %>% filter( Sample == "MP" & Dose == 0.4) %>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.C3_MP <- RDat %>% filter(Study >2) %>% filter( Sample == "MP" & Dose == 1.6) %>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.C4_MP <- RDat %>% filter(Study >2) %>% filter( Sample == "MP" & Dose == 3.2) %>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.C1_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 0.1) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")
OBS.C2_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 0.4) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")
OBS.C3_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 1.6) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")
OBS.C4_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 3.2) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")
OBS.C5_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 0.8) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")
OBS.C6_PP <- RDat %>% filter(Study >2) %>% filter( Sample == "NP" & Dose == 2.0) %>% select(Time = "Time", CPlas_pup = "Conc", SD = "SD")

## Defined the prediction function
pred.C <- function(Gpars, Lpars, DOSE) {
    
    PreG_BW          = 0.24                    ## Female rat body weight (premating); from Luebker et al., 2005b
    tinterval        = 24                     ## Time interval; 
    PreG_TDOSE       = 42 + 14                ## Total dosing/Dose times; Repeat oral dose from beginning 42 days prior to cohabitation; (14 days matting)
    PreG_DOSE        = DOSE                   ## Repeat oral dose (mg/kg/day);  
    PreG_DOSEoral    = PreG_DOSE*PreG_BW      ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    PreG_ex.oral <- ev (ID   = 1,             ## One individual
                        amt  = PreG_DOSEoral, ## Amount of dose 
                        ii   = tinterval,     ## Time interval
                        addl = PreG_TDOSE - 1,## Addtional doseing 
                        cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                        tinf = 0.01,          ## Infusion time;  
                        replicate = FALSE)    ## No replicate
    
    PreG_tsamp       = tgrid(0, tinterval*(PreG_TDOSE - 1) + 24*1, 1) ## Simulation time from 6 weeks prior to mating and during mating (maximum of 14 days)
    
    PreG_out <- PreGmod_R %>% 
        update(atol = 1E-6, maxsteps = 50000) %>%          
        mrgsim_d (data = PreG_ex.oral, tgrid = PreG_tsamp)
    
    PreG_Init <- PreG_out %>% filter (time == (56*24)) %>% select(-c("Plasma", "Liver","Kidney"))
    
    ############### Gestational model
    ## Get out of log domain
    Gpars <- exp(Gpars)                   ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 0.36                   ## Body weight during gestation from Luebker et al., 2005b
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 22                     ## Total dosing/Dose times; Repeat oral dose from GD0 - GD21
    GDOSE        = DOSE                   ## Repeat oral dose 
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSE every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    amt  = GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = GTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)    ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*0, 1) ## Exposure time from GD0 to GD 20 and plus 1 day without dosing
    
    ## Simulation 
    Gout <- 
        Gmod_R %>%
        init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
             AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
             Aefflux = PreG_Init$Aefflux) %>%
        param (Gpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Goutdf = cbind.data.frame(Time      = Gout$time, 
                              CPlas     = Gout$Plasma,
                              CL        = Gout$Liver,
                              CPlas_pup = Gout$Plasma_Fet,
                              CL_pup    = Gout$Liver_Fet) 
    
    GInit <- Gout %>% filter (time == 21*24) 
    
    ### 
    Lpars <- exp(Lpars)          ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.30                  ## Rat body weight during lactation from Luebker et al., 2005b
    tinterval    = 24                    ## Time interval; 
    LTDOSE       = 23                    ## Total dosing/Dose times; Repeat oral dose from PND0 - PND22
    LDOSE        = DOSE                  ## Repeat oral dose of GDOSE mg/kg/day;  
    LDOSEoral    = LDOSE*LBW             ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    Lex.oral <- ev (ID   = 1,            ## One individual
                    amt  = LDOSEoral,    ## Amount of dose 
                    ii   = tinterval,    ## Time interval
                    addl = LTDOSE - 1,   ## Addtional doseing 
                    cmt  = "AST",        ## The dosing comaprtment: AST Stomach  
                    tinf = 0.01,          ## Infusion time;  
                    replicate = FALSE)   ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*2, 1) ## Simulation time from PND0 to PND4 plus 1 day without dosing (PND5)
    
    
    Lout <- Lmod_R %>% 
        init(APlas_free = GInit$APlas_free, APTC = GInit$APTC, AFil = GInit$AFil, AKb = GInit$AKb, ARest = GInit$ARest,
             AL = GInit$AL, AM = GInit$AM, AF = GInit$AF, A_baso = GInit$A_baso, A_apical = GInit$A_apical, Adif = GInit$Adif,
             Aefflux = GInit$Aefflux, APlas_free_pup = GInit$APlas_Fet_free, ARest_pup = GInit$ARest_Fet, AL_pup  = GInit$AL_Fet) %>%
        param (Lpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Lex.oral, tgrid = Ltsamp)
    
    Loutdf = cbind.data.frame(Time   = Lout$time, 
                              CPlas  = Lout$Plasma, 
                              CL     = Lout$Liver, 
                              CPlas_pup = Lout$Plasma_pup,
                              CL_pup = Lout$Liver_pup)
    
    return (list(Goutdf = Goutdf, Loutdf = Loutdf))
}

## Create the dataset for plot
out_0.1_G <- pred.C(GFit_R$par, LFit_R$par, 0.1)[[1]] %>% mutate(Dose_group = 0.1, Stage = "G", GD = Time/24)%>% filter(Time > 0)
out_0.1_L <- pred.C(GFit_R$par, LFit_R$par, 0.1)[[2]] %>% mutate(Dose_group = 0.1, Stage = "L", GD = Time/24)%>% filter(Time > 0)
out_0.1   <- rbind.data.frame (out_0.1_G, out_0.1_L)

out_0.4_G <- pred.C(GFit_R$par, LFit_R$par, 0.4)[[1]] %>% mutate(Dose_group = 0.4, Stage = "G", GD = Time/24)%>% filter(Time > 0)
out_0.4_L <- pred.C(GFit_R$par, LFit_R$par, 0.4)[[2]] %>% mutate(Dose_group = 0.4, Stage = "L", GD = Time/24)%>% filter(Time > 0)
out_0.4   <- rbind.data.frame (out_0.4_G, out_0.4_L)

out_0.8_G <- pred.C(GFit_R$par, LFit_R$par, 0.8)[[1]] %>% mutate(Dose_group = 0.8, Stage = "G", GD = Time/24)%>% filter(Time > 0)
out_0.8_L <- pred.C(GFit_R$par, LFit_R$par, 0.8)[[2]] %>% mutate(Dose_group = 0.8, Stage = "L", GD = Time/24)%>% filter(Time > 0)
out_0.8   <- rbind.data.frame (out_0.8_G, out_0.8_L)

out_1.6_G <- pred.C(GFit_R$par, LFit_R$par, 1.6)[[1]] %>% mutate(Dose_group = 1.6, Stage = "G", GD = Time/24)%>% filter(Time > 0)
out_1.6_L <- pred.C(GFit_R$par, LFit_R$par, 1.6)[[2]] %>% mutate(Dose_group = 1.6, Stage = "L", GD = Time/24)%>% filter(Time > 0)
out_1.6   <- rbind.data.frame (out_1.6_G, out_1.6_L)

out_2.0_G <- pred.C(GFit_R$par, LFit_R$par, 2)[[1]] %>% mutate(Dose_group = 2, Stage = "G", GD = Time/24)%>% filter(Time > 0)
out_2.0_L <- pred.C(GFit_R$par, LFit_R$par, 2)[[2]] %>% mutate(Dose_group = 2, Stage = "L", GD = Time/24)%>% filter(Time > 0)
out_2.0   <- rbind.data.frame (out_2.0_G, out_2.0_L)

### Input the experimental data for maternal plasma
obs_0.1   <- OBS.C1_MP %>% mutate(Dose_group = 0.1, Stage = if_else (Time <= 504, true = "G", false = "L")) %>% 
             mutate(GD = if_else (Stage == "L", true = Time/24 - 22, false = Time/24))
obs_0.4   <- OBS.C2_MP %>% mutate(Dose_group = 0.4, Stage = if_else (Time <= 504, true = "G", false = "L")) %>% 
             mutate(GD = if_else (Stage == "L", true = Time/24 - 22, false = Time/24))
obs_1.6   <- OBS.C3_MP %>% mutate(Dose_group = 1.6, Stage = if_else (Time <= 504, true = "G", false = "L")) %>% 
             mutate(GD = if_else (Stage == "L", true = Time/24 - 22, false = Time/24))
obs_3.2   <- OBS.C4_MP %>% mutate(Dose_group = 3.2, Stage = if_else (Time <= 504, true = "G", false = "L")) %>% 
             mutate(GD = if_else (Stage == "L", true = Time/24 - 22, false = Time/24))

### Input the experimental data for fetal plasma
obs_0.1_pp <- OBS.C1_PP %>% mutate (Dose_group = 0.1, GD = Time/24) 
obs_0.4_pp <- OBS.C2_PP %>% mutate (Dose_group = 0.4, GD = Time/24) 
obs_0.8_pp <- OBS.C5_PP %>% mutate (Dose_group = 0.8, GD = Time/24) 
obs_1.6_pp <- OBS.C3_PP %>% mutate (Dose_group = 1.6, GD = Time/24) 
obs_2.0_pp <- OBS.C6_PP %>% mutate (Dose_group = 2.0, GD = Time/24) 
obs_3.2_pp <- OBS.C4_PP %>% mutate (Dose_group = 3.2, GD = Time/24) 

## Cbomine all data from different dose groups
obs_pp     <- rbind.data.frame(obs_0.1_pp, obs_0.4_pp, obs_0.8_pp, obs_1.6_pp,
                               obs_2.0_pp, obs_3.2_pp)


## Add the font to font database
windowsFonts("Times" = windowsFont("Times New Roman"))

## Define the plot theme
PlotTheme <- theme (
  plot.background         = element_rect (fill="White"),
  text                    = element_text (family = "Times"),   # text front (Time new roman)
  panel.border            = element_rect (colour = "black", fill=NA, size = 2),
  panel.background        = element_rect (fill="White"),
  panel.grid.major        = element_blank(),
  panel.grid.minor        = element_blank(), 
  axis.text.x             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
  axis.text.y             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
  axis.title              = element_blank(),   
  strip.background        = element_blank(), 
  strip.text              = element_blank(),
  legend.position='none') 

###### Plot

## Fig 3
p1 <- ggplot(data = out_0.1_G , aes(x = GD, y = CPlas)) + 
      geom_line    (colour = "steelblue4", size = 1.5) + 
      geom_point   (data   = obs_0.1 %>%filter(Stage == "G"), aes (x = GD, y = CPlas), colour = "steelblue4", size = 2)+  
      geom_errorbar(data = obs_0.1 %>%filter(Stage == "G"), aes(x = GD, ymin = CPlas-SD, ymax = CPlas+SD), 
                    colour = "steelblue4", width=.2,
                    position=position_dodge(0.05)) +
      scale_x_continuous( expand = c (0, 0), limits = c(-1, 22)) +
      scale_y_continuous(limits = c(0, 20)) 



p2 <- ggplot(data = out_0.4_L , aes(x = GD, y = CPlas)) + 
      geom_line  (colour = "steelblue4", size = 1.5, linetype="twodash") + 
      geom_point (data = obs_0.4%>%filter(Stage == "L"), aes (x = GD, y = CPlas), colour = "steelblue4", size = 2)+  
      geom_errorbar(data = obs_0.4%>%filter(Stage == "L"), aes(x = GD, ymin = CPlas-SD, ymax = CPlas+SD), 
                colour = "steelblue4", width=.2,
                position=position_dodge(0.05)) +
      scale_x_continuous(expand = c (0,0), limits = c(0, 22)) + 
      scale_y_continuous(limits = c(5, 50)) 


p3 <- ggplot(data = out_0.8_L , aes(x = GD, y = CPlas_pup)) + 
      geom_line (colour = "steelblue4", size = 1.5, linetype="longdash") + 
      geom_point (data = obs_0.8_pp %>%filter(Time>504), aes (x = GD-22, y = CPlas_pup), colour = "steelblue4", size = 2) + 
      geom_errorbar(data = obs_0.8_pp%>%filter(Time>504), aes(x = GD-22, ymin = CPlas_pup-SD, ymax = CPlas_pup+SD), 
                    colour = "steelblue4", width=.2,
                    position=position_dodge(0.05)) +
      scale_x_continuous(expand = c (0,0), limits = c(0, 22)) +
      scale_y_continuous(limits = c(0, 150)) 

  


p1 <- p1 + PlotTheme
p2 <- p2 + PlotTheme
p3 <- p3 + PlotTheme


# ggsave("Fig.3a.tiff",scale = 1,
#        plot = p1,
#        path = "C:/Users/weichunc/Desktop",
#        width = 9, height = 8, units = "cm", dpi=320)
# 
# ggsave("Fig.3b.tiff",scale = 1,
#        plot = p2,
#        path = "C:/Users/weichunc/Desktop",
#        width = 9, height = 8, units = "cm", dpi=320)
# 
# ggsave("Fig.3c.tiff",scale = 1,
#        plot = p3,
#        path = "C:/Users/weichunc/Desktop",
#        width = 9, height = 8, units = "cm", dpi=320)


### Output the data ploting figure
# write.csv(out_0.1_G, file = 'Data_figure_3a.csv')
# write.csv(out_0.4_L, file = 'Data_figure_3b.csv')
# write.csv(out_0.8_L, file = 'Data_figure_3c.csv')


