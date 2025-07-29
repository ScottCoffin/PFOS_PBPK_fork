#===========================================================================================================================
# Figure 4: Example of illustration of the general features/trends 
# - Author : Wei-Chun Chou
# - Date   : March, 2020
# - Fig.4A : The PFOS time-crouse profiles in maternal palsma during gestation and lactation  
# - Fig.4B : The PFOS time-crouse profiles in human milk during lactation
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
source (file = "HMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code) 
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Loading datasets
RDat <- read.csv(file = "Data_R.csv") 
HDat <- read.csv(file = "Data_H.csv")

## Loading Fit results from rds files
GFit_H <- readRDS(file = "GFit_H.rds")
LFit_H <- readRDS(file = "LFit_H.rds")

## Input mrgsolve-based PBPK Model
source (file = "Pred_Rat.R")
source (file = "Pred_H.R")

##================================================================================================
## Fig. 4A: Time-varying profiles of PFOS concentration in maternal and fetal palsma druing pregnancy
##

## Defined the prediction function from pre-pregnant to gestational exposure
PBPK_H_G <- function(pars, DOSE, pred = FALSE) {
  
  Init <- Pred.preG (DOSE) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
  
  ## Get out of log domain
  Gpars_H <- lapply(pars, exp)          ## Return a list of exp (parametrs) 
  
  ## Exposure scenario for gestational exposure
  GBW          = 67                     ## Body weight before gestation (measrument data if available); 
  tinterval    = 24                     ## Time interval; 
  GTDOSE       = 7*41                   ## Total dosing/Dose times; Repeat oral dose from GA0 - GA40
  GDOSE        = DOSE                   ## Input oral dose  
  GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
  
  # To create a data set of 1 subject receiving DOSE every 24 hours 
  Gex.oral <- ev (ID   = 1,             ## One individual
                  time = 0,             ## Dossed strat time (GA0)
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Addtional doseing 
                  cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time from GA0 to GA 41 (weeks) + 1 day
  
  ## Simulation of exposure scenaior 
  Gout <- 
    Gmod_H %>% ## Gestational PBPK model
    init(APlas_free = Init$APlas_free, APTC = Init$APTC, 
         AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
         AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, 
         A_apical = Init$A_apical, Adif = Init$Adif,
         Aefflux = Init$Aefflux) %>% ## Input the intial concentrations 
    param (Gpars_H) %>% ## Update the parameter list with Gpars
    update(atol = 1E-6,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)   
  
  ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
  Goutdf = cbind.data.frame(Time     = Gout$time/(24*7),   ## Simulation time 
                            CPlas    = Gout$Plasma*1000,   ## Unit change from mg/L to ng/mL; PFOS conc. in maternal palsma (CPlas)
                            CL       = Gout$Liver*1000,    ## PFOS conc. in maternal liver (CL)
                            CPlas_pup  = Gout$CordB*1000,  ## PFOS conc. in cord blood (CPlas_pup)
                            CPla     = Gout$Placenta*1000) ## PFOS conc. in placenta (CPla) 
  
  if (pred) return (Goutdf)
  
  Init_G <- Gout %>% filter (time == 40*24*7) ## Extract the concentration at GA40
  return (Init_G)
}   


## Defined the prediction function for lactation period
PBPK_H_L <- function(Gpars_H, Lpars_H, DOSE) {
  
  Init_G <- PBPK_H_G (Gpars_H, DOSE)    ## Extract the initial concentration from gestational model
  
  ## Get out of log domain
  Lpars_H <- exp(Lpars_H)               ## Return a list of exp (parametrs) 
  
  ## Exposure scenario for Lactational exposure; through posnatal weeks (PNW)0 - PNW41 
  LBW          = 67                     ## human body weight before lactation
  tinterval    = 24                     ## Time interval; 
  LTDOSE       = 41*7                   ## Total dosing/Dose times; Repeat oral dose from PNW0 - PNW40
  LDOSE        = DOSE                   ## Repeat oral dose  (mg/kg/day);  
  LDOSEoral    = LDOSE*LBW              ## Amount of oral dose (mg)
  
  
  # To create a data set of 1 subject receiving DOSE every 24 hours
  Lex.oral <- ev (ID   = 1,             ## One individual
                  amt  = LDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = LTDOSE - 1,    ## Addtional doseing 
                  cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*7, 24) ## Simulation time from PNW0 to PNW41 plus 7 days without dosing
  
  
  Lout <- Lmod_H %>% # Lactational model
    init(APlas_free =  Init_G$APlas_free, APTC =  Init_G$APTC, AFil =  Init_G$AFil, AKb =  Init_G$AKb, ARest =  Init_G$ARest,
         AL =  Init_G$AL, AM =  Init_G$AM, AF =  Init_G$AF, A_baso =  Init_G$A_baso, A_apical = Init_G$A_apical, Adif =  Init_G$Adif,
         Aefflux =  Init_G$Aefflux, APlas_free_neo =  Init_G$APlas_Fet_free, 
         AL_neo =  Init_G$AL_Fet, ARest_neo =  Init_G$ARest_Fet) %>% # Input the intial concentration
    param (Lpars_H) %>% # Update the parameter list with pars
    update(atol = 1E-6, maxsteps = 500000) %>% # Atol: absolute tolerance parameter; hmax: The maximum step size; maxstep: maximum number of steps the solver will take when advancing from one time to the next         
    mrgsim_d (data = Lex.oral, tgrid = Ltsamp) 
  
  ## Extract the concentration; adjustted the unit from mg/l to ng/ml
  Loutdf = cbind.data.frame(Time       = (Lout$time/(24*7)) + 40, 
                            CPlas      = Lout$Plasma*1000,  ## PFOS concentration in maternal plasma 
                            CPlas_pup  = Lout$CPneo*1000,   ## PFOS concentration in neonatal plasma
                            CMilk      = Lout$Milk*1000)    ## PFOS concentration in milk
  
  Loutdf  <- rbind.data.frame (Loutdf)
  return (Loutdf)
  
}

### Loading best-fitted parameters 
bestpar_G <- GFit_H$par
bestpar_L <- LFit_H$par

Gout <- PBPK_H_G (bestpar_G, DOSE = 1.9E-7, pred = TRUE)
Lout <- PBPK_H_L (bestpar_L, bestpar_L, DOSE = 1.9E-7)

Goutdf <- Gout%>%select(Time = Time, CPlas = CPlas, CPlas_pup = CPlas_pup)
Loutdf <- Lout%>%select(Time = Time, CPlas = CPlas, CPlas_pup = CPlas_pup)
LMilk  <- Lout%>%select(Time = Time, CMilk = CMilk)

GL <- rbind (Goutdf, Loutdf) %>% distinct (Time, .keep_all = TRUE)

##================================================================================================
## Fig. 4a: Time-varying profiles of PFOS concentration in maternal and fetal/neolatal plasma druing lactation
##
p1 <- ggplot() + 
  geom_line (data = GL, aes(x = Time, y = CPlas_pup), colour = "black", lwd = 2, linetype = "twodash") + 
  geom_line (data = GL, aes(x = Time, y = CPlas), colour = "black", lwd = 2) +
  scale_x_continuous(breaks = seq(1, 71, 20), limits = c(1, 71), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0),limits = c(-1, 3))  

p1 <- p1 + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size = 2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text.x             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
    axis.text.y             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_blank(),   
    legend.position='none')   

p1

##================================================================================================
## Fig. 4b: Time-varying profiles of PFOS concentration in milk druing lactation
##


pmilk <- ggplot() + 
  geom_line (data = LMilk, aes(x = (Time-40), y = CMilk), colour = "steelblue4", lwd = 2, linetype = "longdash") + 
  scale_x_continuous(breaks = seq(1, 25, 4), limits = c(1, 25), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0),limits = c(0, 0.05))  

pmilk <- pmilk + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text.x             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
    axis.text.y             = element_text (size   = 12, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_blank(),   
    legend.position='none')   

pmilk

## Save the figures
# ggsave("Fig.4a.tiff",scale = 1,
#        plot = p1,
#        path = "C:/Users/weichunc/Desktop",
#        width = 14, height = 9, units = "cm", dpi=320)
# 
# ggsave("Fig.4b.tiff",scale = 1,
#        plot = pmilk,
#        path = "C:/Users/weichunc/Desktop",
#        width = 14, height = 9, units = "cm", dpi=320)


### Output the data ploting figure
# write.csv(GL, file = 'Data_figure_4a.csv')
# write.csv(LMilk, file = 'Data_figure_4b.csv')
# 


