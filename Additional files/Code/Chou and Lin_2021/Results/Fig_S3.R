##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(gridExtra)     ## R-package for function "grid.arrange" 
library(EnvStats)

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")
source (file = "HMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_H <- mcode ("PreGHumanPBPK.code", PreGHumanPBPK.code) 
Gmod_H    <- mcode ("GHumanPBPK.code", GHumanPBPK.code)
Lmod_H    <- mcode ("LHumanPBPK.code", LHumanPBPK.code)

## Loading datasets
HDat <- read.csv(file = "Data_H.csv")

## Loading Fit results from rds files
GFit_H <- readRDS(file = "GFit_H.rds")
LFit_H <- readRDS(file = "LFit_H.rds")

## Input mrgsolve-based PBPK Model
source (file = "Pred_Rat.R")
source (file = "Pred_H.R")

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

## Read the data and later used in model calibration and evluation
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
Pop_pred.G <- function(DOSE, Init, N) {
    
    ## Exposure scenario for gestational exposure
    GBW          = 60                     ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; Repeat oral dose from GA0 - GA40 (weeks)
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    ## Repeat oral dose of GDOSE mg/kg/day;  
    ## Amount of oral dose
    idata <- 
        tibble(ID = 1:N) %>% 
        mutate(BW        = rnormTrunc (N, min = 27.89, max = 107.5, mean = 67.7, sd = 20.3),
               SA_Htc    = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
               SA_QC     = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
               SA_QK_P   = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
               Free      = rlnormTrunc (N, min = exp(-6.53), max = exp(-5.38), meanlog = -5.96, sdlog = 0.29),  
               Free_Fet  = rlnormTrunc (N, min = exp(-6.19), max = exp(-5.04), meanlog = -5.61, sdlog = 0.29),
               KeffluxC  = rlnormTrunc (N, min = exp(-4.82), max = exp(-3.67), meanlog = -4.24, sdlog = 0.29), 
               Ktrans1C  = rlnormTrunc (N, min = exp(-0.85), max = exp(0.29), meanlog = -0.28, sdlog = 0.29), 
               Ktrans2C  = rlnormTrunc (N, min = exp(-0.51), max = exp(0.65), meanlog = 0.07, sdlog = 0.29), 
               QCC       = rnormTrunc  (N, min = 6.76, max = 26.04, mean = 16.4, sd = 4.92), 
               Km_apical = 248,
               PPla      = 0.13,
               PL_Fet    = 0.58,
               PRest_Fet = 2.3,
               GDOSEoral = GDOSE*GBW
        )   
    
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    Gex.oral <- ev (ID   = 1:N,                 ## One individual
                    amt  = idata$GDOSEoral,     ## Amount of dose 
                    ii   = tinterval,           ## Time interval
                    addl = GTDOSE - 1,          ## Addtional doseing 
                    cmt  = "AST",               ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)          ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 24) ## Simulation time from GD0 to GD 40
    
    ## Simulation of exposure scenaior (Repeated oral dose to 1/2/3/5/10 mg/kg)
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>%  ## Input the intial concentration
        idata_set(idata) %>% ## Update the parameter list with Gpars
        update(atol = 1E-3,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim (data = Gex.oral, tgrid = Gtsamp)   
    
    ## Extract the "Time", "CPlas", "CL" , "CPlas_P", "CPla" and "Curine" from Gout
    Goutdf = cbind.data.frame(ID       = Gout$ID,
                              Time     = Gout$time, 
                              CPlas    = Gout$Plasma*1000, 
                              CL       = Gout$Liver*1000,
                              CFPlas   = Gout$CordB*1000,
                              CPla     = Gout$Placenta*1000,
                              CFL      = Gout$Fliver*1000)
    
    return (Goutdf) # Return Goutdf
}


## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)


Plot_G <- rbind.data.frame(
Pop_pred.G (DOSE = DOSE_A1, Init = Init_A1, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_1", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A1, Init = Init_A1, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>%mutate(DOSE_Group = "A_1", Tissue = "CB"),
Pop_pred.G (DOSE = DOSE_A2, Init = Init_A2, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_2", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A2, Init = Init_A2, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_2", Tissue = "CB"),
Pop_pred.G (DOSE = DOSE_A3, Init = Init_A3, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_3", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A3, Init = Init_A3, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_3", Tissue = "CB"),
Pop_pred.G (DOSE = DOSE_A4, Init = Init_A4, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_4", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A4, Init = Init_A4, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_4", Tissue = "CB"),
Pop_pred.G (DOSE = DOSE_A5, Init = Init_A5, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_5", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A5, Init = Init_A5, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPla)   %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_5", Tissue = "Pla"),
Pop_pred.G (DOSE = DOSE_A5, Init = Init_A5, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFL)    %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_5", Tissue = "Fli"),
Pop_pred.G (DOSE = DOSE_A6, Init = Init_A6, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_6", Tissue = "MP"),
Pop_pred.G (DOSE = DOSE_A6, Init = Init_A6, N = 1000)%>% filter (Time == 39*24*7) %>% select (ID = ID, Time = Time, Conc = CFPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "A_6", Tissue = "CB")
)


OBS_G <- rbind.data.frame(
OBS_A1_MP  <- OBS.A1_MP  %>% mutate(DOSE_Group = "A_1", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS_A1_CB  <- OBS.A1_CB  %>% mutate(DOSE_Group = "A_1", Tissue = "CB") %>% filter(row_number()==n()) %>% rename(Conc = CFPlas),
OBS.A2_MP  <- OBS.A2_MP  %>% mutate(DOSE_Group = "A_2", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS.A2_CB  <- OBS.A2_CB  %>% mutate(DOSE_Group = "A_2", Tissue = "CB") %>% filter(row_number()==n()) %>% rename(Conc = CFPlas),
OBS.A3_MP  <- OBS.A3_MP  %>% mutate(DOSE_Group = "A_3", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS.A3_CB  <- OBS.A3_CB  %>% mutate(DOSE_Group = "A_3", Tissue = "CB") %>% filter(row_number()==n()) %>% rename(Conc = CFPlas),
OBS.A4_MP  <- OBS.A4_MP  %>% mutate(DOSE_Group = "A_4", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS.A4_CB  <- OBS.A4_CB  %>% mutate(DOSE_Group = "A_4", Tissue = "CB") %>% filter(row_number()==n()) %>% rename(Conc = CFPlas),
OBS.A5_MP  <- OBS.A5_MP  %>% mutate(DOSE_Group = "A_5", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS.A5_Pla <- OBS.A5_Pla %>% mutate(DOSE_Group = "A_5", Tissue = "Pla")%>% filter(row_number()==n()) %>% rename(Conc = CPla),
OBS.A5_Fli <- OBS.A5_Fli %>% mutate(DOSE_Group = "A_5", Tissue = "Fli")%>% filter(row_number()==n()) %>% rename(Conc = CFL),
OBS.A6_MP  <- OBS.A6_MP  %>% mutate(DOSE_Group = "A_6", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
OBS.A6_CB  <- OBS.A6_CB  %>% mutate(DOSE_Group = "A_6", Tissue = "CB") %>% filter(row_number()==n()) %>% rename(Conc = CFPlas)
)


### Forest-like plot

labe_A <- c("Inoue et al. (2004)", "Fei et al. (2007)", "Kato et al. (2014)", "Pan et al. (2017)", "Mamsen  et al. (2019)", "Midasch et al. (2007)")
names(labe_A)  <- c("A_1", "A_2", "A_3","A_4", "A_5", "A_6")

p1 = ggplot( )+
    geom_pointrange(data = Plot_G,
                    aes(x = Tissue, y = Median, ymin = lb, ymax = ub, shape = Tissue), size = 0.5)+
    geom_point(data = OBS_G, aes(x = Tissue, y = Conc, colour = "darkred"), size = 2.5) +
    xlab('')+ ylab("")+ #+ scale_color_brewer(palette="Dark2") +
    facet_wrap(~DOSE_Group,strip.position="left",nrow=9,scales = "free_y", labeller = labeller(DOSE_Group = labe_A)) +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(size = rel(1.5), face = "bold",colour = "black"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(size = rel(1.5), hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()

p1 = p1 +  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                   labels = scales::trans_format("log10", scales::math_format(10^.x))) 



#==================================================================================================
# Model calibration for lactational PBPK model based on human biomonitoring data
# Abrreviation: MP: Maternal plasma, NP: Neonatal plasma, CB: Cord blood
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
    
    Init <- Pred.preG (DOSE) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ## Get out of log domain
    Gpars <- lapply(GFit_H$par, exp)        ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for gestational exposure
    GBW          = 67                     ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; 
    GDOSE        = DOSE                   ## Input oral dose  
    GDOSEoral    = GDOSE*GBW              ## Amount of oral dose
    
    # To create a data set of 1 subject receiving DOSE every 24 hours 
    Gex.oral <- ev (ID   = 1,             ## One individual
                    time = 0,             ## Dossed strat time (GA0)
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
Pop_pred.L <- function(DOSE, Init_G, N) {
    
    ## Exposure scenario for Lactational exposure;
    LBW          = 67                     ## Rat body weight 
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 30*7                   ## Total dosing/Dose times;
    LDOSE        = DOSE                   ## Repeat oral dose  (mg/kg/day);  
    LDOSEoral    = LDOSE*LBW              ## Amount of oral dose (mg)
    
    
    ## Repeat oral dose of GDOSE mg/kg/day;  
    ## Amount of oral dose
    idata <- 
        tibble(ID = 1:N) %>% 
        mutate(
            SA_BW     = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
            SA_BW_neo = rnormTrunc (N, min = 0.61, max = 1.39, mean = 1, sd = 0.2),
            Free      = rlnormTrunc (N, min = exp(-3.68), max = exp(-2.94), meanlog = -3.32, sdlog = 0.07), 
            Free_neo  = rlnormTrunc (N, min = exp(-4.71), max = exp(-3.87), meanlog = -4.29, sdlog = 0.13),
            KeffluxC  = rlnormTrunc (N, min = exp(-2.52), max = exp(-1.36), meanlog = -1.94, sdlog = 0.19),
            KeffluxC_neo = rlnormTrunc (N, min = exp(-2.52), max = exp(-1.36), meanlog = -1.94, sdlog = 0.19),
            PAMilkC    = rlnormTrunc (N, min = exp(-6.49), max = exp(-5.35), meanlog = -3.6, sdlog = 0.29),
            PRest     = rlnormTrunc (N, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
            PRest_neo = rlnormTrunc (N, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
            QCC       = rnormTrunc (N, min = 6.76, max = 26.04, mean = 16.4, sd = 3.28),
            QKC       = rnormTrunc (N, min = 0.07, max = 0.28, mean = 0.175, sd = 0.035),
            RAFbaso   = rlnormTrunc (N, min = exp(-0.618), max = exp(-0.532), meanlog = -0.043, sdlog = 0.29),
            VKC       = rnormTrunc (N, min = 0.002, max = 0.01, mean = 0.004, sd = 0.001),
            RAFapi   = 0.525,
            LDOSEoral = LDOSE*67
        )   
    
    
    # To create a data set of 1 subject receiving GDOSE every 24 hours for 1 total doses
    Lex.oral <- ev (ID   = 1:N,                 ## One individual
                    amt  = idata$LDOSEoral,     ## Amount of dose 
                    ii   = tinterval,           ## Time interval
                    addl = LTDOSE - 1,          ## Addtional doseing 
                    cmt  = "AST",               ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)          ## No replicate
    
    Ltsamp  = tgrid(0, tinterval*(LTDOSE - 1) + 24*1, 24) ## Simulation time 
    
    
    Lout <- Lmod_H %>% # Lactational model
             init(APlas_free =  Init_G$APlas_free, APTC =  Init_G$APTC, AFil =  Init_G$AFil, AKb =  Init_G$AKb, ARest =  Init_G$ARest,
             AL =  Init_G$AL, AM =  Init_G$AM, AF =  Init_G$AF, A_baso =  Init_G$A_baso, A_apical = Init_G$A_apical, Adif =  Init_G$Adif,
             Aefflux =  Init_G$Aefflux, APlas_free_neo =  Init_G$APlas_Fet_free, 
             AL_neo =  Init_G$AL_Fet, ARest_neo =  Init_G$ARest_Fet) %>% # Input the intial concentration
             idata_set(idata) %>% ## Update the parameter list with Gpars
             update(atol = 1E-3,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
             mrgsim (data = Lex.oral, tgrid = Ltsamp)   
    
    ## Extract the concentration 
    Loutdf = cbind.data.frame(ID        = Lout$ID,
                              Time      = Lout$time + 40*7*24, 
                              CPlas     = Lout$Plasma*1000, 
                              CPlas_pup = Lout$CPneo*1000,
                              CMilk     = Lout$Milk*1000)
    
    outdf  <- rbind.data.frame (Loutdf)
    return (Loutdf)
    
}


## Create a cost fuction and later used in model optimization  
## Estimate the model residual with experimental data by modCost function (from FME package)

Plot_L <- rbind.data.frame(
    Pop_pred.L (DOSE = DOSE_B1, Init_G = Init_G1, N = 1000)%>% filter (Time == 43*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "B_1", Tissue = "MP"),
    Pop_pred.L (DOSE = DOSE_B1, Init_G = Init_G1, N = 1000)%>% filter (Time == 43*24*7) %>% select (ID = ID, Time = Time, Conc = CMilk) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>%mutate(DOSE_Group = "B_1", Tissue = "Milk"),
    Pop_pred.L (DOSE = DOSE_B2, Init_G = Init_G2, N = 1000)%>% filter (Time == 53*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "B_2", Tissue = "MP"),
    Pop_pred.L (DOSE = DOSE_B3, Init_G = Init_G3, N = 1000)%>% filter (Time == 66*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "B_3", Tissue = "MP"),
    Pop_pred.L (DOSE = DOSE_B3, Init_G = Init_G3, N = 1000)%>% filter (Time == 66*24*7) %>% select (ID = ID, Time = Time, Conc = CPlas_pup)  %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "B_3", Tissue = "NP"),
    Pop_pred.L (DOSE = DOSE_B4, Init_G = Init_G4, N = 1000)%>% filter (Time == 53*24*7) %>% select (ID = ID, Time = Time, Conc = CMilk) %>% summarize (Median = quantile(Conc, prob = 0.5), lb = quantile(Conc, prob = 0.025), ub = quantile(Conc, prob = 0.975)) %>% mutate(DOSE_Group = "B_4", Tissue = "Milk")
)


OBS_L <- rbind.data.frame(
    OBS_B1  <- OBS.B1_MP  %>% mutate(DOSE_Group = "B_1", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas) ,
    OBS.B2  <- OBS.B1_Milk  %>% mutate(DOSE_Group = "B_1", Tissue = "Milk") %>% filter(row_number()==n()) %>% rename(Conc = CMilk),
    OBS.B3  <- OBS.B2_MP  %>% mutate(DOSE_Group = "B_2", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
    OBS.B4  <- OBS.B3_MP  %>% mutate(DOSE_Group = "B_3", Tissue = "MP") %>% filter(row_number()==n()) %>% rename(Conc = CPlas),
    OBS.B5  <- OBS.B3_NP  %>% mutate(DOSE_Group = "B_3", Tissue = "NP") %>% filter(Time == 66*24*7) %>% rename(Conc = CPlas_pup),
    OBS.B6  <- OBS.B4_Milk  %>% mutate(DOSE_Group = "B_4", Tissue = "Milk") %>% filter(row_number()==n()) %>% rename(Conc = CMilk)
)


### Forest-like plot

labe_B <- c("Kärrman et al. (2007)", "von Ehrenstein et al. (2009)", "Fromme et al. (2010)", "Lee et al. (2018)")
names(labe_B)  <- c("B_1", "B_2", "B_3","B_4")

p2 = ggplot( )+
    geom_pointrange(data = Plot_L,
                    aes(x = Tissue, y = Median, ymin = lb, ymax = ub, shape = Tissue), size = 0.5)+
    geom_point(data = OBS_L, aes(x = Tissue, y = Conc, colour = "darkred"), size = 2.5) +
    xlab('')+ ylab("")+ scale_shape_manual(values=c(13, 10, 9)) +
    facet_wrap(~DOSE_Group, strip.position="left",nrow=9,scales = "free_y", labeller = labeller(DOSE_Group = labe_B)) +
    theme(plot.title=element_text(size=16,face="bold"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(size = rel(1.5), face = "bold",colour = "black"),
          axis.title=element_text(size=12,face="bold"),
          strip.text.y = element_text(size = rel(1.5), hjust=0,vjust = 1,angle=180,face="bold"))+
    coord_flip()

p2 = p2 +  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) 



### total ##
labe <- c(labe_A, labe_B)
Plotdata <- rbind.data.frame (Plot_G, Plot_L)
OBSdata  <- rbind.data.frame (OBS_G, OBS_L)

p3 = ggplot( )+
    geom_pointrange(data = Plotdata,
                    aes(x = Tissue, y = Median, ymin = lb, ymax = ub, shape = Tissue), size = 0.5)+
    geom_point(data = OBSdata, aes(x = Tissue, y = Conc, colour = "red"), size = 2.5) +
    xlab('')+ ylab("") + 
    facet_wrap(~DOSE_Group,  nrow=10, strip.position="left", scales = "free_y", labeller = labeller(DOSE_Group = labe)) +
    theme_bw() +
    theme(plot.title   = element_text(size=16,face="bold"),
          text         = element_text (family = "Times"),   # text front (Time new roman)
          axis.text.y  = element_blank(),
          #axis.ticks.y = element_blank(),
          axis.text.x  = element_text(size = rel(2), face = "bold",colour = "black"),
          axis.title   = element_text(size = 12,face = "bold"),
          strip.text.y = element_text(size = rel(1.5), hjust = 0,vjust = 1,angle = 180, face="bold"))+
    coord_flip()

p3 = p3 +  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) 

p3

## Save the figures
ggsave ("Fig.S3.tiff",scale = 1,
        plot = p3,
        width = 35, height = 22, units = "cm", dpi = 320)







