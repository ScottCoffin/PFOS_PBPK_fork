##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r
library(magrittr)    ## R-package for the pipe, %>% , 
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for ggplot2
library(FME)         ## R-package for model fitting
library(minpack.lm)  ## R-package for model fitting
library(ggExtra)     ## R-package for function "grid.arrange" 
library(EnvStats)    ## R-package for defined the truncated distribution

## Input mrgsolve-based PBPK Model
source (file = "RMod.R")

## Build mrgsolve-based PBPK Model
PreGmod_R <- mcode ("PreG_RatPBPK.code", PreG_RatPBPK.code)
Gmod_R    <- mcode ("GRatPBPK.code", GRatPBPK.code)
Lmod_R    <- mcode ("LRatPBPK.code", LRatPBPK.code)

## Loading datasets
RDat <- read.csv(file = "Data_R.csv") 

## Loading Fit results from rds files
GFit_R <- readRDS(file = "GFit_R.rds")
LFit_R <- readRDS(file = "LFit_R.rds")

## Input mrgsolve-based PBPK Model
source (file = "Pred_Rat.R")


## Defined the prediction function
pred.C <- function(Gpars, Lpars, DOSE) {
    
    PreG_BW          = 0.20 ## Female rat body weight (premating); 
    tinterval        = 24   ## Time interval; 
    PreG_TDOSE       = 42 + 14  ## Total dosing/Dose times; Repeat oral dose from beginning 42 days prior to cohabitation; (14 days matting)
    PreG_DOSE        = DOSE ## Repeat oral dose (mg/kg/day);  
    PreG_DOSEoral    = PreG_DOSE*PreG_BW ## Amount of oral dose
    
    
    # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
    PreG_ex.oral <- ev (ID   = 1,             ## One individual
                        amt  = PreG_DOSEoral, ## Amount of dose 
                        ii   = tinterval,     ## Time interval
                        addl = PreG_TDOSE - 1,## Addtional doseing 
                        cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                        replicate = FALSE)    ## No replicate
    
    ## Simulation time from 6 weeks prior to mating 
    PreG_tsamp = tgrid(0, tinterval*(PreG_TDOSE - 1) + 24*1, 1) 
    
    PreG_out <- PreGmod_R %>% 
        update(atol = 1E-6, maxsteps = 500000) %>%          
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
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*1, 1) ## Exposure time from GD0 to GD 20 and plus 1 day without dosing
    
    ## Simulation 
    Gout_1 <- 
        Gmod_R %>%
        init(APlas_free = PreG_Init$APlas_free, APTC = PreG_Init$APTC, AFil = PreG_Init$AFil, AKb = PreG_Init$AKb, ARest = PreG_Init$ARest,
             AL = PreG_Init$AL, AM = PreG_Init$AM, AF = PreG_Init$AF, A_baso = PreG_Init$A_baso, A_apical = PreG_Init$A_apical, Adif = PreG_Init$Adif,
             Aefflux = PreG_Init$Aefflux) %>%
        param (Gpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Gout_2 <- 
        Gmod_R %>%
        param (Gpars) %>%
        update(atol = 1E-6, maxsteps = 500000) %>%          
        mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
    
    Goutdf_1 = cbind.data.frame(Time      = Gout_1$time, 
                              CPlas     = Gout_1$Plasma,
                              CL        = Gout_1$Liver,
                              CPlas_pup = Gout_1$Plasma_Fet,
                              CL_pup    = Gout_1$Liver_Fet) 
    
    Goutdf_2 = cbind.data.frame(Time      = Gout_2$time, 
                                CPlas     = Gout_2$Plasma,
                                CL        = Gout_2$Liver,
                                CPlas_pup = Gout_2$Plasma_Fet,
                                CL_pup    = Gout_2$Liver_Fet)
    
    GInit <- Gout_1 %>% filter (time == 21*24) 
    
    ### 
    Lpars <- exp(Lpars)          ## Return a list of exp (parametrs for gestational model) from log scale
    
    ## Exposure scenario for Lactational exposure; through PND0 - PND21 (end of lactation)
    LBW          = 0.26                  ## Rat body weight during lactation from Luebker et al., 2005b
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
                    tinf = 0.01,         ## Infusion time;  
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
    
    return (list(Goutdf_1, Loutdf, Goutdf_2))
}



# Monte Carlo Simulation
# create several empty list ()
## Following dosing 0.1 mg/kg/day during gestaion and lactation
outdf_0.1_G  <-list()
outdf_0.1_L  <-list()
outdf_0.1_G2 <-list()

## Following dosing 0.3 mg/kg/day during gestaion and lactation
outdf_0.3_G  <-list()
outdf_0.3_L  <-list()
outdf_0.3_G2 <-list()

## Following dosing 0.4 mg/kg/day during gestaion and lactation
outdf_0.4_G  <-list()
outdf_0.4_L  <-list()
outdf_0.4_G2 <-list()

## Following dosing 1 mg/kg/day during gestaion and lactation
outdf_1.0_G  <-list()
outdf_1.0_L  <-list()
outdf_1.0_G2 <-list()

## Following dosing 1.6 mg/kg/day during gestaion and lactation
outdf_1.6_G  <-list()
outdf_1.6_L  <-list()
outdf_1.6_G2 <-list()

## Following dosing 3.2 mg/kg/day during gestaion and lactation
outdf_3.2_G  <-list()
outdf_3.2_L  <-list()
outdf_3.2_G2 <-list()


## Monte carlo simulation 
N = 50 # For testing purpose, we use 50 iterations. 
# Use N = 1000 to reproduce results in Table 4. The results using 50 vs. 1000 are similar.

for (i in 1:N) {
    pars_G  <- log(c(
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
    
    pars_L  <- log(c(
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
    
    outdf_0.1_G [[i]] <- pred.C(pars_G, pars_L, 0.1)[[1]] %>% mutate(ID = i, Dose_group = 0.1)
    outdf_0.1_L [[i]] <- pred.C(pars_G, pars_L, 0.1)[[2]] %>% mutate(ID = i, Dose_group = 0.1)
    outdf_0.1_G2 [[i]] <- pred.C(pars_G, pars_L, 0.1)[[3]] %>% mutate(ID = i, Dose_group = 0.1)
    
    outdf_0.3_G [[i]] <- pred.C(pars_G, pars_L, 0.3)[[1]] %>% mutate(ID = i, Dose_group = 0.3)
    outdf_0.3_L [[i]] <- pred.C(pars_G, pars_L, 0.3)[[2]] %>% mutate(ID = i, Dose_group = 0.3)
    outdf_0.3_G2 [[i]] <- pred.C(pars_G, pars_L, 0.3)[[3]] %>% mutate(ID = i, Dose_group = 0.3)
    
    outdf_0.4_G [[i]] <- pred.C(pars_G, pars_L, 0.4)[[1]] %>% mutate(ID = i, Dose_group = 0.4)
    outdf_0.4_L [[i]] <- pred.C(pars_G, pars_L, 0.4)[[2]] %>% mutate(ID = i, Dose_group = 0.4)
    outdf_0.4_G2 [[i]] <- pred.C(pars_G, pars_L, 0.4)[[3]] %>% mutate(ID = i, Dose_group = 0.4)
    
    outdf_1.0_G [[i]] <- pred.C(pars_G, pars_L, 1)[[1]] %>% mutate(ID = i, Dose_group = 1)
    outdf_1.0_L [[i]] <- pred.C(pars_G, pars_L, 1)[[2]] %>% mutate(ID = i, Dose_group = 1)
    outdf_1.0_G2 [[i]] <- pred.C(pars_G, pars_L, 1)[[3]] %>% mutate(ID = i, Dose_group = 1)
    
    outdf_1.6_G [[i]] <- pred.C(pars_G, pars_L, 1.6)[[1]] %>% mutate(ID = i, Dose_group = 1.6)
    outdf_1.6_L [[i]] <- pred.C(pars_G, pars_L, 1.6)[[2]] %>% mutate(ID = i, Dose_group = 1.6)
    outdf_1.6_G2 [[i]] <- pred.C(pars_G, pars_L, 1.6)[[3]] %>% mutate(ID = i, Dose_group = 1.6)
    
    outdf_3.2_G [[i]] <- pred.C(pars_G, pars_L, 3.2)[[1]] %>% mutate(ID = i, Dose_group = 3.2)
    outdf_3.2_L [[i]] <- pred.C(pars_G, pars_L, 3.2)[[2]] %>% mutate(ID = i, Dose_group = 3.2)
    outdf_3.2_G2 [[i]] <- pred.C(pars_G, pars_L, 3.2)[[3]] %>% mutate(ID = i, Dose_group = 3.2)
    
    cat(paste ("iteration = ", i, "\n"))
}

## The rnage of 0.1 mg/kg/day dose group
Range_0.1_GD20 <- do.call (rbind, outdf_0.1_G2)%>% 
    filter (Time == 20*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.1_GD21 <- do.call (rbind, outdf_0.1_G)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.1_PND5 <- do.call (rbind, outdf_0.1_L)%>% 
    filter (Time == 24*5) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.1_PND21 <- do.call (rbind, outdf_0.1_L)%>% 
    filter (Time == 24*21) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

## The rnage of 0.3 mg/kg/day dose group
Range_0.3_GD20 <- do.call (rbind, outdf_0.3_G2)%>% 
    filter (Time == 20*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.3_GD21 <- do.call (rbind, outdf_0.3_G)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.3_PND5 <- do.call (rbind, outdf_0.3_L)%>% 
    filter (Time == 24*5) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.3_PND21 <- do.call (rbind, outdf_0.3_L)%>% 
    filter (Time == 24*21) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))


## The rnage of 0.4 mg/kg/day dose group
Range_0.4_GD20 <- do.call (rbind, outdf_0.4_G2)%>% 
    filter (Time == 20*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.4_GD21 <- do.call (rbind, outdf_0.4_G)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.4_PND5 <- do.call (rbind, outdf_0.4_L)%>% 
    filter (Time == 24*5) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_0.4_PND21 <- do.call (rbind, outdf_0.4_L)%>% 
    filter (Time == 24*21) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

## The rnage of 1 mg/kg/day dose group
Range_1.0_GD20 <- do.call (rbind, outdf_1.0_G2)%>% 
    filter (Time == 24*20) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.0_GD21 <- do.call (rbind, outdf_1.0_G)%>% 
    filter (Time == 24*21) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.0_PND5 <- do.call (rbind, outdf_1.0_L)%>% 
    filter (Time == 5*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.0_PND21 <- do.call (rbind, outdf_1.0_L)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

## The rnage of 1.6 mg/kg/day dose group
Range_1.6_GD20 <- do.call (rbind, outdf_1.6_G2)%>% 
    filter (Time == 20*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.6_GD21 <- do.call (rbind, outdf_1.6_G)%>% 
    filter (Time == 24*21) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.6_PND5 <- do.call (rbind, outdf_1.6_L)%>% 
    filter (Time == 5*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_1.6_PND21 <- do.call (rbind, outdf_1.6_L)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))


## The rnage of 3.2 mg/kg/day dose group
Range_3.2_GD20 <- do.call (rbind, outdf_3.2_G2)%>% 
    filter (Time == 20*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))


Range_3.2_GD21 <- do.call (rbind, outdf_3.2_G)%>% 
    filter (Time == 21*24) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))


Range_3.2_PND5 <- do.call (rbind, outdf_3.2_L)%>% 
    filter (Time == 5*24 & CPlas_pup < 1000 & CL_pup < 1000) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

Range_3.2_PND21 <- do.call (rbind, outdf_3.2_L)%>% 
    filter (Time == 21*24& CPlas_pup < 1000 & CL_pup < 1000) %>% 
    summarize (
        Mean_CPlas      = mean(CPlas), 
        SD_CPlas        = sd (CPlas), 
        Mean_CPlas_pup  = mean(CPlas_pup), 
        SD_CPlas_pup    = sd (CPlas_pup),
        Mean_CL         = mean(CL), 
        SD_CL           = sd(CL),
        Mean_CL_pup     = mean(CL_pup), 
        SD_CL_pup       = sd(CL_pup))

## Make the table
# The mean and SD on the 0.1 dose group on GD20, GD21, PND5 and PND21
Range_0.1 = rbind.data.frame(
    GD20  = Range_0.1_GD20,
    GD21  = Range_0.1_GD21,
    PND5  = Range_0.1_PND5,
    PND21 = Range_0.1_PND21)

Range_0.1  
# The mean and SD on the 0.3 dose group on GD20, GD21, PND5 and PND21
Range_0.3 = rbind.data.frame(
    GD20  = Range_0.3_GD20,
    GD21  = Range_0.3_GD21,
    PND5  = Range_0.3_PND5,
    PND21 = Range_0.3_PND21)

Range_0.3  
# The mean and SD on the 0.4 dose group on GD20, GD21, PND5 and PND21
Range_0.4 = rbind.data.frame(
    GD20  = Range_0.4_GD20,
    GD21  = Range_0.4_GD21,
    PND5  = Range_0.4_PND5,
    PND21 = Range_0.4_PND21)

Range_0.4
# The mean and SD on the 1 dose group on GD20, GD21, PND5 and PND21
Range_1.0 = rbind.data.frame(
    GD20  = Range_1.0_GD20,
    GD21  = Range_1.0_GD21,
    PND5  = Range_1.0_PND5,
    PND21 = Range_1.0_PND21)

Range_1.0

# The mean and SD on the 1.6 dose group on GD20, GD21, PND5 and PND21
Range_1.6 = rbind.data.frame(
    GD20  = Range_1.6_GD20,
    GD21  = Range_1.6_GD21,
    PND5  = Range_1.6_PND5,
    PND21 = Range_1.6_PND21)

Range_1.6

# The mean and SD on the 3.2 dose group on GD20, GD21, PND5 and PND21
Range_3.2 = rbind.data.frame(
    GD20  = Range_3.2_GD20,
    GD21  = Range_3.2_GD21,
    PND5  = Range_3.2_PND5,
    PND21 = Range_3.2_PND21)


Range_3.2








