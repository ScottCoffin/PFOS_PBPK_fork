#===========================================================================================================================
# Figure 5: Histogram of model simulations (right) compared with global measurements of PFOS in maternal plasma (round), cord blood (dimand) and breast milk (square) (left). 
# - Author : Wei-Chun Chou
# - Date   : March, 2020
# - Fig.5A : The PFOS time-crouse profiles in maternal palsma during gestation  
# - Fig.5B : The PFOS time-crouse profiles in maternal palsma during lactation
# - Fig.5C : The PFOS time-crouse profiles in pup palsma during lactation 
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
library(EnvStats)    ## R-package for function "rlnormTrunc"

## Input mrgsolve-based PBPK Model
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
source (file = "Pred_H.R")

## Input dataset
Human_MP_G <- HDat %>% filter (Cal != 1  & GA >= 39 & Matrix == "MP") 
Human_MP_G[is.na(Human_MP_G$SD), "SD"] <- 0
Human_MP_G$no <- seq(1, 9)

Human_CB_G <- HDat %>%  filter (Cal != 1  & GA <= 40 & GA >= 39 & Matrix == "CB")
Human_CB_G[is.na(Human_CB_G$SD), "SD"] <- 0
Human_CB_G$no <- seq(1, 9)

Human_Milk<-HDat %>% filter (Cal != 1  &  GA >= 39 & Matrix == "Milk") 
Human_Milk[is.na(Human_Milk$SD), "SD"] <- 0
Human_Milk$no <- c(10, 11, 12)

## Defined the prediction function from pre-pregnant to gestational exposure
PBPK_H_G <- function(pars, DOSE, pred = FALSE) {
    
    Init <- Pred.preG (DOSE) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ## Get out of log domain
    Gpars_H <- lapply(pars, exp)          ## Return a list of exp (parametrs for gestational model) from log scale
    
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
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>% ## Input the intial concentrations 
        param (Gpars_H) %>% ## Update the parameter list with Gpars
        update(atol = 1E-5,  maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
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

## Define the population prediction function from pre-pregnant to gestational exposure
PBPK_H_G_pop <- function(pars, pred = FALSE) {
    
    Init <- Pred.preG (DOSE = 7e-7) %>% filter (row_number()== n()) %>% select(-c("Plasma", "Liver", "Kidney"))
    
    ## Get out of log domain
    Gpars_H <- exp(pars)                  ## Return a list of exp (parametrs for gestational model) from log scale
    N <- 50000  
    ## Exposure scenario for gestational exposure
    tinterval    = 24                     ## Time interval; 
    GTDOSE       = 7*40                   ## Total dosing/Dose times; 
    
    ## Create "N" number of individal. 
    ## Each individal have a combinations of parmaeters
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
               DOSEoral  = BW*runif(N, min = 1.9e-7, max = 3.7e-6))
    
    # To crete the exposure scenario
    Gex.oral <- ev (ID   = 1:N,             ## Number of N individual
                    time = 0,               ## Dossed strat time 
                    amt  = idata$DOSEoral,  ## Amount of dose 
                    ii   = tinterval,       ## Time interval
                    addl = GTDOSE - 1,      ## Addtional doseing 
                    cmt  = "AST",           ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)      ## No replicate
    
    Gtsamp  = tgrid(0, tinterval*(GTDOSE - 1) + 24*7, 1) ## Simulation time from GD0 to GD 40 (weeks) + 1 day
    
    ## Simulation of exposure scenaior 
    Gout <- 
        Gmod_H %>% ## Gestational PBPK model
        init(APlas_free = Init$APlas_free, APTC = Init$APTC, AFil = Init$AFil, AKb = Init$AKb, ARest = Init$ARest,
             AL = Init$AL, AM = Init$AM, AF = Init$AF, A_baso = Init$A_baso, A_apical = Init$A_apical, Adif = Init$Adif,
             Aefflux = Init$Aefflux) %>% ## Input the intial concentrations 
        idata_set(idata) %>% ## Update the parameter list with Gpars
        update(atol = 1E-3,  maxsteps = 50000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps          
        mrgsim (data = Gex.oral, tgrid = Gtsamp)    
    
    if (pred) {return (as.data.frame(Gout))}
    
    outdf = Gout %>% filter (time == 24*39*7)
    
    ## Extract the concentration 
    outdf = cbind.data.frame( ID        = outdf$ID,
                              Time      = (outdf$time/(24*7)), 
                              CPlas     = outdf$Plasma*1000, 
                              CPlas_pup = outdf$CordB*1000)
    
    return (outdf)
    
}

  
## Exposure sceinario for lactational stage
PBPK_H_L_pop <- function(Lpars_H) {
    
    Init_G <- PBPK_H_G (GFit_H$par,  1.9e-7)
    
    ## Get out of log domain
    Lpars_H <- exp(Lpars_H)           ## Return a list of exp (parametrs for gestational model) from log scale
    N <- 50000
    ## Exposure scenario for gestational exposure
    #GBW          = 67                     ## Body weight before gestation (measrument data if available); Default value adopted from Garner et al. (2015)
    tinterval    = 24                     ## Time interval; 
    LTDOSE       = 7*30                   ## Total dosing/Dose times; 
    
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
            PRest     = rlnormTrunc (N, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
            PRest_neo = rlnormTrunc (N, min = exp(-2.02), max = exp(-1.24), meanlog = -1.63, sdlog = 0.13),
            QCC       = rnormTrunc (N, min = 6.76, max = 26.04, mean = 16.4, sd = 3.28),
            QKC       = rnormTrunc (N, min = 0.07, max = 0.28, mean = 0.175, sd = 0.035),
            RAFbaso   = rlnormTrunc (N, min = exp(-0.618), max = exp(-0.532), meanlog = -0.043, sdlog = 0.29),
            VKC       = rnormTrunc (N, min = 0.002, max = 0.01, mean = 0.004, sd = 0.001),
            RAFapi   = 0.525,
            DOSEoral  = 67*runif(N, min = 1.9e-7, max = 3.7e-6))
    
    
    # To create a data set of 1 subject receiving DOSE every 24 hours
    Lex.oral <- ev (ID   = 1:N,             ## One individual
                    amt  = idata$DOSEoral,     ## Amount of dose 
                    ii   = tinterval,     ## Time interval
                    addl = LTDOSE - 1,    ## Addtional doseing 
                    cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                    replicate = FALSE)    ## No replicate
    
    Ltsamp       = tgrid(0, tinterval*(LTDOSE - 1) + 24*7, 24) ## Simulation time from PND0 to PND74 with dosing of LTDOSE - 1+ 2 days (No dosing)
    
    
    Lout <- Lmod_H %>% # Lactational model
        init(APlas_free =  Init_G$APlas_free, APTC =  Init_G$APTC, AFil =  Init_G$AFil, AKb =  Init_G$AKb, ARest =  Init_G$ARest,
             AL =  Init_G$AL, AM =  Init_G$AM, AF =  Init_G$AF, A_baso =  Init_G$A_baso, A_apical = Init_G$A_apical, Adif =  Init_G$Adif,
             Aefflux =  Init_G$Aefflux, APlas_free_neo =  Init_G$APlas_Fet_free, 
             AL_neo =  Init_G$AL_Fet, ARest_neo =  Init_G$ARest_Fet) %>% # Input the intial concentration
        idata_set(idata) %>% # Update the parameter list with pars
        update(atol = 1E-3, maxsteps = 500000) %>% # Atol: absolute tolerance parameter; hmax: The maximum step size; maxstep: maximum number of steps the solver will take when advancing from one time to the next         
        mrgsim (data = Lex.oral, tgrid = Ltsamp) 
    
    
    ## Extract the concentration 
    Loutdf = cbind.data.frame(ID         = Lout$ID,
                              Time       = (Lout$time/(24*7)) + 40, 
                              CPlas      = Lout$Plasma*1000, 
                              CPlas_pup  = Lout$CPneo*1000,
                              CMilk      = Lout$Milk*1000)
    
    Loutdf  <- rbind.data.frame (Loutdf)
    return (Loutdf)
    
}





## Calculation of average PFOS concentrations in maternal palsma, milk and fetal plasma
H   <- PBPK_H_G_pop(GFit_H$par)%>%select(ID = ID, Time = Time, CPlas = CPlas, CPlas_pup = CPlas_pup)
H_2 <- PBPK_H_L_pop(LFit_H$par)%>%select(ID = ID, Time = Time, CMilk = CMilk)%>%filter (Time == 42)


H_CPlas <- H %>% summarize (median = quantile (CPlas, probs = 0.5), 
                            ci_05 = quantile (CPlas, probs = 0.025),
                            ci_95 = quantile (CPlas, probs = 0.975))

H_CPlas_pup <- H %>% summarize (median = quantile (CPlas_pup , probs = 0.5), 
                                ci_05 = quantile (CPlas_pup, probs = 0.025),
                                ci_95 = quantile (CPlas_pup, probs = 0.975))

H_CMilk <- H_2 %>% summarize (median = quantile (CMilk , probs = 0.5), 
                              ci_05 = quantile (CMilk, probs = 0.025 ),
                              ci_95 = quantile  (CMilk,probs = 0.975 ))

Human_MP_G %>% summarize (median = quantile (Conc , probs = 0.5), 
                          ci_95 = quantile (Conc , probs = 0.975),
                          ci_05 = quantile (Conc , probs = 0.025))

Human_CB_G %>% summarize (median = quantile (Conc , probs = 0.5), 
                          ci_95 = quantile (Conc , probs = 0.975),
                          ci_05 = quantile (Conc , probs = 0.025))

Human_Milk %>% summarize (median = quantile (Conc , probs = 0.5), 
                            ci_95 = quantile (Conc , probs = 0.975),
                            ci_05 = quantile (Conc , probs = 0.025))




### Plot figure 
######################
## Add the font to font database
windowsFonts("Times" = windowsFont("Times New Roman"))
p2 <- ggplot() + 
    geom_pointrange(data = Human_MP_G, 
                    mapping = aes(x = no, y = Conc, ymin = Conc - SD, ymax = Conc + SD),
                    color ="#00AFBB", fill = "#00AFBB", size = 0.8) +
    geom_pointrange(data = Human_CB_G, 
                    mapping = aes(x = no, y = Conc, ymin = Conc - SD, ymax = Conc + SD),
                    color ="gray", fill = "gray", shape = 23, stroke = 1, size = 0.8) +
    geom_pointrange(data = Human_Milk, 
                    mapping = aes(x = no, y = Conc*100, ymin = (Conc - SD)*100, ymax = (Conc + SD)*100),
                    color ="Orange", fill = "Orange", shape = 22, stroke = 1, size = 0.8) + 
    scale_x_continuous(breaks=seq(1:13)) + ylim (0, 40) + labs (x = "", y = "")



p2 <- p2 + 
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, size=2),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(), 
        axis.text.x             = element_blank(),    # tick labels along axes 
        axis.text.y             = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position='none')  



p3 <- ggplot(H, aes(ID, CPlas)) + geom_point() + ylim (0, 40) + theme_bw()
ggMarginal(p3, type = "histogram", margins = "y", size = 4,
           color="black", alpha = 0.4, bins = 30,
           fill = "#00AFBB", position="identity")

p4 <- ggplot(H, aes(ID, CPlas_pup)) + geom_point() + ylim (0, 40) + theme_bw()
ggMarginal(p4, type = "histogram", margins = "y", size = 4,
           color = "black", fill="gray", bins = 30, alpha = 0.4)


p5 <- ggplot(H_2, aes(ID, CMilk*100)) + geom_point() + ylim (0, 40) + theme_bw()
ggMarginal(p5, type = "histogram", margins = "y", size = 4,
           color = "black", fill="Orange", bins = 30, alpha = 0.4)

## Save your figures 

# ggsave("Fig.5a1.tiff",scale = 1,
#        plot = p2,
#        #path = "C:/Users/weichunc/Desktop",
#        width = 20, height = 10, units = "cm", dpi=320)
# 
# ggsave("Fig.5b1.tiff",scale = 1,
#        plot = ggMarginal(p3, type = "histogram", margins = "y", size = 4,
#                          color="black", alpha = 0.4, bins = 30,
#                          fill = "#00AFBB", position="identity"),
#        #path = "C:/Users/weichunc/Desktop",
#        width = 17, height = 12, units = "cm", dpi=320)
# 
# ggsave("Fig.5b2.tiff",scale = 1,
#        plot = ggMarginal(p4, type = "histogram", margins = "y", size = 4,
#                          color = "black", fill="gray", bins = 30, alpha = 0.4),
#        #path = "C:/Users/weichunc/Desktop",
#        width = 17, height = 12, units = "cm", dpi=320)
# 
# ggsave("Fig.5b3.tiff",scale = 1,
#        plot = ggMarginal(p5, type = "histogram", margins = "y", size = 4,
#                          color = "black", fill="Orange", bins = 40, alpha = 0.4),
#        #path = "C:/Users/weichunc/Desktop",
#        width = 17, height = 12, units = "cm", dpi=320)

