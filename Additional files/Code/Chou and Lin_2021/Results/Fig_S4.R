#  Loading requried R package
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

## Define the prediction function 
## Exposure scenario:  Dosing druing pregnnacy from pre-pregnant to gestational stage
pred.preG <- function (N) {
  
  ## Defined the exposure scenario for age specific data
  idata <- tibble(ID = 1:N) %>%
    mutate(BW = rnorm (N, mean = 57, sd = 10.3),
                DOSEoral = BW*runif(N, 0.19e-6, 4.4e-6))
  
  ## Exposure scenario for gestational exposure
  tinterval    = 24                     ## Time interval; 
  TDOSE        = 365*30                   ## Total dosing/Dose times; 
  
  ## To crete the exposure scenario
  ex.oral <- ev (ID   = 1:N,             ## Number of N individual
                  amt  = idata$DOSEoral,  ## Amount of dose 
                  ii   = tinterval,       ## Time interval
                  addl = TDOSE - 1,      ## Addtional doseing 
                  cmt  = "AST",           ## The dosing comaprtment: AST Stomach  
                  replicate = FALSE)      ## No replicate
           
           
           
  tsamp <- tgrid(0, tinterval*(TDOSE - 1) + 24*7, 24) # simulation time from 0 to 30 years old, and the time interval is 0.1
  
  ## Exposure scenario: A dynamic exposure model for PFOS daily intakes 
  
  
  out <- PreGmod_H %>%
        idata_set(idata) %>%
        update(atol = 1E-3, maxsteps = 50000) %>%          
        mrgsim(data = ex.oral, tgrid = tsamp)#%>%
        #filter (time >= 365*24) 
    
    outdf = cbind.data.frame( ID        = out$ID,
                              Time      = out$time, 
                              CPlas     = out$Plasma*1000) 

    
  return (outdf)
}

H <- pred.preG (1000) 

H_Summary <- H %>% filter (Time == 365*30*24)%>%
             summarize (
                 median = quantile (CPlas, probs = 0.5),
                 ci_05 = quantile (CPlas, probs = 0.025),
                 ci_95 = quantile (CPlas, probs = 0.975))


# Two NHANES datasets for female: citation: CDC et al. (2018);https://www.cdc.gov/exposurereport/pdf/FourthReport_UpdatedTables_Volume1_Mar2018.pdf
# 1. 1999-2000: GM: 28 (24.6 - 31.8);
# 2. 2003-2004: GM: 18.4 (17 - 20);
# 3. 2005-2006: GM: 14.4 (13.3-15.4);
# 4. 2007-2008: GM: 10.7 (9.7-11.7);
# 5. 2009-2010: GM: 7.65 (6.73-8.71);
# 6. 2011-2012: GM: 5.1 (4.7-5.53);
# 6. 2013-2014: GM: 3.96 (3.60-4.35);

df <- data.frame(
  trt = factor(c("This study", "NHANES 2005-2006", "NHANES 2007-2008",
                 "NHANES 2009-2010","NHANES 2011-2012", "NHANES 2013-2014")),
  Median = c(7.43, 14.4, 10.7, 7.65, 5.1, 3.96),
  group = factor(c(0.96, 2, 3, 4, 5, 6)),
  upper = c(13.54,15.4, 11.7, 8.71, 5.53, 4.35),
  lower = c(3.01,13.3, 9.7, 6.73, 4.7, 3.60)
)

df$trt <- factor(df$trt, levels = c("This study", "NHANES 2005-2006", "NHANES 2007-2008",
                                    "NHANES 2009-2010","NHANES 2011-2012", "NHANES 2013-2014"))

p <- ggplot(df, aes(trt, Median, shape = group))
p1<-p + geom_pointrange(aes(ymin = lower, ymax = upper),size = 1.2) + xlab("") + ylab("") 

p1 <-p1 + theme_bw() +
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, size=2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 25, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 25, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none')  


# Save the figure
# ggsave("Figure.S4.tiff",scale = 1,
#        plot = p1 ,
#        width = 35, height = 25, units = "cm", dpi = 320)
# 
# 





