##These are model codes for building the route-specific PBTK models, get the prediction, and sensitivity analysis.  The codes were mdified based on Chou et al., 2019 (Environ Int, 2019, 129: 408-422).
## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)    # ggplot is the basic package for creating plots. ggplots is version 2 of ggplot. ggjoy uses some functions in ggplot2.
library(reshape)     # Package for melt function to reshape the table

##Codes for PFOS and PFHxS can be obtained by replacing the parameters, which are summarized in the Table S7

## Build mrgsolve-based PBTK Model

mod.oral <- mcode ("PFOApbtk.oral", PFOAPBTK.oral.code)
mod.nasal <- mcode ("PFOApbtk.nasal", PFOAPBTK.nasal.code)
mod.dermal <- mcode ("PFOApbtk.dermal", PFOAPBTK.dermal.code)
mod.iv <- mcode ("PFOApbtk.iv", PFOAPBTK.iv.code)

########oral#######

## Define the prediction function
pred.PFOA <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                                   ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for oral exposure to 1 mg/kg
  
  BW          = 33                                        # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEoral.L   = 1000                                    # Single oral dose
  #PDOSEoral.L   = 5000                                   #Choose this for high level dosage
  DOSEoral.L    = PDOSEoral.L*BW                          # Amount of oral dose
  ex.oral.L   <- ev (ID = 1, amt= DOSEoral.L, 
                     ii = tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  
  
  ## Get a prediction
  out.A.L <- 
    mod.oral %>%                                                # model object
    param(pars) %>%                                             # to update the parameters in the model subject
    Req(Plasma, Liver, Kidney, Lung,  Skin)%>%                  # select model output
    update(atol = 1E-70, maxsteps=50000) %>%                    # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.L, tgrid=tsamp)                     # Set up the simulation run
  out.A.L<-cbind.data.frame(Time=out.A.L$time/24, 
                            CA=out.A.L$Plasma,
                            CL=out.A.L$Liver,
                            CK=out.A.L$Kidney,
                            CLung=out.A.L$Lung,
                            CSkin=out.A.L$Skin)
  return(list("out.A.L"=out.A.L
  ))
  
}
## The calibrated value using levenberg-marquart algorithm. The code are available from Chou et al.,(2019).
pars.oral <- log(c(
  KbileC                = 0.03223453,                      
  KurineC               = 5.111845,
  Free                  = 2.104425e-02,                        
  PL                    = 1.02858979,                        
  PK                    = 0.04087181,                        
  PLung                 = 0.17605242,                     
  PDer                 = 0.19532406,                       
  PRest                 = 0.25305482,                         
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  K0C                   = 9.193679e-02,               
  Kabsc                 = 2.195354e+1,                       
  KunabsC               = 1.360926e-11                    
))

Sim.A.L = pred.PFOA (pars.oral)$out.A.L 

##Simulated results
dfo.sim.A.L = cbind.data.frame (Time=Sim.A.L$Time, 
                                CA=Sim.A.L$CA, 
                                CL=Sim.A.L$CL,
                                CK=Sim.A.L$CK,
                                CLung=Sim.A.L$CLung,
                                CSkin=Sim.A.L$CSkin)

#########nasal#########

## Define the prediction function
pred.PFOA <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                                        ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for nasal exposure to 1 mg/kg
  
  BW          = 33                                             # Mouse body weight
  tinterval   = 24                                             # Time intervBL
  TDoses      = 1                                              # TotaL dosing/Dose times
  PDOSEnasal.L   = 1000                                         # Single nasal dose
  #PDOSEnasal.L   = 5000                                      #Choose this for high level dosage
  DOSEnasal.L    = PDOSEnasal.L*BW                               # Amount of nasal dose
  ex.nasal.L   <- ev (ID = 1, amt= DOSEnasal.L, 
                      ii = tinterval, addl=TDoses-1, cmt="ALung", replicate = FALSE)
  
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  
  
  ## Get a prediction
  out.B.L <- 
    mod.nasal %>%                                               # model object
    param(pars) %>%                                            # to update the parameters in the model subject
    Req(Plasma, Liver, Kidney, Skin)%>%                        # select model output
    update(atol = 1E-70, maxsteps=50000) %>%                   # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.nasal.L, tgrid=tsamp)                    # Set up the simulation run
  out.B.L<-cbind.data.frame(Time=out.B.L$time/24, 
                            CA=out.B.L$Plasma,
                            CL=out.B.L$Liver,
                            CK=out.B.L$Kidney,
                            CSkin=out.B.L$Skin)
  
  return(list("out.B.L"=out.B.L
  ))
  
}
## The calibrated value using levenberg-marquart algorithm. The code are available from Chou et al.,(2019).
pars.nasal <- log(c(
  KbileC                = 0.03223453,                       
  KurineC               = 5.111845,
  Free                  = 2.104425e-02,                       
  PL                    = 0.82907160,                       
  PK                    = 0.02937613,                         
  PLung                 = 9.832482e-03,                       
  PDer                  = 0.15810824,                       
  PRest                 = 0.35949095,                        
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  PLungA                = 4.268459e+03
))
Sim.B.L = pred.PFOA (pars.nasal)$out.B.L        

##Simulated results
dfo.Sim.B.L = cbind.data.frame (Time=Sim.B.L$Time, 
                                CA=Sim.B.L$CA, 
                                CL=Sim.B.L$CL,
                                CK=Sim.B.L$CK,
                                CSkin=Sim.B.L$CSkin)


#######dermal######

## Define the prediction function
pred.PFOA <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                                   ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for dermal exposue to 1 mg/kg
  BW          = 33                                        # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEdermal.L = 1000*0.083                              # Single dermal dose
  #PDOSEdermal.L   = 5000*0.0719                           #Choose this for high level dosage
  DOSEdermal.L  = PDOSEdermal.L*BW                        # Amount of dermal dose
  ex.dermal.L   <- ev(ID = 1, amt= DOSEdermal.L,              
                      ii = tinterval, addl=TDoses-1, cmt="ACham", replicate = FALSE)
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)  
  ## Get a prediction  
  out.C.L <- 
    mod.dermal %>%                                          # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(Plasma, Liver, Kidney, Lung)%>%                   # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.dermal.L, tgrid=tsamp)               # Set up the simulation run
  out.C.L<-cbind.data.frame(Time=out.C.L$time/24, 
                            CA=out.C.L$Plasma,
                            CL=out.C.L$Liver,
                            CK=out.C.L$Kidney,
                            CLung=out.C.L$Lung)
  
  
  return(list("out.C.L"=out.C.L
  ))
  
}
## The calibrated value using levenberg-marquart algorithm. The code are available from Chou et al.,(2019).
pars.dermal <- log(c(
  KbileC                = 0.03223453,                       
  KurineC               = 5.111845,                        
  Free                  = 2.104425e-02,                       
  PL                    = 1.11999145,                       
  PK                    = 0.04465616,                         
  PLung                 = 0.17447577,                       
  PDer                 = 0.32097265,                        
  PRest                 = 0.33042595,                         
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  PDerCha                = 5.26928447,
  Kvol                  = 6.993407e-07,
  Kp                    = 1.471556e-03,
  SA                    = 0.164
))

Sim.C.L = pred.PFOA (pars.dermal)$out.C.L         

##Simulated results
dfo.Sim.C.L = cbind.data.frame (Time=Sim.C.L$Time, 
                                CA=Sim.C.L$CA, 
                                CL=Sim.C.L$CL,
                                CK=Sim.C.L$CK,
                                CLung=Sim.C.L$CLung)

######IV######

## Define the prediction function
pred.PFOA <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                                   ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for IV exposue to 1 mg/kg
  
  BW          = 33                                        # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEiv.L = 1000                                        # Single IV dose 
  #PDOSEdermal.L   = 5000                                 #Choose this for high level dosage
  DOSEiv.L = PDOSEiv.L*BW                                 # Amount of IV dose
  ex.iv.L <- ev(ID = 1, amt= DOSEiv.L, 
                ii = tinterval, addl=TDoses-1, cmt="APlas_free", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  ## Get a prediction
  out.D.L <- 
    mod.iv %>%                                            # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(Plasma, Liver, Kidney, Lung, Skin)%>%             # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.iv.L, tgrid=tsamp)                 # Set up the simulation run
  out.D.L<-cbind.data.frame(Time=out.D.L$time/24, 
                            CA=out.D.L$Plasma,
                            CL=out.D.L$Liver,
                            CK=out.D.L$Kidney,
                            CLung=out.D.L$Lung,
                            CSkin=out.D.L$Skin)
  return(list(
    "out.D.L"=out.D.L))
  
}
## The calibrated value using levenberg-marquart algorithm. The code are available from Chou et al.,(2019).
pars.IV <- log(c(
  KbileC                = 0.03223453,                       
  KurineC               = 5.111845,                        
  Free                  = 2.104425e-02,                        
  PL                    = 0.92283219,                        
  PK                    = 0.04401475,                        
  PLung                 = 1.367528e-01,                      
  PDer                 = 0.13323227,                       
  PRest                 = 0.35531232,                         
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04
))
Sim.D.L = pred.PFOA (pars.IV)$out.D.L         

##Simulated results
dfo.sim.D.L = cbind.data.frame (Time=Sim.D.L$Time, 
                                CA=Sim.D.L$CA, 
                                CL=Sim.D.L$CL,
                                CK=Sim.D.L$CK,
                                CLung=Sim.D.L$CLung,
                                CSkin=Sim.D.L$CSkin)

##########NSC for oral#########

pred.PFOA.oral <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Define the four exposure scenario
  ## Exposure scenario for oral exposue to 1 mg/kg
  
  BW          = 33                                     # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEoral.L   = 1000                                         # Single oral dose
  #PDOSEoral.L   = 5000                                      #Choose this for NSC of high level dosage
  DOSEoral.L    = PDOSEoral.L*BW                              # Amount of oral dose
  ex.oral.L   <- ev (ID = 1, amt= DOSEoral.L, 
                     ii = tinterval, addl=TDoses-1, cmt="AST", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  
  
  ## Get a prediction
  out.A.L <- 
    mod.oral %>%                                               # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(AUC_CA, AUC_CL, AUC_CK, AUC_CLung, AUC_CSkin, AUC_CST, AUC_CSI)%>%                        # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.oral.L, tgrid=tsamp)%>%             # Set up the simulation run
    filter(time!=0)                     #The code can produce time-dependent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  out.A.L<-cbind.data.frame(Time=out.A.L$time/24, 
                            AUC_CA=out.A.L$AUC_CA,
                            AUC_CL=out.A.L$AUC_CL,
                            AUC_CK=out.A.L$AUC_CK,
                            AUC_CLung=out.A.L$AUC_CLung,
                            AUC_CSkin=out.A.L$AUC_CSkin,
                            AUC_CST=out.A.L$AUC_CST,
                            AUC_CSI=out.A.L$AUC_CSI)
  return(list("out.A.L"=out.A.L))
  
}

pars.oral <- log(c(
  KbileC                = 0.03223453,                      
  #KurineC               = 1.111845e+00,                         
  KurineC               =  5.111845,
  Free                  = 2.104425e-02,                       
  PL                    = 1.02858979,                        
  PK                    = 0.04087181,                         
  PLung                 = 0.17605242,                       
  PDer                 = 0.19532406,                        
  PRest                 = 0.25305482,                        
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  K0C                   = 9.193679e-02,                
  Kabsc                 = 2.195354e+1,                        
  KunabsC               = 1.360926e-11,                     
  BW=33, QCC=16.5, QLC=0.161, QKC=0.091, QLungC=0.005, QSkinC=0.058, QFilC=0.045, 
  Htc=0.48, VLC=0.0589, VKC=0.0169, VPlasC=0.049, VFilC=0.0017, VLungC=0.00769, 
  VSkinC=0.1653, VSTC=0.0075, VSIC=0.038, Qp=6.99, GEC=0.188, VCham=0.00005
))
Sim.oral = pred.PFOA.oral (pars.oral)

## Create the matrix for normalized sensitivity coefficient data
NSC_oral   = matrix(nrow=32,ncol=7)


for (i in 1:32) {
  pars.oral.new      <- log(c(exp(pars.oral[i])*1.01,exp(pars.oral[-i]))) # Each cycle, generate a new value of parameter i (e.g., 10.0a), and delete parameter i, so that you can proceed to the next parameter i+1
  Simnew.oral             <- pred.PFOA.oral (pars.oral.new)
  delta.pars.oral       <- exp(pars.oral[i])/(exp(pars.oral[i])*0.01) # Basically, the ratio is 100 for each parameter.
  
  ## Estimated the AUC
  AUC.CA.new   =  Simnew.oral$out.A.L %>% filter (Time == 32) %>% select (AUC_CA) # AUC_CA, AUC_CL, and AUC_CK were defined in the main mrgsolve code
  AUC.CA.ori   =  Sim.oral$out.A.L %>% filter (Time == 32) %>% select (AUC_CA)
  AUC.CL.new   =  Simnew.oral$out.A.L %>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CL.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CK.new   =  Simnew.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CK.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CLung.new   =  Simnew.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CLung.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CSkin.new   =  Simnew.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CSkin.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CST.new   =  Simnew.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CST.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CSI.new   =  Simnew.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  AUC.CSI.ori   =  Sim.oral$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  
  delta.AUC.CA    =  AUC.CA.new - AUC.CA.ori
  delta.AUC.CL    =  AUC.CL.new - AUC.CL.ori
  delta.AUC.CK    =  AUC.CK.new - AUC.CK.ori
  delta.AUC.CLung    =  AUC.CLung.new - AUC.CLung.ori
  delta.AUC.CSkin    =  AUC.CSkin.new - AUC.CSkin.ori
  delta.AUC.CST    =  AUC.CST.new - AUC.CST.ori
  delta.AUC.CSI    =  AUC.CSI.new - AUC.CSI.ori
  
  NSC_oral   [i, 1]   <- as.numeric((delta.AUC.CA/AUC.CA.ori) * delta.pars.oral)
  NSC_oral   [i, 2]   <- as.numeric((delta.AUC.CL/AUC.CL.ori) * delta.pars.oral)
  NSC_oral   [i, 3]   <- as.numeric((delta.AUC.CK/AUC.CK.ori) * delta.pars.oral)
  NSC_oral   [i, 4]   <- as.numeric((delta.AUC.CLung/AUC.CLung.ori) * delta.pars.oral)
  NSC_oral   [i, 5]   <- as.numeric((delta.AUC.CSkin/AUC.CSkin.ori) * delta.pars.oral)
  NSC_oral   [i, 6]   <- as.numeric((delta.AUC.CST/AUC.CST.ori) * delta.pars.oral)
  NSC_oral   [i, 7]   <- as.numeric((delta.AUC.CSI/AUC.CSI.ori) * delta.pars.oral)
  
}


colnames (NSC_oral)  = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK","NSC_AUC_CLung","NSC_AUC_CSkin","NSC_AUC_CST","NSC_AUC_CSI") 
rownames(NSC_oral)   = c("KbileC",  "KurineC", "Free", "PL", "PK", "PLung", "PDer",
                         "PRest", "K0C", "Kabsc", "KunabsC", "Tmc", "Kt", "BW", "QCC", "QLC", "QKC", "QLungC", "QSkinC", "QFilC", 
                         "Htc", "VLC", "VKC", "VPlasC", "VFilC", "VLungC", "VSkinC", "VSTC", "VSIC", "Qp", "GEC", "VCham")

##########NSC for nasal#######

pred.PFOA.nasal <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for nasal exposue to 1 mg/kg
  
  BW          = 33                                     # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEnasal.L   = 1000                                         # Single nasal dose
  #PDOSEnasal.L   = 5000                                      #Choose this for NSC of high level dosage
  DOSEnasal.L    = PDOSEnasal.L*BW                              # Amount of nasal dose
  ex.nasal.L   <- ev (ID = 1, amt= DOSEnasal.L, 
                      ii = tinterval, addl=TDoses-1, cmt="ALung", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  
  
  ## Get a prediction
  ## Out A: nasal exposure to 1 mg/kg-d, matrix: Plasma, liver, kidney, lung, stomach, intestine, skin
  out.A.L <- 
    mod.nasal %>%                                               # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(AUC_CA, AUC_CL, AUC_CK, AUC_CLung, AUC_CSkin, AUC_CST, AUC_CSI)%>%                        # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.nasal.L, tgrid=tsamp)%>%             # Set up the simulation run
    filter(time!=0)                     #The code can produce time-dependent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  out.A.L<-cbind.data.frame(Time=out.A.L$time/24, 
                            AUC_CA=out.A.L$AUC_CA,
                            AUC_CL=out.A.L$AUC_CL,
                            AUC_CK=out.A.L$AUC_CK,
                            AUC_CLung=out.A.L$AUC_CLung,
                            AUC_CSkin=out.A.L$AUC_CSkin,
                            AUC_CST=out.A.L$AUC_CST,
                            AUC_CSI=out.A.L$AUC_CSI)
  return(list("out.A.L"=out.A.L))
  
}

pars.nasal <- log(c(
  KbileC                = 0.03223453,                       
  KurineC               =  5.111845,
  Free                  = 2.104425e-02,                        
  PL                    = 0.82907160,                        
  PK                    = 0.02937613,                         
  PLung                 = 9.832482e-03,                      
  PDer                 = 0.15810824,                       
  PRest                 = 0.35949095,                        
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  PLungA                = 4.268459e+03,
  BW=33, QCC=16.5, QLC=0.161, QKC=0.091, QLungC=0.005, QSkinC=0.058, QFilC=0.045, 
  Htc=0.48, VLC=0.0589, VKC=0.0169, VPlasC=0.049, VFilC=0.0017, VLungC=0.00769, 
  VSkinC=0.1653, VSTC=0.0075, VSIC=0.038, QpC=69.9, GEC=0.188, VCham=0.00005
))
Sim.nasal = pred.PFOA.nasal (pars.nasal)

## Create the matrix for normalized sensitivity coefficient data
NSC_nasal   = matrix(nrow=30,ncol=7)


for (i in 1:30) {
  pars.nasal.new      <- log(c(exp(pars.nasal[i])*1.01,exp(pars.nasal[-i]))) # Each cycle, generate a new value of parameter i (e.g., 10.0a), and delete parameter i, so that you can proceed to the next parameter i+1
  Simnew.nasal             <- pred.PFOA.nasal (pars.nasal.new)
  delta.pars.nasal       <- exp(pars.nasal[i])/(exp(pars.nasal[i])*0.01) # Basically, the ratio is 100 for each parameter.
  
  ## Estimated the AUC
  AUC.CA.new   =  Simnew.nasal$out.A.L %>% filter (Time == 32) %>% select (AUC_CA) 
  AUC.CA.ori   =  Sim.nasal$out.A.L %>% filter (Time == 32) %>% select (AUC_CA)
  AUC.CL.new   =  Simnew.nasal$out.A.L %>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CL.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CK.new   =  Simnew.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CK.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CLung.new   =  Simnew.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CLung.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CSkin.new   =  Simnew.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CSkin.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CST.new   =  Simnew.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CST.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CSI.new   =  Simnew.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  AUC.CSI.ori   =  Sim.nasal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  
  delta.AUC.CA    =  AUC.CA.new - AUC.CA.ori
  delta.AUC.CL    =  AUC.CL.new - AUC.CL.ori
  delta.AUC.CK    =  AUC.CK.new - AUC.CK.ori
  delta.AUC.CLung    =  AUC.CLung.new - AUC.CLung.ori
  delta.AUC.CSkin    =  AUC.CSkin.new - AUC.CSkin.ori
  delta.AUC.CST    =  AUC.CST.new - AUC.CST.ori
  delta.AUC.CSI    =  AUC.CSI.new - AUC.CSI.ori
  
  NSC_nasal   [i, 1]   <- as.numeric((delta.AUC.CA/AUC.CA.ori) * delta.pars.nasal)
  NSC_nasal   [i, 2]   <- as.numeric((delta.AUC.CL/AUC.CL.ori) * delta.pars.nasal)
  NSC_nasal   [i, 3]   <- as.numeric((delta.AUC.CK/AUC.CK.ori) * delta.pars.nasal)
  NSC_nasal   [i, 4]   <- as.numeric((delta.AUC.CLung/AUC.CLung.ori) * delta.pars.nasal)
  NSC_nasal   [i, 5]   <- as.numeric((delta.AUC.CSkin/AUC.CSkin.ori) * delta.pars.nasal)
  NSC_nasal   [i, 6]   <- as.numeric((delta.AUC.CST/AUC.CST.ori) * delta.pars.nasal)
  NSC_nasal   [i, 7]   <- as.numeric((delta.AUC.CSI/AUC.CSI.ori) * delta.pars.nasal)
  
}


colnames (NSC_nasal)  = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK","NSC_AUC_CLung","NSC_AUC_CSkin","NSC_AUC_CST","NSC_AUC_CSI") 
rownames(NSC_nasal)   = c("KbileC",  "KurineC", "Free", "PL", "PK", "PLung", "PDer",
                          "PRest", "PLungA",  "Tmc", "Kt","BW", "QCC", "QLC", "QKC", "QLungC", "QSkinC", "QFilC", 
                          "Htc", "VLC", "VKC", "VPlasC", "VFilC", "VLungC", "VSkinC", "VSTC", "VSIC", "QpC", "GEC", "VCham")

########NSC for dermal######

pred.PFOA.dermal <- function(pars) {
  
  ## Get out of log domain
  pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  
  ## Exposure scenario for dermal exposue to 1 mg/kg
  
  BW          = 33                                     # Mouse body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  PDOSEdermal.L   = 1000                                         # Single dermal dose
  #PDOSEdermal.L   = 5000                                    #Choose this for NSC of high level dosage
  DOSEdermal.L    = PDOSEdermal.L*BW                              # Amount of dermal dose
  ex.dermal.L   <- ev (ID = 1, amt= DOSEdermal.L, 
                       ii = tinterval, addl=TDoses-1, cmt="ACham", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0,tinterval*(TDoses-1)+24*40,1)
  
  
  ## Get a prediction
  out.A.L <- 
    mod.dermal %>%                                               # model object
    param(pars) %>%                                       # to update the parameters in the model subject
    Req(AUC_CA, AUC_CL, AUC_CK, AUC_CLung, AUC_CSkin, AUC_CST, AUC_CSI)%>%                        # select model output
    update(atol = 1E-70, maxsteps=50000) %>%              # solver setting, atol: Absolute tolerance parameter
    mrgsim_d(data = ex.dermal.L, tgrid=tsamp)%>%             # Set up the simulation run
    filter(time!=0)                     #The code can produce time-dependent NSC values, but at time = 0, NSC cannot be calculated, so data at time = 0 needs to be filtered out.
  
  out.A.L<-cbind.data.frame(Time=out.A.L$time/24, 
                            AUC_CA=out.A.L$AUC_CA,
                            AUC_CL=out.A.L$AUC_CL,
                            AUC_CK=out.A.L$AUC_CK,
                            AUC_CLung=out.A.L$AUC_CLung,
                            AUC_CSkin=out.A.L$AUC_CSkin,
                            AUC_CST=out.A.L$AUC_CST,
                            AUC_CSI=out.A.L$AUC_CSI)
  return(list("out.A.L"=out.A.L))
  
}

pars.dermal <- log(c(
  KbileC                = 0.03223453,                       
  KurineC               = 5.111845,                         
  Free                  = 2.104425e-02,                       
  PL                    = 1.11999145,                       
  PK                    = 0.04465616,                        
  PLung                 = 0.17447577,                       
  PDer                 = 0.32097265,                       
  PRest                 = 0.33042595,                         
  Tmc                   = 2.242729e+04,
  Kt                    = 7.344000e+04,
  PDerCha                = 5.26928447,
  Kvol                  = 6.993407e-07,
  Kp                    = 1.471556e-03,
  SA                    = 0.164,
  BW=33, QCC=16.5, QLC=0.161, QKC=0.091, QLungC=0.005, QSkinC=0.058, QFilC=0.045, 
  Htc=0.48, VLC=0.0589, VKC=0.0169, VPlasC=0.049, VFilC=0.0017, VLungC=0.00769, 
  VSkinC=0.1653, VSTC=0.0075, VSIC=0.038, Qp=6.99, GEC=0.188, VCham=0.00005
))
Sim.dermal = pred.PFOA.dermal (pars.dermal)

## Create the matrix for normalized sensitivity coefficient data
NSC_dermal   = matrix(nrow=33,ncol=7)


for (i in 1:33) {
  pars.dermal.new      <- log(c(exp(pars.dermal[i])*1.01,exp(pars.dermal[-i]))) # Each cycle, generate a new value of parameter i (e.g., 10.0a), and delete parameter i, so that you can proceed to the next parameter i+1
  Simnew.dermal             <- pred.PFOA.dermal (pars.dermal.new)
  delta.pars.dermal       <- exp(pars.dermal[i])/(exp(pars.dermal[i])*0.01) # Basically, the ratio is 100 for each parameter.
  
  ## Estimated the AUC
  AUC.CA.new   =  Simnew.dermal$out.A.L %>% filter (Time == 32) %>% select (AUC_CA) # AUC_CA, AUC_CL, and AUC_CK were defined in the main mrgsolve code
  AUC.CA.ori   =  Sim.dermal$out.A.L %>% filter (Time == 32) %>% select (AUC_CA)
  AUC.CL.new   =  Simnew.dermal$out.A.L %>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CL.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CL)
  AUC.CK.new   =  Simnew.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CK.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CK)
  AUC.CLung.new   =  Simnew.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CLung.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CLung)
  AUC.CSkin.new   =  Simnew.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CSkin.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSkin)
  AUC.CST.new   =  Simnew.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CST.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CST)
  AUC.CSI.new   =  Simnew.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  AUC.CSI.ori   =  Sim.dermal$out.A.L%>% filter (Time == 32) %>% select (AUC_CSI)
  
  delta.AUC.CA    =  AUC.CA.new - AUC.CA.ori
  delta.AUC.CL    =  AUC.CL.new - AUC.CL.ori
  delta.AUC.CK    =  AUC.CK.new - AUC.CK.ori
  delta.AUC.CLung    =  AUC.CLung.new - AUC.CLung.ori
  delta.AUC.CSkin    =  AUC.CSkin.new - AUC.CSkin.ori
  delta.AUC.CST    =  AUC.CST.new - AUC.CST.ori
  delta.AUC.CSI    =  AUC.CSI.new - AUC.CSI.ori
  
  NSC_dermal   [i, 1]   <- as.numeric((delta.AUC.CA/AUC.CA.ori) * delta.pars.dermal)
  NSC_dermal   [i, 2]   <- as.numeric((delta.AUC.CL/AUC.CL.ori) * delta.pars.dermal)
  NSC_dermal   [i, 3]   <- as.numeric((delta.AUC.CK/AUC.CK.ori) * delta.pars.dermal)
  NSC_dermal   [i, 4]   <- as.numeric((delta.AUC.CLung/AUC.CLung.ori) * delta.pars.dermal)
  NSC_dermal   [i, 5]   <- as.numeric((delta.AUC.CSkin/AUC.CSkin.ori) * delta.pars.dermal)
  NSC_dermal   [i, 6]   <- as.numeric((delta.AUC.CST/AUC.CST.ori) * delta.pars.dermal)
  NSC_dermal   [i, 7]   <- as.numeric((delta.AUC.CSI/AUC.CSI.ori) * delta.pars.dermal)
  
}


colnames (NSC_dermal)  = c("NSC_AUC_CA","NSC_AUC_CL","NSC_AUC_CK","NSC_AUC_CLung","NSC_AUC_CSkin","NSC_AUC_CST","NSC_AUC_CSI") 
rownames(NSC_dermal)   = c("KbileC",  "KurineC", "Free", "PL", "PK", "PLung", "PDer", "PDerCha",
                           "PRest",  "Tmc", "Kt",  "Kvol", "Kp", "SA","BW", "QCC", "QLC", "QKC", "QLungC", "QSkinC", "QFilC", 
                           "Htc", "VLC", "VKC", "VPlasC", "VFilC", "VLungC", "VSkinC", "VSTC", "VSIC", "Qp", "GEC", "VCham")

#########plot function######
Circle.plot <- function (melt.data){ # melt.data is an argument of Circle.plot function.
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=3
  to_add = data.frame(matrix(NA, empty_bar*nlevels(as.factor(melt.data$group)), ncol(melt.data)) )
  colnames(to_add) = colnames(melt.data)
  to_add$group=rep(levels(as.factor(melt.data$group)), each=empty_bar)
  melt.data=rbind(melt.data, to_add)
  melt.data=melt.data %>% arrange(group)
  melt.data$id=seq(1, nrow(melt.data)) # id is the number of rows. In total, there were 68 rows.
  
  # Get the name and the y position of each label
  label_data=melt.data
  number_of_bar=nrow(label_data) 
  angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust<-ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data=melt.data %>%
    group_by(group) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start = grid_data$start - 1
  grid_data=grid_data[-1,]
  
  # Make the plot
  windowsFonts(Times=windowsFont("Times New Roman"))
  
  p.cir.plot <- 
    ggplot(melt.data, aes(x=as.factor(id), y = abs(value*100), fill = group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
    scale_fill_brewer(palette = "Dark2")+
    # Add a val=80/60/40/20 lines. I do it at the beginning to make sure barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 160, xend = start, yend = 160), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 140, xend = start, yend = 140), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 120, xend = start, yend = 120), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=0.5, size=0.8 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 20/40/60/80/100/120/140/160 lines
    annotate("text", x = rep(max(melt.data$id),8), y = c(20, 40, 60, 80, 100, 120, 140, 160), label = c("20%", "40%", "60%", "80%", "100%", "120%", "140%", "160%") , color="black", size=3 , angle=0, fontface="bold", hjust=1) +
    geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.9) +
    ylim(-100,200) +
    theme_minimal() +
    theme(
      legend.position         = "none",
      text                    = element_text (family = "Times"),
      panel.background        = element_blank (),
      plot.background         = element_blank (),
      axis.text               = element_blank(),
      axis.title              = element_blank(),
      panel.grid              = element_blank(),
      plot.margin             = unit(rep(-2,4), "cm")
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=id, y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha = 1, size = 3.5, angle= label_data$angle, inherit.aes = FALSE) +
    
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=1.0 , inherit.aes = FALSE )  +
    geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(1,0.5,0), colour = "black", alpha=0.8, size= 4, fontface="bold", inherit.aes = FALSE)
  
  return (p.cir.plot)
}

###################### Figure ####################
melt.oral.CA        = melt(NSC_oral[,1]) 
melt.oral.CA$group  = c("Oral") 
melt.nasal.CA          = melt(NSC_nasal[,1])
melt.nasal.CA$group    = c("Nasal")
melt.dermal.CA       = melt(NSC_dermal[,1])
melt.dermal.CA$group = c("Dermal")

melt.data.CA         = rbind (melt.oral.CA,melt.nasal.CA,melt.dermal.CA)
melt.data.CA$par     = c(rownames(melt.oral.CA),rownames(melt.nasal.CA),rownames(melt.dermal.CA)) 

p1.CA                = Circle.plot (melt.data.CA%>%filter(abs(value)>0.1))
p1.CA