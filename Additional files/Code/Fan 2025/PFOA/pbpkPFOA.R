if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here())

rm(list=ls())
library(openxlsx)
library(stringr)
library(ggplot2)
library(sp)
library(mrgsolve)
library(magrittr)
library(dplyr)
library(FME)
library(truncnorm)
library(invgamma)
library(foreach)
library(doParallel)
library(ggplot2)
library(EnvStats)
library(plyr)
library(dplyr)


# The permeability-limited PBPK model for male rats
PBPK <- "
$PARAM @annotated
BW       : 3e5    : mg
PVB      : 5.4e-2   : 5.4/100 unitless
PVplasma : 3.12e-2  : 3.12/100 unitless
PVK      : 7.3e-3   : 0.73/100 unitless
PVKB     : 0.16     : unitless
PVKF     : 0.13     : unitless
VFil     : 2.5e2    : 0.25*1e3 uL
PVL      : 0.0366   : 3.66/100 unitless
PVLB     : 0.21     : unitless
PVLF     : 0.049    : unitless
PVbile   : 0.004    : unitless
PVG      : 0.0269   : 2.69/100 unitless
PVGB     : 0.034    : unitless
PVGF     : 0.28     : unitless
PVGL     : 0.045    : 4.5/100 unitless
PVM      : 0.4043   : 40.43/100 unitless
PVMB     : 0.04     : unitless
PVMF     : 0.054    : unitless
PVA      : 0.07     : 7.0/100 unitless
PVAB     : 0.02     : unitless
PVAF     : 0.174    : unitless
PVRB     : 0.036    : unitless
PVRF     : 0.18     : unitless
PAK      : 35       : mm2/mg BW
n        : 5        : unitless
PAKG     : 6.89     : mm2/mg BW
PAL      : 25       : mm2/mg BW
PAG      : 10       : mm2/mg BW
PAGL     : 4.14     : mm2/mg BW
PAM      : 7        : mm2/mg BW
PAA      : 7        : mm2/mg BW
PAR      : 10       : mm2/mg BW
PQBK     : 0.141    : 14.1/100 unitless
PQBG     : 0.151    : 15.1/100 unitless
PQBL     : 0.024    : 2.4/100 unitless
PQBM     : 0.278    : 27.8/100 unitless
PQBA     : 0.07     : 7.0/100 unitless
Qfeces   : 234.58   : 5.63*1e3/24 uL/h
PQbile   : 3.75e-3  : 90.0*1e3/24/1e6 uL/h/mg BW
PQGFR    : 0.6444   : 10.74*1e3*60/1e6 uL/h/mg BW
PQurine  : 8.33e-3  : 200.0*1e3/24/1e6 uL/h/mg BW
PeffB    : 0.1793   : 4.98e-5*3600 mm/h permeability
PeffK    : 0.1577   : 4.38e-5*3600 mm/h permeability
PeffL    : 0.1854   : 5.15e-5*3600 mm/h permeability
PeffG    : 9.54e-2  : 2.65e-5*3600 mm/h permeability
PeffA    : 9.54e-2  : 2.65e-5*3600 mm/h permeability
PeffM    : 9.54e-2  : 2.65e-5*3600 mm/h permeability
PeffR    : 9.54e-2  : 2.65e-5*3600 mm/h permeability
CRssG    : 3.75     : cell-water concentration ratio 
CRssL    : 7.28     : cell-water concentration ratio
CRssK    : 6.19     : cell-water concentration ratio
Pbclear  : 0.9936   : 2.76e-4*3600 mm/h renal clearance rate
Pbreab   : 0.4248   : 1.18e-4*3600 mm/h renal reabsorption rate
Pbabs    : 0.6408   : 1.78e-4*3600 mm/h absorption rate of liver
Pbefflux : 0.4968   : 1.38e-4*3600 mm/h renal efflux rate
CalbB    : 281e-3   : unit mol/m3
CalbKF   : 243e-3   : unit mol/m3
CalbLF   : 243e-3   : unit mol/m3
CalbGF   : 146e-3   : unit mol/m3
CalbMF   : 146e-3   : unit mol/m3
CalbAF   : 73e-3    : unit mol/m3
CalbRF   : 73e-3    : unit mol/m3
Ca2uKT   : 110e-3   : unit mol/m3
CLfabpKT : 2.65e-3  : unit mol/m3
CLfabpLT : 133e-3   : unit mol/m3
Ka       : 32.894    : 3.1*10.578 m3/mol multiplying by number of binding sites
KLfabp   : 81.260      : 81.260 m3/mol
Ka2u     : 0.5      : m3/mol

$MAIN
// parameter manipulation
double VB = PVB * BW;
double Vplasma = PVplasma * BW;
double VK = PVK * BW;
double VKB = PVKB * PVK * BW;
double VKF = PVKF * PVK * BW;
double VKT = VK - VKF;
double VL = PVL * BW;
double VLB = PVLB * PVL * BW;
double VLF = PVLF * PVL * BW;
double VLT = VL - VLF;
double Vbile = PVbile * PVL * BW;
double VG = PVG * BW;
double VGB = PVGB * PVG * BW;
double VGF = PVGF * PVG * BW;
double VGT = VG - VGF;
double VGL = PVGL * BW;
double VM = PVM * BW;
double VMB = PVMB * PVM * BW;
double VMF = PVMF * PVM * BW;
double VMT = VM - VMF;
double VA = PVA * BW;
double VAB = PVAB * PVA * BW;
double VAF = PVAF * PVA * BW;
double VAT = VA - VAF;
double PVR = 1.0 - PVB - PVK - PVL - PVG - PVM - PVA;
double VR = PVR * BW;
double VRB = PVRB * PVR * BW;
double VRF = PVRF * PVR * BW;
double VRT = VR - VRF;

double AK = PAK * VK;
double AKG = PAKG * VK;
double AL = PAL * VL;
double AG = PAG * VG;
double AGL = PAGL * BW;
double AM = PAM * VM;
double AA = PAA * VA;
double AR = PAR * VR;

double Qcardiac = 0.235*1e6*60 * pow((BW / 1e6), 0.75); // unit uL/h
double QBK = PQBK * Qcardiac;
double QBG = PQBG * Qcardiac;
double QBL = (PQBL + PQBG) * Qcardiac;
double QBM = PQBM * Qcardiac;
double QBA = PQBA * Qcardiac;
double PQBR = 1.0 - PQBK - PQBG - PQBL - PQBM - PQBA;
double QBR = PQBR * Qcardiac;
double Qbile = PQbile * BW;
double QGFR = PQGFR * BW;
double Qurine = PQurine * BW;

// operations on parameters to generate coefficients for differential equation
double kBKF = 1. / ((1. / QBK) + 1./(PeffB * AK));
double kBF =  PeffB * AKG;
double kKFKT = PeffK * AK;
double kFKT = PeffK * AK * n;
double kBLF = 1. / ((1./QBL) + 1./(PeffB * AL));
double kLFLT = PeffL * AL;
double kBGF = 1. / ((1./QBG) + 1./(PeffB * AG));
double kGFGT = PeffG * AG;
double kGLGT = PeffG * AGL;
double kBMF = 1. / ((1./QBM) + 1./(PeffB * AM));
double kMFMT = PeffM * AM;
double kBAF = 1. / ((1./QBA) + 1./(PeffB * AA));
double kAFAT = PeffA * AA;
double kBRF = 1. / ((1./QBR) + 1./(PeffB * AR));
double kRFRT = PeffR * AR;
double kbileLT = PeffL * AL;

double bBKF = kBKF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bKFB = kBKF / VKF;
double bKFKT = kKFKT / VKF;
double bKTKF = kKFKT / VKT;
double bFKT = kFKT / VFil;
double bKTF = kFKT / (VKT * CRssK);
double bBF = QGFR / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bFB = kBF / VFil;

double bBLF = kBLF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bLFB = kBLF / VLF;
double bLFLT = kLFLT / VLF;
double bLTLF = kLFLT / VLT;
double bbileLT = kbileLT / Vbile;
double bLTbile = kbileLT / (VLT * CRssL);

double bBGF = kBGF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bGFB = kBGF / VGF;
double bGFGT = kGFGT / VGF;
double bGTGF = kGFGT / VGT;
double bGLGT = kGLGT / VGL;
double bGTGL = kGLGT / (VGT * CRssG);

double bBMF = kBMF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bMFB = kBMF / VMF;
double bMFMT = kMFMT / VMF;
double bMTMF = kMFMT / VMT;

double bBAF = kBAF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bAFB = kBAF / VAF;
double bAFAT = kAFAT / VAF;
double bATAF = kAFAT / VAT;

double bBRF = kBRF / (VB+VLB+VKB+VGB+VMB+VAB+VRB);
double bRFB = kBRF / VRF;
double bRFRT = kRFRT / VRF;
double bRTRF = kRFRT / VRT;

double bclear = Pbclear * AK / VKF;
double breab = n * Pbreab * AK / VFil;
double babs = Pbabs * AL / VLF;
double befflux = Pbefflux * AK / VKT;

$CMT yB yKF yKT yFil yLF yLT yBile yGF yGT yMF yMT yAF yAT yRF yRT yGL

$ODE
double yfreeB = yB * 1.0 / (1.0 + CalbB * Ka);
double yfreeKF = yKF * 1.0 / (1.0 + CalbKF * Ka);
double yfreeLF = yLF * 1.0 / (1.0 + CalbLF * Ka);
double yfreeGF = yGF * 1.0 / (1.0 + CalbGF * Ka);
double yfreeMF = yMF * 1.0 / (1.0 + CalbMF * Ka);
double yfreeAF = yAF * 1.0 / (1.0 + CalbAF * Ka);
double yfreeRF = yRF * 1.0 / (1.0 + CalbRF * Ka);
double yfreeKT = yKT * 1.0 / (1.0 + Ca2uKT * Ka2u + CLfabpKT * KLfabp);
double yfreeLT = yLT * 1.0 / (1.0 + CLfabpLT * KLfabp);

dxdt_yB = bKFB*yfreeKF + bLFB*yfreeLF - bBLF*yfreeB - bBKF*yfreeB - bBF*yfreeB 
+ bFB*yFil + bGFB*yfreeGF - bBGF*yfreeB + bMFB*yfreeMF - bBMF*yfreeB + 
bAFB*yfreeAF - bBAF*yfreeB + bRFB*yfreeRF - bBRF*yfreeB;
    
dxdt_yKF = bBKF*yfreeB - bKFB*yfreeKF + befflux*yfreeKT + bKTKF*yfreeKT - 
bKFKT*yfreeKF - bclear*yfreeKF;
dxdt_yKT = bKFKT*yfreeKF + bFKT*yFil + breab*yFil + bclear*yfreeKF - 
befflux*yfreeKT - bKTKF*yfreeKT - bKTF*yfreeKT;
dxdt_yFil = bKTF*yfreeKT - breab*yFil + bBF*yfreeB - bFB*yFil - bFKT*yFil - 
Qurine*yFil / VFil;

dxdt_yLF = bBLF*yfreeB - bLFB*yfreeLF + bLTLF*yfreeLT - bLFLT*yfreeLF - 
babs*yfreeLF;
dxdt_yLT = bLFLT*yfreeLF + babs*yfreeLF + bbileLT*yBile - bLTbile*yfreeLT - 
bLTLF*yfreeLT;
dxdt_yBile = bLTbile*yfreeLT - bbileLT*yBile - Qbile*yBile / Vbile;

dxdt_yGF = bBGF*yfreeB - bGFB*yfreeGF + bGTGF*yGT - bGFGT*yfreeGF;
dxdt_yGT = bGFGT*yfreeGF - bGTGF*yGT + bGLGT*yGL - bGTGL*yGT;
dxdt_yMF = bBMF*yfreeB - bMFB*yfreeMF + bMTMF*yMT - bMFMT*yfreeMF;
dxdt_yMT = bMFMT*yfreeMF - bMTMF*yMT;
dxdt_yAF = bBAF*yfreeB - bAFB*yfreeAF + bATAF*yAT - bAFAT*yfreeAF;
dxdt_yAT = bAFAT*yfreeAF - bATAF*yAT;
dxdt_yRF = bBRF*yfreeB - bRFB*yfreeRF + bRTRF*yRT - bRFRT*yfreeRF;
dxdt_yRT = bRFRT*yfreeRF - bRTRF*yRT;
dxdt_yGL = bGTGL*yGT - bGLGT*yGL + Qbile*yBile / Vbile - Qfeces*yGL / VGL;

$TABLE
// assume the density of tissue is 1 g/mL or 1 mg/uL, so final unit is ng/g
capture blood = yB / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VB / Vplasma * 1e3;
capture kidney = (yB / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VKB + yKF + yKT) / 
(VKB+VKT+VKF) * 1e3;
capture liver = ((yB) / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VLB + yLF + yLT) / 
(VLB+VLT+VLF) * 1e3;
capture gut = ((yB) / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VGB + yGF + yGT) / 
(VGB+VGT+VGF) * 1e3;
capture muscle = ((yB) / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VMB + yMF + yMT) / 
(VMB+VMT+VMF) * 1e3;
capture adipose = ((yB) / (VB+VLB+VKB+VGB+VMB+VAB+VRB) * VAB + yAF + yAT) / 
(VAB+VAT+VAF) * 1e3;
"
model <- mcode("PBPK", PBPK)

#  model parameters
theta <- log(c(
  PeffB = 0.1793, # population mean
  PeffK = 0.1577,
  PeffL = 0.1854,
  CRssL = 7.28,
  CRssK = 6.19,
  Pbclear = 0.9936,
  Pbreab = 0.4248,
  Pbabs = 0.6408,
  Pbefflux = 0.4968
  
))

# set working directory to base project directory (root relative path)
getwd()


#sig_idx <- grep("sig", names(theta))
## input data set for model calibration/ oral
Obs.A1 <-read.csv(file="Additional files/Code/Fan 2025/PFOA/epa_blood_iv_ex.csv")        
names(Obs.A1) = c("time", "blood")

Obs.B1 <-read.csv(file="Additional files/Code/Fan 2025/PFOA/epa_blood_oral_ex.csv")      
names(Obs.B1) = c("time", "blood")

Obs.B2 <-read.csv(file="Additional files/Code/Fan 2025/PFOA/epa_new_oral_ex.csv")        
names(Obs.B2) = c("time", "blood")

## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred.ratA1 <- function(pars) {
  
  ## Get out of log domain
 # pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  pars %<>% lapply(exp) 
   
  ## Define the three exposure scenario
  ## 1mg/kg/d Exposure scenario for IV 
  
  BW          = 0.29                                       # Rat body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  Ka = 32.894
  KLfabp = 81.260
  
  # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
  ex.oral.A1   <- ev (ID = 1, amt= 290000*1, 
                     ii = tinterval, addl=TDoses-1, cmt="yB", replicate = FALSE)
 
    ## set up the exposure time
  tsamp=tgrid(0, 30*20, 0.5)
  
  
  ## Get a prediction
  ## Out A1: oral exposure to 1 mg/kg-d, matrix: Plasma, urine and feces
  out.A1 <-model %>%
    param(pars) %>%
    Req(blood, kidney, liver) %>%
    update(atol = 1E-50, maxsteps = 50000) %>%
    mrgsim_d(data=ex.oral.A1, tgrid=tsamp) %>%
    filter(time!=0) 
  out.A1 <- cbind.data.frame(time=out.A1$time, 
                                 blood=out.A1$blood,
                                 Kidney=out.A1$kidney,
                                 Liver=out.A1$liver)  
  
  return(list("out.A1"=out.A1 ))
  
}


pred.ratB1 <- function(pars) {
  
  ## Get out of log domain
  # pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  pars %<>% lapply(exp) 
  
  ## Define the three exposure scenario
  ## 1mg/kg/d Exposure scenario for IV 
  
  BW          = 0.29                                       # Rat body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  Ka = 32.894
  KLfabp = 81.260
  
  # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
  ex.oral.B1   <- ev (ID = 1, amt= 290000*1, 
                      ii = tinterval, addl=TDoses-1, cmt="yGL", replicate = FALSE)
  
  ## set up the exposure time
  tsamp=tgrid(0, 30*20, 0.5)
  
  
  ## Get a prediction
  ## Out B1: oral exposure to 1 mg/kg-d, matrix: Plasma, urine and feces
  out.B1 <-model %>%
    param(pars) %>%
    Req(blood, kidney, liver) %>%
    update(atol = 1E-50, maxsteps = 50000) %>%
    mrgsim_d(data=ex.oral.B1, tgrid=tsamp) %>%
    filter(time!=0) 
  out.B1 <- cbind.data.frame(time=out.B1$time, 
                             blood=out.B1$blood,
                             Kidney=out.B1$kidney,
                             Liver=out.B1$liver)  

  return(list("out.B1"=out.B1))
  
}

pred.ratB2 <- function(pars) {
  
  ## Get out of log domain
  # pars %<>% lapply(exp)                 ## return a list of exp (parameters) from log domain
  pars %<>% lapply(exp) 
  
  ## Define the three exposure scenario
  ## 0.1mg/kg/d Exposure scenario for oral 
  
  BW          = 0.29                                       # Rat body weight
  tinterval   = 24                                        # Time interval
  TDoses      = 1                                         # Total dosing/Dose times
  Ka = 32.894
  KLfabp = 81.260
  
  # To create a data set of 1 subject receiving DOSEoral.A every 24 hours for 1 total doses
  ex.oral.B2   <- ev (ID = 1, amt= 290000*0.1, 
                      ii = tinterval, addl=TDoses-1, cmt="yGL", replicate = FALSE)  
  
  ## set up the exposure time
  tsamp=tgrid(0, 30*20, 0.5)
  
  
  ## Get a prediction
  ## Out B2: oral exposure to 1 mg/kg-d, matrix: Plasma, urine and feces
  out.B2 <-model %>%
    param(pars) %>%
    Req(blood, kidney, liver) %>%
    update(atol = 1E-50, maxsteps = 50000) %>%
    mrgsim_d(data=ex.oral.B2, tgrid=tsamp) %>%
    filter(time!=0) 
  out.B2 <- cbind.data.frame(time=out.B2$time, 
                             blood=out.B2$blood,
                             Kidney=out.B2$kidney,
                             Liver=out.B2$liver)  
  
  return(list("out.B2"=out.B2))
  
}


## Cost fuction (FME) 
## Estimate the model residual by modCost function
MCcostA1<-function (pars){
  out<-  pred.ratA1 (pars)
  cost<- modCost  (model=out$out.A1, obs= Obs.A1, x="time")
    return(cost)
}

MCcostB1<-function (pars){
  out<-  pred.ratB1 (pars)
  cost<- modCost  (model=out$out.B1, obs= Obs.B1, x="time")
  return(cost)
}

MCcostB2<-function (pars){
  out<-  pred.ratB2 (pars)
  cost<- modCost  (model=out$out.B2, obs= Obs.B2, x="time")
  return(cost)
}

FitA1<- modFit(f=MCcostA1, p=theta, method ="Marq")
FitB1<- modFit(f=MCcostB1, p=theta, method ="Marq")
FitB2<- modFit(f=MCcostB2, p=theta, method ="Marq")

Sim.fit.A1 = pred.ratA1 (FitA1$par)$out.A1         ## Time-course PFOA concentration profiles using estimated parameters under exposure senaior A
df.sim.A1 = cbind.data.frame (Time=Sim.fit.A1$time, CA=Sim.fit.A1$blood, Kidney=Sim.fit.A1$Kidney, Liver=Sim.fit.A1$Liver)
Sim.fit.B1 = pred.ratB1 (FitB1$par)$out.B1         
df.sim.B1 = cbind.data.frame (Time=Sim.fit.B1$time, CA=Sim.fit.B1$blood, Kidney=Sim.fit.B1$Kidney, Liver=Sim.fit.B1$Liver)
Sim.fit.B2 = pred.ratB2 (FitB2$par)$out.B2         
df.sim.B2 = cbind.data.frame (Time=Sim.fit.B2$time, CA=Sim.fit.B2$blood, Kidney=Sim.fit.B2$Kidney, Liver=Sim.fit.B2$Liver)



################################################################################################################
###################################### ADDITIONAL CODE #############################################################
#################################################################################################################
# predict on PFOA concentrations in male rats exposed for 28 days at varying doses (t = 24 hr after last dose)
predict_serum_at_696h <- function(dose_mg_per_kg_day, pars, BW = 0.29) {
  # Dose in total mg
  total_dose_mg = dose_mg_per_kg_day * BW * 1000  # mg/kg * kg * 1000 = mg
  total_dose_ug = total_dose_mg * 1000            # convert to µg (model uses µg)
  
  # Dosing event
  dosing <- ev(ID = 1, amt = total_dose_ug, ii = 24, addl = 27, cmt = "yGL", replicate = FALSE)
  
  # Time grid up to 720 hrs to safely catch 696 hr
  tsamp <- tgrid(0, 720, 0.5)
  
  # Simulate
  out <- model %>%
    param(lapply(pars, exp)) %>%
    Req(blood) %>%
    mrgsim_d(data = dosing, tgrid = tsamp) %>%
    as.data.frame()
  
  # Get blood concentration at 696 hr
  blood_ng_per_ml <- out[out$time == 696, "blood"]
  blood_mg_per_L <- blood_ng_per_ml / 1e3  # convert ng/mL to mg/L
  
  return(blood_mg_per_L)
}

# Fitted parameters (choose one of your fits)
pars <- FitB1$par  # or FitA1$par / FitB2$par if more appropriate

# Doses in mg/kg-day
doses <- c(0.625, 1.25, 2.5, 5, 10)

# Predict for each dose
results <- sapply(doses, predict_serum_at_696h, pars = pars)

# Format and show results
pred_df <- data.frame(
  Dose_mg_kg_day = doses,
  Serum_mg_L_at_696hr = round(results, 4)
)
print(pred_df)
