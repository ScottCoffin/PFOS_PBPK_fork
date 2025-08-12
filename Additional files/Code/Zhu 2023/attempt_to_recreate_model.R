PFOA_mouse_pbpk_code <- '
$PARAM @annotated
BW         : 33       : Body weight (g)
QCC        : 16.5     : Cardiac output (L/h/kg^0.75)
QLC        : 0.161    : Fractional liver blood flow
QKC        : 0.091    : Fractional kidney blood flow
QLungC     : 0.005    : Fractional lung blood flow
QSkinC     : 0.058    : Fractional skin blood flow
QFilC      : 0.045    : Fractional kidney to filtrate
Htc        : 0.48     : Hematocrit

VLC        : 0.0589   : Liver volume fraction
VKC        : 0.0169   : Kidney volume fraction
VPlasC     : 0.049    : Plasma volume fraction
VFilC      : 0.0017   : Filtrate volume fraction
VLungC     : 0.00769  : Lung volume fraction
VSkinC     : 0.1653   : Skin volume fraction
VSTC       : 0.0075   : Stomach volume fraction
VSIC       : 0.038    : Small intestine volume fraction

PL         : 1.0286   : Liver/plasma partition coefficient
PK         : 0.0409   : Kidney/plasma PC
PLung      : 0.1760   : Lung/plasma PC
PDer       : 0.1953   : Skin/plasma PC
PRest      : 0.2531   : Rest/plasma PC

Free       : 0.02104  : Free fraction in plasma
Tmc        : 22427.29 : Tubular max capacity
Kt         : 73440    : Transporter affinity

K0C        : 0.0919   : Stomach to liver absorption rate
Kabsc      : 21.9535  : Small intestine to liver absorption rate
KunabsC    : 1.3609e-11 : Unabsorbed elimination

KbileC     : 0.0322   : Biliary elimination
KurineC    : 5.1118   : Urinary elimination

PDerCha    : 5.2693   : Dermal chamber-to-skin partition
Kvol       : 6.9934e-07 : Dermal volume transfer rate
Kp         : 0.001472 : Skin permeability
SA         : 0.164    : Skin surface area (m²)

GEC        : 0.188    : Gastric emptying constant

$CMT 
APlas_free ALiver AKidney ALung ASkin ARest AST ASI 
Afiltrate APTC 
AUCCA_free AUCLiver AUCKidney AUCLung AUCSkin AUCRest 
Afeces Aurine AabsST AabsSI ALung_dose ACham

$ODE
// Volumes (L)
double BW_kg = BW / 1000;
double VPlas = VPlasC * BW_kg;
double VL = VLC * BW_kg;
double VK = VKC * BW_kg;
double VLung = VLungC * BW_kg;
double VSkin = VSkinC * BW_kg;
double VRest = 0.93*BW_kg - VPlas - VL - VK - VLung - VSkin;

// Flows (L/h)
double QC = QCC * pow(BW_kg, 0.75) * (1 - Htc);
double QL = QLC * QC;
double QK = QKC * QC;
double QLu = QLungC * QC;
double QSk = QSkinC * QC;
double QRest = QC - QL - QK - QLu - QSk;

double GE = GEC * pow(BW_kg, -0.25);  // gastric emptying

// Tissue concentrations (mg/L)
double CPlas = APlas_free / VPlas / Free;
double CL = ALiver / VL;
double CK = AKidney / VK;
double CLu = ALung / VLung;
double CSk = ASkin / VSkin;
double CR = ARest / VRest;

double CVL = CL / PL;
double CVK = CK / PK;
double CVLu = CLu / PLung;
double CVSk = CSk / PDer;
double CVR = CR / PRest;

// Plasma dynamics
dxdt_APlas_free = QL * (CVL - CPlas) + QK * (CVK - CPlas) +
                  QLu * (CVLu - CPlas) + QSk * (CVSk - CPlas) +
                  QRest * (CVR - CPlas) + PDerCha * ASkin + AabsST + AabsSI + ALung_dose;

// Liver
dxdt_ALiver = QL * (CPlas - CVL) - KbileC * pow(BW_kg, -0.25) * ALiver + Kabsc * ASI + K0C * AST;
dxdt_AUCLiver = CL;

// Kidney
dxdt_AKidney = QK * (CPlas - CVK) - KurineC * pow(BW_kg, -0.25) * Afiltrate;
dxdt_AUCKidney = CK;

// Lung
dxdt_ALung = QLu * (CPlas - CVLu);
dxdt_AUCLung = CLu;

// Skin
dxdt_ASkin = QSk * (CPlas - CVSk) + Kvol * ACham - Kp * SA * CSk;
dxdt_AUCSkin = CSk;

// Rest of body
dxdt_ARest = QRest * (CPlas - CVR);
dxdt_AUCRest = CR;

// GI tract - oral dosing
dxdt_AST = -GE * AST - K0C * AST;
dxdt_AabsST = K0C * AST;

dxdt_ASI = GE * AST - Kabsc * ASI - KunabsC * ASI;
dxdt_AabsSI = Kabsc * ASI;

dxdt_Afeces = KunabsC * ASI + KbileC * pow(BW_kg, -0.25) * ALiver;
dxdt_Aurine = KurineC * pow(BW_kg, -0.25) * Afiltrate;

// Dermal chamber
dxdt_ACham = -Kvol * ACham;

// Nasal
dxdt_ALung_dose = -QLu * ALung_dose;

// Placeholder equations (no dynamics)
dxdt_Afiltrate = 0;
dxdt_APTC = 0;

// AUC
dxdt_AUCCA_free = CPlas;

$TABLE
capture Plasma = CPlas;
capture Liver = CL;
capture Kidney = CK;
capture Lung = CLu;
capture Skin = CSk;

capture AUC_CA = AUCCA_free;
capture AUC_CL = AUCLiver;
capture AUC_CK = AUCKidney;
capture AUC_CLung = AUCLung;
capture AUC_CSkin = AUCSkin;
'

mod.pfoa <- mcode("PFOA_mouse_pbpk", PFOA_mouse_pbpk_code)

ev_oral <- ev(amt = 1000 * 0.033, cmt = "AST")  # Oral dose
out <- mod.pfoa %>% ev(ev_oral) %>% mrgsim()
plot(out)

## test to reproduce Excel Table S13-1 The predicted concentrations of PFAS in mice plasma and liver dosing at 0.023 mg/kg/d for 7 days.
library(mrgsolve)
library(dplyr)
library(ggplot2)

# Compile model
mod.pfoa <- mcode("PFOA_mouse_pbpk", PFOA_mouse_pbpk_code)

# Simulation parameters
BW <- 33                         # grams
dose_mg_per_kg <- 0.023          # mg/kg/day
dose_mg <- dose_mg_per_kg * (BW / 1000)

ev_oral <- ev(
  amt = dose_mg,
  ii = 24,
  addl = 6,
  cmt = "AST"
)

# Run simulation
out <- mod.pfoa %>%
  ev(ev_oral) %>%
  mrgsim(end = 24*1, delta = 1/6) %>%  # simulate 1 day with 4-h resolution
  as_tibble()

# Convert units
converted <- out %>%
  transmute(
    Time_day = time / 24,
    Plasma_ng_mL = Plasma * 1000,   # mg/L → ng/mL
    Liver_ng_g = Liver * 1000       # mg/L → ng/g (≈ density 1)
  )

print(converted, n = 25)

