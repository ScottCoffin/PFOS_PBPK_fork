library(mrgsolve)

# Mouse PBPK Model for PFOS and PFHxS in mrgsolve, based on Chou et al. (2019) with new publication modifications
mouse_PFOS_PFHxS_PBPK.code <- '
$PARAM @annotated
BW                  : 0.025   : kg, Bodyweight (mice)
QCC                 : 16.5    : L/h/kg^0.75, Cardiac output
QLC                 : 0.161   : Fraction, Fraction blood flow to liver
QKC                 : 0.091   : Fraction, Fraction blood flow to kidney
Htc                 : 0.48    : Hematocrit, Hematocrit for mice
VLC                 : 0.055   : Fraction, Fractional liver tissue
VKC                 : 0.017   : Fraction, Fractional kidney tissue
VPlasC              : 0.049   : L/kg BW, Fractional plasma
VfilC               : 0.0017  : L/kg BW, Fraction vol. of filtrate
Free                : 0.09    : Fraction, Free fraction (initial for PFOS)
PL                  : 3.720   : Unitless, Liver/plasma partition coefficient
PK                  : 0.80    : Unitless, Kidney/plasma partition coefficient
PRest               : 0.20    : Unitless, Rest of body/plasma partition coefficient
Kabsc               : 2.12    : 1/(h*BW^-0.25), Absorption rate from small intestine to liver
KunabsC             : 7.05e-5 : 1/(h*BW^-0.25), Rate of unabsorbed dose to feces
GEC                 : 0.540   : 1/(h*BW^0.25), Gastric emptying time
K0C                 : 1.000   : 1/(h*BW^-0.25), Rate of uptake from the stomach into the liver
KbileC              : 0.004   : 1/(h*BW^-0.25), Biliary elimination rate
KurineC             : 1.60    : 1/(h*BW^-0.25), Urinary elimination rate
Kp                  : 0.001   : cm/h, Permeability coefficient of the dermis
SA                  : 0.164   : cm^2, Exposed surface area for dermal exposure
PDerCha             : 5.27    : Unitless, Dermis-to-chamber partition coefficient
PDer                : 0.321   : Unitless, Dermis-to-plasma partition coefficient
Kvol                : 6.99e-7 : 1/h, Loss rate from the chamber

$MAIN
// Physiological parameters based on body weight
double QC = QCC*pow(BW, 0.75)*(1-Htc);          // L/h, Cardiac output (adjusted for plasma)
double QL = QLC*QC;                             // L/h, Plasma flow to liver
double QK = QKC*QC;                             // L/h, Plasma flow to kidney
double QRest = QC - QK - QL;                    // L/h, Plasma flow to rest of body

// Tissue volumes
double VL = VLC*BW;                             // L, Volume of liver
double VK = VKC*BW;                             // L, Volume of kidney
double VPlas = VPlasC*BW;                       // L, Volume of plasma
double Vfil = VfilC*BW;                         // L, Volume of filtrate

// GI Absorption rates
double Kabs = Kabsc*pow(BW, -0.25);             // 1/h, Absorption rate from small intestine to liver
double Kunabs = KunabsC*pow(BW, -0.25);         // 1/h, Rate of unabsorbed dose to feces
double GE = GEC*pow(BW, -0.25);                 // 1/h, Gastric emptying time
double K0 = K0C*pow(BW, -0.25);                 // 1/h, Rate of uptake from stomach to liver

$CMT AST ASI APlas_free AL AK Afeces Aurine ACha ADer

$ODE
// Concentrations in plasma, liver, kidney, and rest of body
double CA_free = APlas_free / VPlas; // Free PFOS or PFHxS concentration in plasma
double CA = CA_free / Free;          // Total PFOS or PFHxS concentration in plasma
double CL = AL / VL;                 // PFOS or PFHxS concentration in liver
double CK = AK / VK;                 // PFOS or PFHxS concentration in kidney

// Plasma compartment dynamics
dxdt_APlas_free = QRest * (CA - CL/PL) + QK * (CA - CK/PK) + QL * (CA - CL) * Free;

// Liver compartment dynamics with absorption from small intestine
dxdt_AL = QL * (CA - CL / PL) * Free - KbileC * AL + Kabs * ASI + K0 * AST;

// Kidney compartment with urinary elimination
dxdt_AK = QK * (CA - CK / PK) * Free - KurineC * AK;

// GI Tract compartments for oral exposure
dxdt_AST = -K0 * AST - GE * AST;         // Stomach compartment for oral dose
dxdt_ASI = GE * AST - Kabs * ASI - Kunabs * ASI;  // Small intestine compartment for oral dose

// Fecal and Urine elimination
dxdt_Afeces = KbileC * AL + Kunabs * ASI; // Fecal elimination route
dxdt_Aurine = KurineC * AK;               // Urinary elimination route

// Dermal absorption dynamics
double RCha = Kp * SA * (ADer / PDerCha - ACha) - Kvol * ACha;
dxdt_ACha = RCha; // Chamber compartment
dxdt_ADer = QL * (CA - CL) * Free - KbileC * ADer + RCha; // Dermis compartment

$TABLE
// Output for capturing plasma, liver, and kidney concentrations
capture Plasma = CA_free / Free;
capture Liver = CL;
capture Kidney = CK;
capture Feces = Afeces;
capture Urine = Aurine;
'

# Load model into mrgsolve
mod_PFOS_PFHxS <- mcode("mouse_PFOS_PFHxS_PBPK", mouse_PFOS_PFHxS_PBPK.code)

# Example: Single oral dose event for PFOS
oral_dose <- ev(amt = 0.1, cmt = "AST", ii = 24, addl = 0)  # 0.1 mg single dose
out_oral <- mod_PFOS_PFHxS %>% ev(oral_dose) %>% mrgsim(end = 168, delta = 1)
# Define a single dermal dose event of 0.1 mg
dermal_dose <- ev(amt = 0.1, cmt = "ACha", ii = 24, addl = 0)  # 0.1 mg dose
# Run the simulation with the dermal dose event
output_dermal <- mod_PFOS_PFHxS %>%
  ev(dermal_dose) %>%
  mrgsim(end = 168, delta = 1)  # Run for 168 hours with 1-hour intervals

# Plot the concentration results for a dermal dose
plot(output_dermal, xlab = "Time (hours)", ylab = "Concentration (mg/L)")


# Run the simulation with the oral dose event
output_oral <- mod_PFOS_PFHxS %>%
  ev(oral_dose) %>%
  mrgsim(end = 168, delta = 1)  # Run for 168 hours with 1-hour intervals

# Plot the concentration results for an oral dose
plot(output_oral, xlab = "Time (hours)", ylab = "Concentration (mg/L)")

# Convert output to a data frame
df_oral <- as.data.frame(output_oral)
df_dermal <- as.data.frame(output_dermal)

# Access plasma, liver, and kidney concentrations for oral dose
df_oral[, c("time", "Plasma", "Liver", "Kidney")]

# Access plasma, liver, and kidney concentrations for dermal dose
df_dermal[, c("time", "Plasma", "Liver", "Kidney")]

