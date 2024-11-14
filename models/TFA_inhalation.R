library(deSolve)

#### model code in
#### https://www.nist.gov/document/vinegar-toxicological-assessment-human-health-consequences-associated-inhalation-halon


# Define parameters
params <- c(
  BW = 70,     # Body weight (kg)
  H = 1.8,     # Height (m)
  A = 21,      # Age (yr)
  QPC = 17.4,  # Cardiac output per kg
  QLC = 0.0885,  # Liver blood flow proportion
  QGC = 0.2192,  # Gut blood flow proportion
  QFC = 0.0288,  # Fat blood flow proportion
  QSC = 0.2019,  # Slowly perfused tissues blood flow proportion
  QRC = 0.4616,  # Rapidly perfused tissues blood flow proportion
  VMAXC = 0.0,   # Michaelis-Menten Vmax (mg/hr/kg BW)
  KM = 1,        # Michaelis-Menten Km (mg/L)
  KFC = 0.0,     # First order metabolism rate constant (1/hr/kg BW)
  PLA = 0.038,   # Liver/air partition coefficient
  PGA = 0.056,   # Gut/air partition coefficient
  PFA = 0.745,   # Fat/air partition coefficient
  PSA = 0.074,   # Slowly perfused tissues/air partition coefficient
  PRA = 0.038,   # Richly perfused tissues/air partition coefficient
  PB = 0.033,    # Blood/air partition coefficient
  MW = 168.02,   # Molecular weight (g/mol)
  CI = 105000 * 168.02 / 24450, # Inhaled concentration (mg/L)
  F = 900,       # Breathing frequency /hr
  VGD = 0.2,     # Dead space gas volume
  VGP = 3.0,     # Pulmonary region gas volume
  VPR = 1.0,     # Ventilation perfusion ratio
  VCAPP = 0.169, # Volume of pulmonary capillaries
  VTP = 0.270,   # Volume of pulmonary tissues
  PI = 3.1415927 # Pi constant
)

# Initial state
init_state <- c(
  AS = 0,  # Amount in slowly perfused tissues (mg)
  AR = 0,  # Amount in rapidly perfused tissues (mg)
  AF = 0,  # Amount in fat (mg)
  AG = 0,  # Amount in gut (mg)
  AL = 0,  # Amount in liver (mg)
  CVGDD = 0,  # Concentration in distal deadspace (mg/L)
  CVGDP = 0,  # Concentration in proximal deadspace (mg/L)
  AI = 10,   # Initial amount inhaled (mg)
  AX = 0,   # Amount exhaled (mg)
  AP = 0,   # Amount in pulmonary region (mg)
  AM = 0    # Amount metabolized in liver (mg)
)

# Define the model
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Physiological parameters
    QP <- QPC * BW^0.74  # Cardiac output (L/hr)
    QPM <- QP / 60       # Cardiac output (L/min)
    QS <- QSC * QP       # Blood flow to slowly perfused tissues
    QR <- QRC * QP       # Blood flow to rapidly perfused tissues
    QF <- QFC * QP       # Blood flow to fat
    QL <- QLC * QP       # Blood flow to liver
    QG <- QGC * QP       # Blood flow to gut
    
    # Tissue volumes
    VL <- 0.027 * BW  # Liver volume
    VG <- 0.022 * BW  # Gut volume
    VF <- 0.215 * BW  # Fat volume
    VR <- 0.041 * BW  # Rapidly perfused volume
    VS <- 0.575 * BW  # Slowly perfused volume
    
    # Partition coefficients
    PL <- PLA / PB
    PG <- PGA / PB
    PF <- PFA / PB
    PS <- PSA / PB
    PR <- PRA / PB
    
    # Inhalation parameters
    VENTLM <- QPM / 0.7  # L/min
    VT <- (VENTLM * 60) / F  # Tidal volume (L)
    VTA <- VT - VGD  # Max volume of fresh air into pulmonary region
    VDOTA <- F * VTA  # Alveolar ventilation (L/hr)
    VDOTE <- VT * F   # Total ventilation (L/hr)
    
    # Gas exchange
    QCI <- VDOTA / VPR
    QL <- QLC * QCI
    QG <- QGC * QCI
    QF <- QFC * QCI
    QS <- QSC * QCI
    QR <- QRC * QCI
    QC <- QL + QG + QF + QS + QR  # Total cardiac output
    
    # Concentrations
    CA <- AP / (VGP / PB + VTP * PR + VCAPP)  # Arterial concentration
    CS <- AS / VS  # Concentration in slowly perfused tissues
    CR <- AR / VR  # Concentration in rapidly perfused tissues
    CF <- AF / VF  # Concentration in fat
    CL <- AL / VL  # Concentration in liver
    CG <- AG / VG  # Concentration in gut
    
    # Metabolism in liver using Michaelis-Menten kinetics
    RAM <- (VMAXC * CL) / (KM + CL)  # Rate of metabolism in liver
    RAL <- QL * (CA - CL) - RAM  # Rate of change in liver amount
    RAG <- QG * (CA - CG)        # Rate of change in gut amount
    RAF <- QF * (CA - CF)        # Rate of change in fat amount
    RAR <- QR * (CA - CR)        # Rate of change in rapidly perfused tissues
    RAS <- QS * (CA - CS)        # Rate of change in slowly perfused tissues
    
    # Inhalation adds to pulmonary region
    RAI <- VDOTE * CI           # Rate of inhalation
    RAP <- VDOTA * CI           # Transfer of inhaled to pulmonary region
    RAX <- VDOTA * CVGDD        # Rate of exhalation
    
    # Rate of change in pulmonary region amount
    RAP <- VDOTA * (CI - CVGDD)
    
    # Return the rate of change
    list(c(
      dAS = RAS,
      dAR = RAR,
      dAF = RAF,
      dAG = RAG,
      dAL = RAL,
      dCVGDD = VDOTE * (CA - CVGDD),  # Dead space dynamics
      dCVGDP = VDOTE * (CA - CVGDP),
      dAI = RAI,
      dAX = RAX,
      dAP = RAP,
      dAM = RAM
    ))
  })
}

# Time sequence
times <- seq(0, 300, by = 1)

# Solve the ODEs
out <- ode(y = init_state, times = times, func = model, parms = params)

# Output results
out_df <- as.data.frame(out)

# Display results
print(out_df)

# You can also visualize the results
library(ggplot2)
ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = AS, color = "Slowly Perfused")) +
  geom_line(aes(y = AR, color = "Rapidly Perfused")) +
  geom_line(aes(y = AF, color = "Fat")) +
  geom_line(aes(y = AG, color = "Gut")) +
  geom_line(aes(y = AL, color = "Liver")) +
  labs(y = "Amount (mg)", x = "Time (s)", color = "Compartment") +
  theme_minimal()

