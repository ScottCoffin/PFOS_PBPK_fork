PreG_RatPBPK.code <- '
$PROB
## Pre-pregnant PFOS PBPK model for female SD rat
- Author    : Wei-Chun Chou
- Date      : March, 2020
- Strucutre : GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC
- Note.1    : PL, PRest, Vmax_apical_invitro, Km_api, RAF_api, RAF_baso,  KeffluxC, Kdif, KbileC were optimized by Chou and Lin, 2019
- Note.2    : BW, Free, GFRC, PF, PM were based on female parameters value

$PARAM @annotated
BW                  : 0.185   : kg,                  Bodyweight (measrument data if available)     ; Value obtained from Garner et al., 2015
Htc                 : 0.46    : Unitless,            Hematocrit for Rat                            ; Value obtained from Davies and Morris, 1993; Brown et al., 1997
QCC                 : 14.000  : L/h/kg^0.75,         Cardiac output                                ; Value obtained from Brown et al., 1997
QLC                 : 0.183   : Unitless,            Fraction blood flow to liver                  ; Value obtained from Brown et al., 1997
QKC                 : 0.141	  : Unitless,            Fraction blood flow to kidney                 ; Value obtained from Brown et al., 1997
QMC                 : 0.002   : Unitless,            Fraction blood flow to Mammary gland          ; Value obtained from Hanwell and Linzell, 1973
QFC                 : 0.07    : Unitless,            Fraction blood flow to Fat                    ; Value obtained from Brown et al., 1997
VLC                 : 0.035   : Unitless,            Fraction of liver volume                      ; Value obtained from Brown et al., 1997
VKC                 : 0.0084  : Unitless,            Fraction of kidney volume                     ; Value obtained from Brown et al., 1997
VMC                 : 0.01    : Unitless,            Fraction of Mammary gland tissue              ; Value obtained from Hanwell and Linzell, 1973
VFC                 : 0.07    : Unitless,            Fraction of fat tissue                        ; Assume 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0466  : L/kg BW,             Fraction of plasma volume                     ; Value obtained from Davies 1993
VFilC               : 8.4e-4  : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obtained from Worley et al., 2017
VPTCC               : 1.35e-4 : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6  : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtained from Addis et al., 1936
PL                  : 3.66    : Unitless,            Liver-to-plasma partition coefficient         ; Value was optimized by Chou and Lin, 2019 
PK                  : 0.80    : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from Chou and Lin, 2019
PM                  : 0.16    : Unitless,            Mammary gland-to-plasma partition coefficient ; Value obtained from Loccisano et al., 2012
PF                  : 0.13    : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from Loccisano et al., 2012
PRest               : 0.26    : Unitless,            Rest of body partition coefficient            ; Value was optimized by Chou and Lin, 2019 
MW                  : 500.126 : g/mol,               PFOS molecular mass                           ; Value obtained from Worley and Fisher, 2015
Free                : 0.09    : Unitless,            Female free fraction                          ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA 
KbileC              : 0.0026  : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value was optimized by Chou and Lin, 2019 
KurineC             : 1.60    : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
K0C                 : 1.000   : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KabsC               : 2.12    : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KunabsC             : 7.05e-5 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
GFRC                : 41.04   : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value obtained from Corley, 2005
GEC                 : 1.4     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value obtained from Yang et al., 2013
Vmax_baso_invitro   : 393.45  : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
Km_baso             : 27.2    : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al., 2007
Vmax_apical_invitro : 1808    : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimzed from Chou and Lin, 2019 
Km_apical           : 278     : mg/L,                Km of apical transporter                       ; Value optimzed from Chou and Lin, 2019
RAFbaso             : 1.90    : Unitless             Relative activity factor                      ; Value optimzed from Chou and Lin, 2019 
RAFapi              : 4.15    : Unitless             Relative acitivty factor                      ; Value optimzed from Chou and Lin, 2019 
Kdif                : 5.1e-4  : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value optimzed from Chou and Lin, 2019 
KeffluxC            : 2.09    : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimzed from Chou and Lin, 2019  

$MAIN
// #+ Blood flows
// #+ QC     : L/h, Cardiac output (adjusted for plasma)
// #+ QK     : L/h, Plasma flow to kidney
// #+ QL     : L/h, Plasma flow to liver
// #+ QM     : L/h, Plasma flow to mammary gland
// #+ QF     : L/h, Plasma flow to Fat
// #+ QRest  : L/h, Plasma flow to the rest of body

double QC    = QCC*pow(BW, 0.75)*(1-Htc);                  
double QK    = QKC*QC;                                    
double QL    = QLC*QC;                                     
double QM    = QMC*QC;                                     
double QF    = QFC*QC;                                     
double QRest = QC - QK - QL - QM - QF;                          

// #+ Tissue volumes
// #+ VL     : L, Volume of liver
// #+ VK     : L, Volume of kidney
// #+ VM     : L, Volume of mammary gland
// #+ VF     : L, Volume of fat
// #+ VPlas  : L, Volume of plasma
// #+ VFil   : L, Volume of filtrate
// #+ VPTC   : L, Volume of proximal tubule cells
// #+ VKb    : L, Volume of blood in the kidney; fraction blood volume of kidney (0.16) from Brown, 1997
// #+ VRest  : L, Rest of body; volume of remaining tissue (L); 
// #+ MK     : g, Kidney weight in gram

double VL    = VLC*BW;                                     
double VK    = VKC*BW;                                     
double VM    = VMC*BW;                                     
double VF    = VFC*BW;                                     
double VPlas = VPlasC*BW;                               
double VFil  = VFilC*BW;                                  
double VPTC  = VK*VPTCC;                                  
double VKb   = VK*0.16;                                   
double VRest = (0.92*BW) - VL - VK - VM - VF - VPlas;   

// #+ Kidney realted parameters
// #+ MK     : g, kidney weight in gram
// #+ PTC    : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney); Revised original equation (PTC = VKC*6e7) from Worley and Fisher, 2015
// #+ MPTC   : g, mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC: mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3) equation from Worley and Fisher, 2015
// #+ Vmax_apicalC: mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1) equation from Worley and Fisher, 2015
// #+ Kbile  : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine : 1/h, Urinary elimination
// #+ Kefflux: 1/h, Efflux clearance rate from PTC to blood
// #+ GFR    : L/h, Glomerular filtration rate, scaled to mass of kidney 

double MK           = VK*1000;                                    
double PTC          = VKC*6e7*1000;                              
double MPTC         = VPTC*1000;                                
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);         
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);      
double Vmax_baso    = Vmax_basoC*pow(BW,0.75);             
double Vmax_apical  = Vmax_apicalC*pow(BW,0.75);         
double Kbile        = KbileC*pow(BW,(-0.25));                  
double Kurine       = KurineC*pow(BW,(-0.25));                 
double Kefflux      = KeffluxC*pow(BW,(-0.25));               
double GFR          = GFRC*(MK/1000);                            

// #+ GI tract parameters
// #+ Kabs          : 1/h, Rate constant of absorption of PFOS from small intestine to liver
// #+ Kunabs        : 1/h, Rate constant of unabsorbed dose to appear in feces
// #+ GE            : 1/h, Gastric emptying rate
// #+ K0            : 1/h, Rate constant of uptake from the stomach into the liver

double Kabs         = KabsC*pow(BW,(-0.25));                    
double Kunabs       = KunabsC*pow(BW,(-0.25));                
double GE           = GEC*pow(BW,(-0.25));                          
double K0           = K0C*pow(BW,(-0.25));  

$CMT A_baso A_apical Adif Aefflux ACI APlas_free APTC AFil Aurine AKb ARest AST
AabsST ASI AabsSI Afeces AL AF AM // ADOSE 

$ODE
// #+ Concentrations in the tissues and in the venous plasma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for dam; CX indicate the concentration in the tissue (e.g., CL); CVX represent the concentration of chemical in tissue leaving tissue 
// #+ CPlas_free   : mg/L, Free PFOS concentration in the plasma
// #+ CPlas        : mg/L, Concentration of total PFOS in the plasma
// #+ CL           : mg/L, Concentration of PFOS in the liver compartment
// #+ CKb          : mg/L, Concetraitons of PFOS in venous plasma leaving kidney
// #+ CK           : mg/L, Concetraitons of PFOS in Kidney compartment
// #+ CM           : mg/L, Concentration of PFOS in the mammary gland compartment
// #+ CF           : mg/L, Concentration of PFOS in the fat compartment
// #+ CRest        : mg/L, Concentration of PFOS in the rest of the body
// #+ CPTC         : mg/L, Concentration of PFOS in proximal tubule cells (PTC)
// #+ CFil         : mg/L, Concentration of PFOS in filtrate (fil)
// #+ CVL          : mg/L, Concentration of PFOS in plasma leaving liver
// #+ CVK          : mg/L, Concentration of PFOS in plasma leaving kidney
// #+ CVM          : mg/L, Concentration of PFOS in plasma leaving mammary gland
// #+ CVF          : mg/L, Concentration of PFOS in plasma leaving fat
// #+ CVRest       : mg/L, Concentration of PFOS in plasma leaving rest of body


double CPlas_free  = APlas_free/VPlas;                      
double CPlas       = CPlas_free/Free;                               
double CL          = AL/VL;                                      
double CKb         = AKb/VKb;                                   
double CK          = CKb*PK;                                     
double CM          = AM/VM;                                      
double CF          = AF/VF;                                       
double CRest       = ARest/VRest;                             
double CPTC        = APTC/VPTC;                                
double CFil        = AFil/VFil;                                
double CVL         = CL/PL;                                     
double CVK         = CKb;                                       
double CVM         = CM/PM;                                     
double CVF         = CF/PF;                                     
double CVRest      = ARest/(VRest*PRest);                    

// #+ Equation for estimation of the rate of each compartment in the PBPK model for female rat
// #+ RA_baso      : mg/h, Rate of basolateral transporters
// #+ RA_apical    : mg/h, Rate of apical transporter
// #+ Rdif         : mg/h, Rate of diffusion from blood into the PTC
// #+ RAefflux     : mg/h, Rate of efflux clearance from PTC to blood
// #+ RCI          : mg/h, Rate of clerance (CL) to urine via glomerular filtration (GFR)
// #+ RPTC         : mg/h, Rate of change in PTC
// #+ RFil         : mg/h, Rate of change in Fil
// #+ RKb          : mg/h, Rate of change in Kidney-serum compartment
// #+ RST          : mg/h, Rate of change in stomach compartment
// #+ RSI          : mg/h, Rate of change in small intestines
// #+ RabsST       : mg/h, Rate of absorption in Stomach
// #+ RabsSI       : mg/h, Rate of absorption in small intestines
// #+ RL           : mg/h, Rate of change in liver compartment
// #+ RF           : mg/h, Rate of chnage in fat comaprtment
// #+ RM           : mg/h, Rate of change in mammary tissues
// #+ RRest        : mg/h, Rate of change in rest of body
// #+ RPlas_free   : mg/h, Rate of free PFOS change in the dam plasma
// #+ Rurine       : mg/h, Rate of change in urine
// #+ Rbile        : mg/h, Rate of change in bile compartment 

double RA_baso     = (Vmax_baso*CKb)/(Km_baso+CKb);                                             
double RA_apical   = (Vmax_apical*CFil)/(Km_apical + CFil);                                   
double Rdif        = Kdif*(CKb - CPTC);                                                            
double RAefflux    = Kefflux*APTC;                                                             
double RCI         = CPlas*GFR*Free;                                                                   
double RPTC        = Rdif + RA_apical + RA_baso - RAefflux;                                        
double RFil        = CPlas*GFR*Free - RA_apical - AFil*Kurine;                                        
double RKb         = QK*(CPlas-CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;                               
double RST         = - K0*AST - GE*AST;                                                             
double RSI         = GE*AST - Kabs*ASI - Kunabs*ASI;                                                
double RabsST      = K0*AST;                                                                     
double RabsSI      = Kabs*ASI;                                                                   
double RL          = QL*(CPlas-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;                                
double RF          = QF*(CPlas-CVF)*Free;                                                               
double RM          = QM*(CPlas-CVM)*Free;                                                               
double RRest       = QRest*(CPlas-CVRest)*Free;                                                      
double RPlas_free  = (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) + (QM*CVM*Free) + (QF*CVF*Free) - (QC*CPlas*Free) + RAefflux;  
double Rurine      = Kurine*AFil;                                                                
double Rfeces      = Kbile*AL + Kunabs*ASI;                                                      

// #+ ODE equations for compartments in the female rat
dxdt_A_baso        = RA_baso;                                                                      
dxdt_A_apical      = RA_apical;                                                                  
dxdt_Adif          = Rdif;                                                                           
dxdt_Aefflux       = RAefflux;                                                                    
dxdt_ACI           = RCI;                                                                             
dxdt_APTC          = RPTC;                                                                           
dxdt_AFil          = RFil;                                                                          
dxdt_AKb           = RKb;                                                                             
dxdt_AST           = RST; 
dxdt_ASI           = RSI;                                                                             
dxdt_AabsST        = RabsST;
dxdt_AabsSI        = RabsSI;
dxdt_AL            = RL;
dxdt_AF            = RF;
dxdt_AM            = RM;
dxdt_ARest         = RRest;                                                                         
dxdt_APlas_free    = RPlas_free;                                                               
dxdt_Aurine        = Rurine;
dxdt_Afeces        = Rfeces;

// #+ Virtural compartment; input dose
// dxdt_ADOSE         = 0;

// #+ Mass Balance check 
// double ATissue     = APlas_free + ARest + AKb + AFil + APTC + AL + AM + AF + AST + ASI;
// double ALoss       = Aurine + Afeces;
// double ATotal      = ATissue + ALoss;
// double MBal        = ADOSE - ATotal; 

$TABLE
capture Plasma     = CPlas_free/Free;
capture Liver      = AL/VL;
capture Kidney     = CK;
'

GRatPBPK.code <- '
$PROB
# Gestational PFOS PBPK model for SD rat
- Author    : Wei-Chun Chou
- Date      : March, 2020
- Structure : GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, rest of Fetuss, fetal liver, amniotic fluid  
- Note.1    : Initial physiological parameters and optimized parameters ares matched with the values in non-pregnnat PBPK models
- Note.2    : Growth equation of physioloical parameters values were taken from Loccisano et al. (2012): https://doi.org/10.1016/j.reprotox.2011.07.003

$PARAM @annotated
// #+ Dams parameters
BW                  : 0.185   : kg,                  Bodyweight (pre-pregnant body weight)         ; Value obtained from Garner et al., 2015
Htc                 : 0.46    : Unitless,            Hematocrit for Rat                            ; Value obtained from Davies and Morris, 1993; Brown et al., 1997
QLC                 : 0.183   : Unitless,            Fraction blood flow to liver                  ; Value obtained from Brown et al., 1997
QKC                 : 0.141	  : Unitless,            Fraction blood flow to kidney                 ; Value obtained from Brown et al., 1997
QMC                 : 0.002   : Unitless,            Fraction blood flow to Mammary gland          ; Value obtained from Hanwell and Linzell, 1973
QFC                 : 0.07    : Unitless,            Fraction blood flow to Fat                    ; Value obtained from Brown et al., 1997
VLC                 : 0.035   : Unitless,            Fraction of liver volume                      ; Value obtained from Brown et al., 1997
VKC                 : 0.0084  : Unitless,            Fraction of kidney volume                     ; Value obtained from Brown et al., 1997
VMC                 : 0.01    : Unitless,            Fraction of Mammary gland tissue              ; Value obtained from Hanwell and Linzell, 1973
VFC                 : 0.07    : Unitless,            Fraction of fat tissue                        ; Assume 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0466  : L/kg BW,             Fraction of plasma volume                     ; Value obtained from Davies 1993
VFilC               : 8.4e-4  : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obtained from Worley et al., 2017
VPTCC               : 1.35e-4 : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6  : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtained from Addis et al., 1936 
PL                  : 3.66    : Unitless,            Liver-to-plasma partition coefficient         ; Value was optimized by Chou and Lin, 2019 
PK                  : 0.80    : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from Loccisano et al., 2012
PM                  : 0.16    : Unitless,            Mammary gland-to-plasma partition coefficient ; Value obtained from Loccisano et al., 2012
PF                  : 0.13    : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from Loccisano et al., 2012
PRest               : 0.22    : Unitless,            Rest of body partition coefficient            ; Value obtained from Loccisano et al., 2012 
PPla                : 0.41    : Unitless,            placenta-to-plasmapartition coefficient       ; Value obtained from Loccisano et al., 2012
MW                  : 500.126 : g/mol,               PFOS molecular mass                           ; Value obtained from Worley and Fisher, 2015
Free                : 0.022   : Unitless,            Free fraction                                 ; value obtained from Loccisano et al., (2012)
KbileC              : 0.0026  : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value was optimized by Chou and Lin, 2019 
KurineC             : 1.60    : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
K0C                 : 1.000   : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KabsC               : 2.12    : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KunabsC             : 7.05e-5 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
GFRC                : 41.04   : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value obtained from Corley, 2005
GEC                 : 1.4     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value obtained from Yang et al., 2013
Vmax_baso_invitro   : 393.45  : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
Km_baso             : 27.2    : mg/L,                Km of basolateral transporter                 ; Value calculated from Worley et al., 2007
Vmax_apical_invitro : 1808    : pmol/mg protein/min, Vmax of apical transporter                    ; Value optimzed from Chou and Lin, 2019 
Km_apical           : 278     : mg/L,                Km of apical transporter                      ; Value optimzed from Chou and Lin, 2019
RAFbaso             : 1.90    : Unitless             Relative activity factor                      ; Value optimzed from Chou and Lin, 2019 
RAFapi              : 4.15    : Unitless             Relative acitivty factor                      ; Value optimzed from Chou and Lin, 2019 
Kdif                : 5.1e-4  : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value optimzed from Chou and Lin, 2019 
KeffluxC            : 2.09    : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimzed from Chou and Lin, 2019   
Ktrans1C            : 0.46    : L/h/kg^0.75,         Mother-to-fetus placental transfer rate       ; Value from Loccisano et al., 2013; assumed to be same as PFOA
Ktrans2C            : 1.00    : L/h/kg^0.75,         Fetus-to-mother placental transfer rate       ; Value from Loccisano et al., 2013; assumed to be same as PFOA
Ktrans3C            : 0.008   : L/h/kg^0.75,         Fetus-to-amniotic fluid transfer rate         ; Value from Loccisano et al., 2012; assumed to be same as PFOA 
Ktrans4C            : 0.001   : L/h/kg^0.75,         Amniotic fluid-to-fetus transfer rate         ; Value from Loccisano et al., 2012; assumed to be same as PFOA 

// #+ Fetal parameters
N                   : 8       : number,              Number of fetus                               ; Value from Loccisano et al. (2012)
VPlasC_Fet          : 0.047   : Unitless,            Fraction volume of fetal plasma               ; Value was assumed to be same as the mother
Free_Fet            : 0.022   : Unitless,            Free fraction of PFOS in plasma for Fetus     ; Value was assumed to be same as the mother
PL_Fet              : 3.66    : Unitless,            Liver-to-plasma partition coefficient in Fetus; Value was assumed to be same as the mother
PRest_Fet           : 0.22    : Unitless,            Rest of body-to-plasma Partition coefficient  ; Value was assumed to be same as the mother
QLC_Fet             : 0.061   : Unitless,            Blood flows to fetal liver                    ; Value from Itskovitz, 1987 and Loccisano et al. (2012) (Sheep); 

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ GD             : Day, Gestational days
// #+ GA             : Week, Gestational age

double GD            = TIME/ 24;          
double GA            = TIME/ 168;          

// #+ Growth equations for the dam during pregnacy; Equation from Loccisano et al. (2012)
// #+ Placenta volume (L); all equations from Loccisano et al. (2012); OFlaherty, et al. 1992 and Gentry, et al. 2002.

double VPla          = 1e-7;
if (GD <= 6)       {VPla = 1e-7;} 
else if (GD <= 10) {VPla = (N*(8*(GD - 6)))/(1.0e6);} 
else if (GD > 10)  {VPla = ((N*(32*(exp(-0.23*(GD-10))))) + ((40*(exp(0.28*(GD-10))-1))))/1.0e6;}

// #+ Blood flow to rat placenta (L/h); equations from Loccisano et al. (2012)
// #+ Qdec1          : L/h, Plasma flow to the placenta during GD6 - 10; yolk sac placenta predominates
// #+ Qdec2          : L/h, Plasma flow to the placenta during GD10 -12; yolk sac placenta disappears
// #+ Qcap           : L/h, Chorioallantoic placenta appears
// #+ QPla1          : L/h, 25% of incr. in maternal blood flow assoc. with
// #+ QPla           : L/h, Plasma flow to placenta

double Qdec1         = 0;
double Qdec2         = 0;
double Qcap          = 0;

if (GD <= 6)       {QPla1 = 1e-7;} 
else if (GD <= 10) {Qdec1 = 0.55*(GD - 6);}                                 
else if (GD <= 12) {Qdec2 = (2.2*exp(-0.23*(GD-10)));}                         
else if (GD > 12)  {Qcap  = pow((0.1207*(GD-12)), 4.36);}   

double Qdec          = Qdec1 + Qdec2;
double QPla1         = (N*(0.02*Qdec + Qcap))/24;   
double QPla          = QPla1*(1-Htc); 


// #+ Changes of tissue volumes during pregnancy
// #+ VM             : L, Volume of mammary gland at GD0
// #+ VM_P           : L, Volume of mammary gland during pregnancy; Equation from Yoon et al. (2009), Loccisano et al. (2012) and Lin et al. (2013)
// #+ VF             : L, Volume of fat at GD0
// #+ VF_P           : L, Volume of fat gland during pregnancy; Equation from Yoon et al. (2009), Loccisano et al. (2012) and Lin et al. (2013)
// #+ VL             : L, Volume of liver
// #+ VK             : L, Volume of kidney
// #+ VKb            : L, Volume of kidney serum
// #+ VPlas          : L, Volume of plasma
// #+ VFil           : L, Volume of flitrate

double VM            = VMC*BW;  
double VM_P          = GD > 3? VM*(1 + 0.2*GD) : VM; 
double VF            = VFC*BW; 
double VF_P          = VF*(1 + 0.0182*GD); 
double VL            = VLC*BW;
double VK            = VKC*BW;
double VPlas         = VPlasC*BW;

// #+ Volumes of kidney and related compartmnet
double MK            = VKC*BW*1000;                                                          
double ML            = VLC*BW*1000;                     
double VPTC          = MK*VPTCC;                       
double VKb           = VK*0.16;          
double VFil          = VFilC*BW; 

// #+ Growth equations of tissues for the fetus druing pregnancy
// #+ VFet_1,        : kg, Fetus volume for one fetus; Equlation from Yoon et al. (2009)
// #+ VFet,          : kg, Fetus volume for whole litter; Equlation from Loccisano et al. (2012)
// #+ VAmX,          : L,  Amniotic fluid volume for one fetus; Equlation from Loccisano et al. (2012)
// #+ VAm,           : L,  Amniotic fluid volume for whole litter; Equlation from Loccisano et al. (2012)
// #+ VPlas_Fet,     : L,  Plasma volume for fetus; Equlation from Loccisano et al. (2012)
// #+ VLC_Fet,       : %,  Fractional liver tissue for fetus; Equlation from Loccisano et al. (2012); data from Scheidereit, 1985
// #+ VRest_Fet,     : L,  Rest of body volume for fetus; Equlation from Loccisano et al. (2012) 
// #+ VBal_Fet,      : L,  Check the balance for fetus volume

double VFet_1        = (0.1089 + (16*exp(-exp(5.515-0.2565*GD))))/1000;
double VFet          = VFet_1*N;
double VAmX          = GD >= 10 ? ((-4E-6)*pow(GD, 3) + 0.0002*pow(GD, 2) - 0.0023*GD + 0.0099) : 1E-7;
double VAm           = VAmX*N;
double VPlas_Fet     = VPlasC_Fet*VFet;                      
double VLC_Fet       = GD < 17 ? 0 : (-0.0013)*pow(GD, 3) + 0.0731*pow(GD, 2) - 1.375*GD + 8.6997;
double VL_Fet        = VLC_Fet*VFet;
double VRest_Fet     = 0.93*VFet - VPlas_Fet - VL_Fet;                
double VBal_Fet      = 0.93*VFet - (VRest_Fet + VPlas_Fet + VL_Fet);    

// #+ Fetal blood flows for fetus druing pregnancy
// #+ QC_Fet,        : L/h, Plasma flow to the fetus
// #+ QL_Fet,        : L/h, Plasma flow to the fetal liver
// #+ QRest_Fet,     : L/h, Plasma flow to the fetal rest of body
// #+ QFetBal,       : L/h, Check the balance of fetal blood flow

double QC_Fet        = (QPla/(1 + 20000*exp(-0.55*GD)));     
double QL_Fet        = QC_Fet*QLC_Fet;
double QRest_Fet     = QC_Fet - QL_Fet; 
double QFetBal       = QC_Fet - (QRest_Fet + QL_Fet);

// #+ Changes of maternal body weight during pregnancy
// #+ BW_P           : Dam BW during pregnancy
// #+ BWinc          : BW increased during pregnancy
// #+ VRest          : Volume of rest of body in Dam
// #+ VBal           : Volume Balance check 

double BW_P          = BW + (VF_P - VF) + (VM_P - VM) + VPla + VFet + VAm; 
double BWinc         = (VF_P - VF) + (VM_P - VM) + VPla + VFet + VAm; 
double VRest         = (0.93*BW_P) - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 
double VBal          = (0.93*BW_P) - VRest -(VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 

// #+ Equations for dam physiological parameters; growth equation obtained from Loccisano et al. (2012)
// #+ QCl_P          : L/h/kg, Cardiac index during gestation (Dowell, 1997)
// #+ QC_P           : L/h, Cardiac output during gestation (adjusted for plasma)
// #+ QCl            : L/h/kg, Cardiac index at GD0 (Dowell, 1997)
// #+ QC             : L/h, Cardiac output at GD0 (adjusted for plasma)
// #+ QF             : L/h, Plasma flow to fat at GD0
// #+ QF_P           : L/h, Plasma flow to fat during pregnacy
// #+ QM             : L/h, Plasma flow to mammary gland at GD0
// #+ QM_P           : L/h, Plasma flow to mammary gland during pregnacy
// #+ QL             : L/h, Plasma flow to liver
// #+ QK             : L/h, Plasma flow to kidney
// #+ QK_P           : L/h, Plasma flow to kidney during pregnacy 
// #+ QRest          : L/h, Plasma flow to the rest of body
// #+ QBal           : L/h, Blood flow balance check

double QCl_P         = 24.56 - 0.1323*GD;                        
double QC_P          = QCl_P*BW_P*(1-Htc);                        
double QCl           = 24.56 - 0.1323*0;                            
double QC            = QCl*BW*(1-Htc);                               
double QF            = QFC*QC;                                      
double QF_P          = QF*(VF_P/VF);                              
double QM            = QMC*QC;                                      
double QM_P          = QM*(VM_P/VM);                              
double QL            = QLC*QC;                                      
double QK            = QKC*QC;
double QK_P          = (-0.0014)*pow(GD, 2) + 0.0308*GD + 0.449;
double QRest         = QC_P - (QK_P + QL + QM_P + QF_P + QPla);                
double QBal          = QC_P - (QK_P + QL + QRest + QPla + QM_P + QF_P);

// #+ Pregnant rat Kinetics parameters
// #+ GFR            : L/h, Glomerular Filtration Rate (GFR); GFR was assumed 50% of renal plasma flow
// #+ GE             : 1/h, Gastric emptying rate
// #+ K0             : 1/h, Rate of uptake from the stomach into the liver
// #+ Kbile          : 1/h, Biliary elimination rate constant, liver to feces storage
// #+ Kurine         : 1/h, Urinary elimination rate constant
// #+ Kabs           : 1/h, Rate constant of absorption of PFOS from small intestine to liver
// #+ Kunabs         : 1/h, Rate constant of unabsorbed dose to appear in feces
// #+ Kinetics parmaters for kidney and its subcompartment
// #+ PTC            : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC           : g, mass of the proximal tubule cells (assuming density 1 kg/L) 
// #+ Vmax_basoC     : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_baso      : 
// #+ Vmax_apicalC   : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_apical    : 
// #+ Kefflux        : 1/h, Efflux clearance rate from PTC to blood 
// #+ Ktrans_1       : L/h, Rate constant for placental transfer; dam to fetus
// #+ Ktrans_2       : L/h, Rate constant for placental transfer; fetus to dam
// #+ Ktrans_3       : L/h, Amniotic fluid transfer rate; fetus to fluid
// #+ Ktrans_4       : L/h, Amniotic fluid transfer rate; fluid to fetus

double GFR           = GFRC*(MK/1000);
double GE            = GEC*pow(BW_P,(-0.25));                
double K0            = K0C*pow(BW_P,(-0.25));               
double Kbile         = KbileC*pow(BW_P,(-0.25));                   
double Kurine        = KurineC*pow(BW_P,(-0.25));                 
double Kabs          = KabsC*pow(BW_P,(-0.25));                    
double Kunabs        = KunabsC*pow(BW_P,(-0.25));                
double PTC           = VKC*6e7*1000;                                       
double MPTC          = VPTC*1000;                                         
double Vmax_basoC    = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_baso     = Vmax_basoC*pow(BW_P,0.75);                          
double Vmax_apicalC  = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);   
double Vmax_apical   = Vmax_apicalC*pow(BW_P,0.75);                      
double Kefflux       = KeffluxC*pow(BW_P,(-0.25));                        
double Ktrans_1      = Ktrans1C*(pow(VFet_1,0.75)*N);                       
double Ktrans_2      = Ktrans2C*(pow(VFet_1,0.75)*N);                      
double Ktrans_3      = Ktrans3C*(pow(VFet_1,0.75)*N);                       
double Ktrans_4      = Ktrans4C*(pow(VFet_1,0.75)*N);                      


// #+ Mass balance adjusted factor; avoding the negative occur at time = 0
//double KDOSE        = (TIME==0)?0:1;     

$INIT @annotated
// #+ Set up the initial concentration; 
//ADOSE            : 0   : mg, Amount of input dose; assumed a virtual compartment for validating the model mass balance
APlas_free         : 0   : mg, Amount of free PFOS in the plasma compartment
APTC               : 0   : mg, Amount of PFOS in the proximal tubule cells subcompartment
AFil               : 0   : mg, Amount of PFOS in the filtrate subcompartment
Aurine             : 0   : mg, Amount of PFOS in the urine virtual compartment
AKb                : 0   : mg, Amount of PFOS in the kidney blood compartment
ARest              : 0   : mg, Amount of PFOS in the rest of body compartment
Afeces             : 0   : mg, Amount of PFOS in the feces virtual compartment 
AL                 : 0   : mg, Amount of PFOS in the liver virtual compartment
AM                 : 0   : mg, Amount of PFOS in the mammary gland
AF                 : 0   : mg, Amount of PFOS in the fat compartment
A_baso             : 0   : mg, Amount of PFOS in the baso subcompartment 
A_apical           : 0   : mg, Amount of PFOS in the apical subcompartment 
Adif               : 0   : mg, Amount of PFOS in the diffusion virtual compartment 
Aefflux            : 0   : mg, Amount of PFOS in the efflux virtual compartment; simulation of PFOS from PTC to blood
Atrans_1           : 0   : mg, Amount of PFOS in the placental transfer from dam to Fetus  
Atrans_2           : 0   : mg, Amount of PFOS in the placental transfer from Fetus to dam 
Atrans_3           : 0   : mg, Amount of PFOS in the Amniotic fluid transfer from Amniotic fluid to fetus  
Atrans_4           : 0   : mg, Amount of PFOS in the Amniotic fluid transfer from fetus to Amniotic fluid 
APla               : 0   : mg, Amount of PFOS in the placenta compartment
ASI                : 0   : mg, Amount of PFOS in the small intestine compartment
AST                : 0   : mg, Amount of PFOS in the stomach compartment
AabsST             : 0   : mg, Amount of absorbed PFOS in the stomach compartment
AabsSI             : 0   : mg, Amount of absorbed PFOS in the small intestine compartment
AAm                : 0   : mg, Amount of PFOS in the Amniotic fluid compartment
AL_Fet             : 0   : mg, Amount of PFOS in the Fetal liver compartment
ARest_Fet          : 0   : mg, Amount of PFOS in the Fetal rest of Fetus body compartment
APlas_Fet_free     : 0   : mg, Amount of PFOS in the Fetal plasma compartment
AUC_CPlas          : 0   : mg/L*hr, Area under curve of PFOS in maternal plasma
AUC_CPlas_Fet      : 0   : mg/L*hr, Area under curve of PFOS in fetal plasma
  
$ODE 
// #+ Concentrations in the tissues and in the venous plasma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for dam; CX indicate the concentration in the tissue (e.g., CL); CVX represent the concentration of chemical in tissue leaving tissue 
// #+ CPlas_free   : mg/L, Free PFOS concentration in the plasma
// #+ CPlas        : mg/L, Concentration of total PFOS in the plasma
// #+ CL           : mg/L, Concentration of PFOS in the liver compartment
// #+ CKb          : mg/L, Concetraitons of PFOS in venous plasma leaving kidney
// #+ CK           : mg/L, Concetraitons of PFOS in Kidney compartment
// #+ CM           : mg/L, Concentration of PFOS in the mammary gland compartment
// #+ CF           : mg/L, Concentration of PFOS in the fat compartment
// #+ CAm          : mg/L, Concentration of PFOS in the amniotic fluid compartment
// #+ CRest        : mg/L, Concentration of PFOS in the rest of the body
// #+ CPTC         : mg/L, Concetraitons of PFOS in proximal tubule cells (PTC)
// #+ CFil         : mg/L, Concetraitons of PFOS in filtrate (Fil)
// #+ CPla         : mg/L, Concetraitons of PFOS in placenta

double CPlas_free = APlas_free/ VPlas;    
double CPlas       = CPlas_free/Free;                                       
double CL          = AL/VL;                                             
double CKb         = AKb/VKb;                                            
double CK          = CVK*PK;                                            
double CM          = AM/VM;                                              
double CF          = AF/VF;                                               
double CAm         = AAm/(VAm + 1E-7);
double CRest       = ARest/VRest;                                     
double CPTC        = APTC/VPTC;                                       
double CFil        = AFil/VFil;                                       
double CPla        = APla/VPla;
double CVL         = CL/PL;                                            
double CVK         = CKb;                                              
double CVM         = CM/PM;                                            
double CVF         = CF/PF;                                            
double CVRest      = CRest/PRest;                                   
double CVPla       = CPla/PPla;

// #+ Concentrations for fetus; 
// #+ CPlas_Fet_free   : mg/L, Free PFOS concentration in the plasma
// #+ CRest_Fet        : mg/L, Free PFOS concentration in the rest of body
// #+ CL_Fet           : mg/L, Free PFOS concentration in the liver

double CPlas_Fet_free  = APlas_Fet_free /(VPlas_Fet + 1e-7); 
double CPlas_Fet       = CPlas_Fet_free/Free_Fet;
double CRest_Fet       = ARest_Fet/(VRest_Fet+ 1e-7);
double CVRest_Fet      = CRest_Fet/ PRest_Fet;
double CL_Fet          = AL_Fet/ (VL_Fet+ 1e-7);
double CVL_Fet         = CL_Fet/ PL_Fet;


// #+ Equation for estimation the rates of compartments in the PBPK model for the pregnant rat
// #+ RA_baso      : mg/h, Rate of basolateral transporters
// #+ RA_apical    : mg/h, Rate of apical transporter
// #+ Rdif         : mg/h, Rate of diffusion from blood into the PTC
// #+ RAefflux     : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI          : mg/h, Rate of clerance (CL) to via glomerular filtration (GFR)
// #+ RPTC         : mg/h, Rate of change in PTC
// #+ RFil         : mg/h, Rate of change in Fil
// #+ RKb          : mg/h, Rate of change in Kidney-serum compartment
// #+ RST          : mg/h, Rate of change in stomach compartment
// #+ RSI          : mg/h, Rate of change in small intestines
// #+ RabsST       : mg/h, Rate of absorption in Stomach
// #+ RabsSI       : mg/h, Rate of absorption in small intestines
// #+ RL           : mg/h, Rate of change in liver compartment
// #+ RF           : mg/h, Rate of chnage in fat comaprtment
// #+ RM           : mg/h, Rate of change in mammary tissues
// #+ RRest        : mg/h, Rate of change in rest of body
// #+ RPla         : mg/h, Rate of change in placenta compartment
// #+ RPlas_free   : mg/h, Rate of free PFOS change in the dam plasma
// #+ Rurine       : mg/h, Rate of change in urine
// #+ Rtrans_1     : mg/h, Rate of change in placenta tranfer from mother to fetus
// #+ Rtrans_2     : mg/h, Rate of change in placenta tranfer from fetus to mother
// #+ Rtrans_3     : mg/h, Rate of change in amniotic fluid tranfer from fetus to fluid
// #+ Rtrans_4     : mg/h, Rate of change in amniotic fluid tranfer from fluid to fetus
// #+ RAm          : mg/h, Rate of change in amniotic fluid compartment 
// #+ Rbile        : mg/h, Rate of change in bile compartment 
// #+ RL_Fet       : mg/h, Rate of change in fetal liver
// #+ RRest_Fet    : mg/h, Rate of change in fetal rest of body
// #+ RPlas_Fet    : mg/h, Rate of free PFOS change in the fetal plasma

double RA_baso        = (Vmax_baso*CKb)/(Km_baso + CKb);                 
double RA_apical      = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif           = Kdif*(CKb - CPTC);                               
double RAefflux       = Kefflux*APTC;                                
double RCI            = CPlas*GFR*Free;                             
double RPTC           = Rdif + RA_apical + RA_baso - RAefflux;           
double RFil           = RCI - RA_apical - AFil*Kurine;          
double RKb            = QK_P*(CPlas - CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RST            = -K0*AST - GE*AST;
double RabsST         = K0*AST;
double RSI            = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsSI         = Kabs*ASI;
double RL             = QL*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;   
double RF             = QF_P*(CPlas - CVF)*Free;            
double RM             = QM_P*(CPlas - CVM)*Free;                            
double RRest          = QRest*(CPlas - CVRest)*Free;                         
double Rurine         = Kurine*AFil;                                   
double Rfeces         = Kbile*AL + Kunabs*ASI;                         
double Rtrans_1       = Ktrans_1*CVPla*Free;                          
double Rtrans_2       = Ktrans_2*CPlas_Fet*Free;                           
double Rtrans_3       = Ktrans_3*CVRest_Fet*Free_Fet;                              
double Rtrans_4       = Ktrans_4*CAm;
double RAm            = Rtrans_3 - Rtrans_4;                               
double RPla           = QPla*(CPlas - CVPla)*Free + Rtrans_2 - Rtrans_1;          
double RPlas_free     = (QRest*CVRest*Free) + (QK_P*CVK*Free) + (QL*CVL*Free) + (QM_P*CVM*Free) + (QF_P*CVF*Free) + (QPla*CVPla*Free) - 
                        (QC_P*CPlas*Free) + RAefflux;  
double RL_Fet         = QL_Fet*(CPlas_Fet - CVL_Fet)*Free_Fet;     
double RRest_Fet      = QRest_Fet*(CPlas_Fet - CVRest_Fet)*Free_Fet - Rtrans_3 + Rtrans_4;     
double RPlas_Fet      = (QRest_Fet*CVRest_Fet*Free_Fet) + (QL_Fet*CVL_Fet*Free_Fet) - (QC_Fet*CPlas_Fet*Free_Fet) + Rtrans_1 - Rtrans_2;

// #+ ODE equation for compartments in the pregnant rat
dxdt_A_baso           = RA_baso;                                           
dxdt_A_apical         = RA_apical;                                       
dxdt_Adif             = Rdif;                                              
dxdt_Aefflux          = RAefflux;                                      
dxdt_APTC             = RPTC;                                            
dxdt_AFil             = RFil;                                            
dxdt_AKb              = RKb;                                                
dxdt_AST              = RST;
dxdt_ASI              = RSI;
dxdt_AabsST           = RabsST;
dxdt_AabsSI           = RabsSI;
dxdt_AL               = RL;                                                  
dxdt_AF               = RF;                                            
dxdt_AM               = RM;                                            
dxdt_ARest            = RRest;                                            
dxdt_Aurine           = Rurine;                                          
dxdt_Afeces           = Rfeces;                                          
dxdt_Atrans_1         = Rtrans_1;                                
dxdt_Atrans_2         = Rtrans_2;                                
dxdt_Atrans_3         = Rtrans_3;                                 
dxdt_Atrans_4         = Rtrans_4;
dxdt_APla             = RPla;                                        
dxdt_APlas_free       = RPlas_free;                                  
dxdt_AAm              = RAm;                                          
dxdt_AUC_CPlas        = CPlas;
dxdt_AUC_CPlas_Fet    = CPlas_Fet;

// #+ ODE equation for fetus compartment
dxdt_ARest_Fet        = RRest_Fet;                                 
dxdt_AL_Fet           = RL_Fet;                                 
dxdt_APlas_Fet_free   = RPlas_Fet;                         

// #+ Virtural compartment; input dose
// dxdt_ADOSE         = 0;

// #+ Mass Balance check (pregnant rat)
// double ATissue     = AF + AM + ARest + APlas_free + AKb + AL + APla + AFil + APTC + AST + ASI;
// double ALoss       = Aurine + Atrans_1 - Atrans_2 + Afeces;
// double ATotal      = ATissue + ALoss;
// double Mbal        = ADOSE - ATotal; 

// #+ Mass Balance check (Fetus)
// double ATissueF    = APlas_Fet_free + ARest_Fet + AL_Fet;
// double ALossF      = Atrans_2 + Atrans_3 - Atrans_4;
// double ATotalF     = ATissueF + ALossF;
// double DoseF       = Atrans_1;
// double MbalF       = DoseF - ATotalF; 

$TABLE
capture Plasma        = CPlas;
capture Liver         = CL;
capture Liver_Fet     = CL_Fet;
capture Plasma_Fet    = CPlas_Fet;
'

LRatPBPK.code <- '
$PROB
# Lactational PFOS PBPK model for SD rat
- Author: Wei-Chun Chou
- Date  : March, 2020
- Dam model strucutre: GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, rest of body 
- Pup model structure: Gut, liver, Kidney, Filtrate, PTC, rest of body
- Note1 : Initial physiological parameters and conditions is matched with the values at GD 22 in gestational PBPK models
- Note2 : Growth equation of physioloical parameters values was taken from Loccisano et al. (2012), yoon et al. (2011) 

$PARAM @annotated
// #+ Dam Parameters
BW0                  : 0.247    : kg,                  Body weight on PND0 (GD22)                    ; Value obtained from Shireley, 1984
Htc                  : 0.46     : Unitless,            Hematocrit                                    ; Value obtained from Davies 1993
QLC0                 : 0.196    : %QC,                 Fraction blood flow to liver on PND0          ; Value obtained Loccisano et al., 2012
QKC0                 : 0.1161   : %QC,                 Fraction blood flow to kidney on PND0         ; Value obtained from Loccisano et al., 2012
QMC0                 : 0.0887   : %QC,                 Fraction blood flow to mammary on PND0        ; Value obtained from Loccisano et al., 2012
QFC                  : 0.07     : %QC,                 Fraction of blood flow of fat tissue          ; Value obtained from Loccisano et al., 2012
QCl0                 : 28.645   : L/hr/kg,             Cardiac Output on PND0                        ; Value obtained from Loccisano et al., 2012
VMilk                : 0.002    : L,                   Volume of milk compartment; Fisher 1990       ; Value obtained from Loccisano et al., 2012
VFC0                 : 0.1245   : %BW,                 Volume of fat tissue on PND0                  ; Value obtained from Loccisano et al., 2012
VMC0                 : 0.049    : %BW,                 Volume of mammary tissue on PND0              ; Value obtained from Loccisano et al., 2012
VLC0                 : 0.046    : %BW,                 Volume of liver tissue on PND0                ; Value obtained from Loccisano et al., 2012
VKC0                 : 0.011    : %BW,                 Volume of kidney tissue on PND0               ; Value obtained from Loccisano et al., 2012
VPlasC               : 0.0466   : L/kg BW,             Fraction of plasma volume                     ; Value obtained from Davies 1993
VFilC                : 8.4e-4   : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obtained from Worley et al., 2017
VPTCC                : 1.35e-4  : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein              : 2.0e-6   : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtained from Addis et al., 1936
PL                   : 3.66     : Unitless,            Liver-to-plasma partition coefficient         ; Value was optimized by Chou and Lin, 2019 
PK                   : 0.80     : Unitless,            Kidney-to-plasma partition coefficient        ; Value obtained from Loccisano et al., 2012
PM                   : 0.16     : Unitless,            Mammary gland-to-plasma partition coefficient ; Value obtained from Loccisano et al., 2012
PF                   : 0.13     : Unitless,            Fat-to-plasma partition coefficient           ; Value obtained from Loccisano et al., 2012
PRest                : 0.26     : Unitless,            Rest of body partition coefficient            ; Value was optimized by Chou and Lin, 2019 
PMilkM               : 1.9      : Unitless,            Milk-to-mammary gland partition coefficient   ; Value obtained from  (mouse, Fenton 2009)
PMilkP               : 0.11     : Unitless,            Milk-to-plasma partition coefficient          ; Value obtained from  (rat, Hinderliter 2005)
PAMilkC              : 0.5      : L/h/kg,              Permeability area cross product (mammary to milk); Value obtained from Loccisano et al., 2012; assumed to be same as PFOA
MW                   : 500.126  : g/mol,               PFOS molecular mass                           ; Value obtained from Worley and Fisher, 2015
Free                 : 0.022    : Unitless,            Female free fraction                          ; Value obtained from Loccisano et al., 2012 
KbileC               : 0.0026   : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value was optimized by Chou and Lin, 2019 
KurineC              : 1.60     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
K0C                  : 1.000    : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KabsC                : 2.12     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KunabsC              : 7.05e-5  : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
KMilk0               : 0.0268   : L/hr/kg,             Zero-order Milk suckling rate (for one individual pup); from Loccisano et al., 2012 
GFRC                 : 41.04    : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value obtained from Corley, 2005
GEC                  : 1.4      : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value obtained from Yang et al., 2013
Vmax_baso_invitro    : 393.45   : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value obtained from Worley and Fisher, 2015; assumed to be same as PFOA
Km_baso              : 27.2     : mg/L,                Km of basolateral transporter                 ; Value calculated from Worley et al., 2007
Vmax_apical_invitro  : 1808     : pmol/mg protein/min, Vmax of apical transporter                    ; Value optimzed from Chou and Lin, 2019 
Km_apical            : 278      : mg/L,                Km of apical transporter                      ; Value optimzed from Chou and Lin, 2019
RAFbaso              : 1.90     : Unitless             Relative activity factor                      ; Value optimzed from Chou and Lin, 2019 
RAFapi               : 4.15     : Unitless             Relative acitivty factor                      ; Value optimzed from Chou and Lin, 2019 
Kdif                 : 5.1e-4   : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value optimzed from Chou and Lin, 2019 
KeffluxC             : 2.09     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimzed from Chou and Lin, 2019  


// #+ Pups Parameters
N                    : 8        : Unitless,            Number of pup                                 ; Value obtained from Loccisano et al. (2012)
Free_pup             : 0.022    : Unitless,            Free fraction                                 ; Value was assumed to be same with dam              
PL_pup               : 3.66     : Unitless,            Liver-to-plasma partition coefficient         ; Value was assumed to be same with dam   
PK_pup               : 0.80     : Unitless,            Kidney-to-plasmaPartition coefficient         ; Value was assumed to be same with dam   
PRest_pup            : 0.22     : Unitless,            Rest of body-to-plasmaPartition coefficient   ; Value was assumed to be same with dam 
KabsC_pup            : 2.12     : 1/h,                 Oral absorption rate constant
KbileC_pup           : 0.0026   : 1/(h*BW^-0.25),      Biliary elimination rate constant             ; Value collected from the value of PFOA from Loccisano et al. (2012);
KurineC_pup          : 1.6      : 1/(h*BW^-0.25),      Urinary elimination rate constant             ; Value was assumed to be same with dam 
Kdif_pup_0           : 0.001    : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value was assumed to be same with dam 
KeffluxC_p           : 2.09     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value was assumed to be same with dam 
Vmax_baso_invitro_p  : 393.45   : pmol/mg,             Protein/min, Vmax of basolateral transporter  ; Value was assumed to be same with dam 
Km_baso_p            : 27.2     : mg/L,                Km of basolateral transporter                 ; Value was assumed to be same with dam 
Vmax_apical_invitro_p : 1808    : pmol/mg protein/min, Vmax of apical transporter                    ; Value was assumed to be same with dam 
Km_apical_p          : 278      : mg/L,                Km of apical transporter                      ; Value was assumed to be same with dam 

// #+ Growth parameters for Body weight, Liver, Kidney and QC in SD rat pup (Mirfazaelian et al., 2007)
Wt0                  : 0.00731  : kg,                  Body weight at birth (GD22); for one pup
K                    : 63.21    : days,                Half maximal growth for pup body weight
g                    : 2.01     : Unitless,            Hill coefficient for body weight growth
Wtmax                : 0.52     : kg,                  Maximal body weight
Wt_LIV0              : 2.96e-4  : kg,                  Liver weight at birth (GD22); for one pup
K_LIV                : 43.49    : days,                Half maximal growth for pup liver weight
g_LIV                : 2.76     : Unitless,            Hill coefficient for liver weight growth
Wtmax_LIV            : 1.54e-2  : kg,                  Maximal liver weight
Wt_KID0              : 6.01e-5  : kg,                  Kidney weight at birth (GD22); for one pup
K_KID                : 50.83    : Days,                Half maximal growth for pup kidney weight
g_KID                : 1.78     : Unitless,            Hill coefficient for kidney weight growth
Wtmax_KID            : 3.8e-3   : kg,                  Maximal Kidney weight
QCmax                : 8.72     : L/hr,                Maximal cardia output
BW50                 : 0.189    : kg,                  50% of bod weight


$MAIN
// #+ Defined time parameters
// #+ PND            : postnatal day

double PND           = TIME/24;

// #+ Growth equations for dam parameters; all equations from OFlaherty, et al. 1992 and Gentry, et al. 2002.
// #+ BW             : day; Calcuate the days from GD0 to PND 
// #+ QCl            : L/h/kg; Cardiac index for total BW (Hanwell & Linzell 1997 (PND1-22) and Dowell 1997) 
// #+ QMC            : Unitless; Fraction blood flow to Mammary gland 
// #+ QLC            : Unitless; Fraction blood flow to liver tissue
// #+ QKC            : Unitless; Fraction blood flow to Kidney tissue 
// #+ VFC            : Unitless; Fraction volume of fat tissue  
// #+ VMC            : Unitless; Fraction volume of mammary gland 
// #+ VLC            : Unitless; Fraction volume of liver 
// #+ VKC            : Unitless; Fraction volume of kidney 
// #+ KMilkC         : L/hr/kg BW; milk production is assumed as suckling rate; for individual pup 

double BW            = PND > 0? 0.0021*PND + BW0 : BW0;
double QCl           = PND > 0? 0.0123*pow(PND,3) - 0.4059*pow(PND,2) + 3.9661*PND + QCl0 : QCl0;
double QMC           = PND > 0? (-0.0001)*pow(PND,2) + 0.0051*PND + QMC0 : QMC0;
double QLC           = PND > 0? 0.0079*PND + 0.2208 : QLC0;
double QKC           = PND > 0? (-0.0014)*PND + QKC0 : QKC0;
double VFC           = PND > 16? 0.07: (-0.0012)*pow(PND,2) + 0.0162*PND + 0.1245;
double VMC           = PND > 0? (-1e-5)*pow(PND,3) + 0.0004*pow(PND,2) + 0.0027*PND + 0.049 : VMC0;
double VLC           = PND > 0? 0.0384*pow(PND,(0.096)) : VLC0;
double VKC           = PND > 0? 0.0685*pow(PND,(0.0514)) : VKC0;
double KMilkC        = PND > 0? (-7e-6)*pow(PND,3) + 0.0003*pow(PND,2) - 0.0032*PND + KMilk0 : KMilk0;

// #+ Dam physiological parameters
// #+ Blood flows
// #+ QC1            : L/h,      Dam cardiac output 
// #+ QPlas          : Unitless, Adjustment for plasma flow using hematocrit 
// #+ QL             : L/h,      Plasma flow to liver 
// #+ QM             : L/h,      Plasma flow to mammary gland  
// #+ QK             : L/h,      Plasma flow to kidney   
// #+ QF             : L/h,      Plasma flow to fat tissue   
// #+ QRest          : L/h,      Plasma flow to rest of body   
// #+ QBal           : L/h,      Balance check for plasma flow   

double QC1           = QCl*BW;
double QPlas         = (1-Htc);
double QC            = QC1*QPlas;
double QL            = QLC*QC;
double QM            = QMC*QC;
double QK            = QKC*QC;
double QF            = QFC*QC;
double QRest         = QC - (QL + QK + QM + QF);
double QBal          = QC - (QRest + QL + QK + QM + QF);

// #+ Tissue volumes in dam
double VPlas         = VPlasC*BW;
double VM            = VMC*BW;
double VL            = VLC*BW;
double VK            = VKC*BW;
double VF            = VFC*BW;

// #+ Volumes of kidney and related compartmnet
double MK            = VKC*BW*1000;                                                                              
double ML            = VLC*BW*1000;                                    
double VPTC          = MK*VPTCC;                                        
double VKb           = VK*0.16;                            
double VFil          = VFilC*BW;                                     
double VRest         = (0.93*BW) - (VL + VK + VF + VM + VPlas);
double VBal          = (0.93*BW) - (VRest + VL + VK + VF + VM + VPlas);

// #+ Kinetics parmaters 
// #+ GFR            : L/h, Glomerular Filtration Rate (GFR); GFR was assumed 50% of renal plasma flow
// #+ GE             : 1/h, Gastric emptying rate constant
// #+ K0             : 1/h, Rate of uptake from the stomach into the liver
// #+ Kbile          : 1/h, Biliary elimination rate constant, liver to feces storage
// #+ Kurine         : 1/h, Urinary elimination rate constant
// #+ Kabs           : 1/h, Rate constant of absorption of PFOS from small intestine to liver
// #+ Kunabs         : 1/h, Rate constant of unabsorbed dose to appear in feces
// #+ PTC            : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC           : g, mass of the proximal tubule cells (assuming density 1 kg/L) 
// #+ Vmax_basoC     : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_baso      : mg/h/kg, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC   : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_apical    : mg/h/kg, Vmax of apical transporters in in vitro studies
// #+ Kefflux        : 1/h, Efflux clearance rate from PTC to blood 

double GFR           = GFRC*(VKC/1000);
double GE            = GEC*pow(BW,(-0.25));                                   
double K0            = K0C*pow(BW,(-0.25));                                   
double Kbile         = KbileC*pow(BW,(-0.25));
double Kurine        = KurineC*pow(BW,(-0.25));
double Kabs          = KabsC*pow(BW,(-0.25));
double Kunabs        = KunabsC*pow(BW,(-0.25));
double PTC           = VKC*6e7*1000;                                       
double MPTC          = VPTC*1000;                                         
double Vmax_basoC    = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_baso     = Vmax_basoC*pow(BW,0.75);                          
double Vmax_apicalC  = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_apical   = Vmax_apicalC*pow(BW,0.75);                      
double Kefflux       = KeffluxC*pow(BW,(-0.25));                         


// #+ Growth equations of tissues for the pup
// #+ BW_pup,        : kg, Body weight for one pup 
// #+ BW_pup_W,      : Kg, Body weight for whole pup
// #+ QC_pup_i,      : L/h, Plasma flow to one pup 
// #+ QC_pup_W,      : L/h, Plasma flow to whole pup 
// #+ QLC_pup,       : L/h, Plasma flow to pup liver  
// #+ QKC_pup,       : L/h, Plasma flow to pup kidney
// #+ Htc_pup,       : L/h, Hematocrit of pup; equation from fit to data; data from Garcia, 1957 and from Loccisano et al. (2012)
// #+ VPlasC_pup,    : %, Fraction of pup plasma volume
// #+ KMilk1,        : L/h, Milk production rate for one pup
// #+ KMilk,         : L/h, Milk production rate for whole litter
// #+ PAMilk,        : L/h, Diffusion rate from mammary tissue to milk
// #+ VL_pup,        : L, Liver volume for one pup
// #+ VK_pup,        : L, Kidney volume for one pup
// #+ VL_pup_w,      : L, Liver volume for whole litter
// #+ VK_pup_w,      : L, Kidney volume for whole litter
// #+ MK_pup         : g, Kidney weight in gram for one pup
// #+ ML_pup         : g, Liver weight in gram for one pup
// #+ VPTC_pup_w     : L, Volume of proximal tubule cells for whole pups
// #+ VKb_pup_W      : L, Volume of kidney serum for whole pups
// #+ VFil_pup_W     : L, Volume of filtrate for whole pups
// #+ VPlas_pup_W    : L, Volume of plasma for whole pups
// #+ VRest_pup_W    : L, Volume of rest of body for whole pups
// #+ VBal_pup       : L, Balance check for pup volume

double BW_pup        = (Wt0*pow(K, g) + Wtmax*(pow(PND,g)))/(pow(K,g)+(pow(PND,g))); 
double BW_pup_W      = BW_pup*N;
double QC_pup_i      = (QCmax*BW_pup)/(BW50 + BW_pup);    
double QC_pup_W      = QC_pup_i*N;                       
double QLC_pup       = (4e-6)*pow(PND, 3) - 0.0005*pow(PND, 2) + 0.0175*PND + 0.0571;
double QKC_pup       = (-4e-5)*pow(PND, 2) + 0.0033*PND + 0.0126;
double Htc_pup       = (-0.0136)*PND + 0.3814;
double VPlasC_pup    = (-6e-6)*pow(PND, 3) + (1e-4)*pow(PND,2) + 0.0006*PND + 0.0454;
double KMilk1        = KMilkC*BW_pup;
double KMilk         = KMilkC*BW_pup_W;
double PAMilk        = PAMilkC*(pow(BW_pup,0.75))*N;
double VL_pup        = (Wt_LIV0*(pow(K_LIV,g_LIV)) + Wtmax_LIV*(pow((PND),g_LIV)))/(pow(K_LIV,g_LIV) + (pow((PND),g_LIV)));
double VL_pup_W      = VL_pup*N;
double VK_pup        = (Wt_KID0*pow(K_KID,g_KID) + Wtmax_KID*(pow((PND),g_KID)))/(pow(K_KID,g_KID) + (pow(PND,g_KID)));
double VK_pup_W      = VK_pup*N;
double MK_pup        = VKC*BW_pup*1000;                                      
double ML_pup        = VLC*BW_pup*1000;  
double VPTC_pup_W    = MK_pup*VPTCC*N;     
double VKb_pup_W     = VK_pup_W*0.16*N;    
double VFil_pup_W    = VFilC*BW_pup_W;   
double VPlas_pup_W   = VPlasC_pup*BW_pup_W;
double VRest_pup_W   = (0.92*BW_pup_W) - (VL_pup_W + VPlas_pup_W + VK_pup_W + VPTC_pup_W + VFil_pup_W);       
double VBal_pup      = (0.92*BW_pup_W) - (VL_pup_W + VPlas_pup_W + VK_pup_W + VRest_pup_W + VPTC_pup_W + VFil_pup_W);    

// #+ Pup scaled rat constant
// #+ GFR_pup_W      : L/h, Glomerular Filtration Rate (GFR) for whole pup;
// #+ Kbile_pup_W    : 1/h, Biliary elimination rate for whole pup
// #+ Kurine_pup_W   : 1/h, Urinary elimination rate for whole pup

double GFR_pup_W     = GFRC*(MK_pup/1000)*N;
double Kbile_pup_W   = KbileC_pup*pow(BW_pup,(-0.25))*N;
double Kurine_pup_W  = KurineC_pup*pow(BW_pup,(-0.25))*N;
double Kabs_pup_W    = KabsC_pup*pow(BW,(-0.25))*N;

// #+ Blood flows for pups
// #+ QPlas_pup,     : L/h, Plasma flow for whole pup 
// #+ QC_pup,        : L/h, Adjust pup CO for plasma flow
// #+ QL_pup,        : L/h, Liver blood flow for pup
// #+ QK_pup,        : L/h, Kidney bloow flow for pup
// #+ QRest_pup,     : L/h, Rest of body bloow flow for pup
// #+ QBal_pup,      : L/h, Balance check for bloow flow for pup

double QPlas_pup       = (1-Htc_pup);
double QC_pup          = QC_pup_W * QPlas_pup;   
double QL_pup          = QC_pup*QLC_pup;
double QK_pup          = QC_pup*QKC_pup;
double QRest_pup       = QC_pup - (QL_pup + QK_pup);
double QBal_pup        = QC_pup - (QL_pup + QK_pup + QRest_pup);

// #+ Kinetics parameters of kidney for PUP
double Vmax_basoC_p    = (Vmax_baso_invitro_p*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_apicalC_p  = (Vmax_apical_invitro_p*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_baso_pup   = Vmax_basoC_p*pow(BW_pup,0.75);                          
double Vmax_apical_pup = Vmax_apicalC_p*pow(BW_pup,0.75);                     
double Kefflux_pup     = KeffluxC_p*pow(BW_pup,(-0.25));                         
double Km_baso_pup     = Km_baso_p;
double Km_apical_pup   = Km_apical_p;
double Kdif_pup        = Kdif_pup_0;


$INIT @annotated
// #+ Set up the initial concentration; the intial concentration was obtained from the PFOS concentration at GD22 from gestatnioal model; Depedent on study desing of animal studies
// #+ Dam compartment
ADOSE           : 0   : mg, Amount of PFOS; input dose; virtual compartment
APlas_free      : 0   : mg, Amount of PFOS in the plasma 
APTC            : 0   : mg, Amount of PFOS in proximal tubule cells
AFil            : 0   : mg, Amount of PFOS in filtrate
Aurine          : 0   : mg, Amount of PFOS in urine
AKb             : 0   : mg, Amount of PFOS in the kidney blood 
ARest           : 0   : mg, Amount of PFOS in rest of body
Afeces          : 0   : mg, Amount of PFOS in the feces  
AL              : 0   : mg, Amount of PFOS in the liver
AM              : 0   : mg, Amount of PFOS in the mammary gland
AF              : 0   : mg, Amount of PFOS in fat compartment
A_baso          : 0   : mg, Amount of PFOS in the baso subcompartment 
A_apical        : 0   : mg, Amount of PFOS in the apical subcompartment
Adif            : 0   : mg, Amount of PFOS in the diffusion virtual compartment
Aefflux         : 0   : mg, Amount of PFOS in the efflux virtual compartment; simulation of PFOS from PTC to blood 
ASI             : 0   : mg, Amount of PFOS in the small intestine compartment
AST             : 0   : mg, Amount of PFOS in the stomach compartment
AabsST          : 0   : mg, Amount of absorbed PFOS in the stomach compartment
AabsSI          : 0   : mg, Amount of absorbed PFOS in the small intestine compartment
Atrans          : 0   : mg, Amount of PFOS from milk to pup
AMilk           : 0   : mg, Amount of PFOS in milk comaprtmnet

// #+ Pup compartment
Abaso_pup       : 0   : mg, Amount of PFOS in the pup rest of body compartment
Aapical_pup     : 0   : mg, Amount of PFOS in the pup plasma compartment
Adif_pup        : 0   : mg, Amount of PFOS in the pup diffusion virtual compartment
Aefflux_pup     : 0   : mg, Amount of PFOS in the efflux virtual compartment; simulation of PFOS from PTC to blood 
APTC_pup        : 0   : mg, Amount of PFOS in proximal tubule cells
AFil_pup        : 0   : mg, Amount of PFOS in the kidney blood 
AKb_pup         : 0   : mg, Amount of PFOS in the kidney blood
Aurine_pup      : 0   : mg, Amount of PFOS in urine
ARest_pup       : 0   : mg, Amount of PFOS in rest of body compartment
AGI_pup         : 0   : mg, Amount of PFOS in the small intestine compartment
AL_pup          : 0   : mg, Amount of PFOS in the small intestine compartment
Afeces_pup      : 0   : mg, Amount of PFOS in the small intestine compartment
APlas_free_pup  : 0   : mg, Amount of PFOS in the plasma 
AUC_CPlas       : 0   : mg/L*hr, Area under curve of PFOS in maternal plasma
AUC_CPlas_pup   : 0   : mg/L*hr, Area under curve of PFOS in pup plasma

$ODE
// #+ Concentrations in the tissues and in the venous plasma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for dam; CX indicate the concentration in the tissue (e.g., CL); CVX represent the concentration of chemical in venous plasma leaving tissue 
double CPlas_free     = APlas_free/VPlas;
double CPlas          = CPlas_free/Free;    
double CKb            = AKb/VKb;
double CPTC           = APTC/VPTC;                                      
double CFil           = AFil/(VFil + 1E-7);                                      
double CMilk          = AMilk/VMilk;
double CL             = AL/ VL;
double CM             = AM/ VM;
double CF             = AF/ VF;
double CRest          = ARest/VRest;
double CVL            = CL/ PL;                                            
double CVK            = AKb/VKb;
double CVM            = CM/ PM;
double CVF            = CF/ PF;
double CVRest         = CRest/ PRest;

// #+ Concentrations for pups
double CPlas_free_pup = APlas_free_pup/(VPlas_pup_W + 1E-7);
double CPlas_pup      = CPlas_free_pup/Free_pup;    
double CL_pup         = AL_pup/ (VL_pup_W + 1E-7);
double CKb_pup        = AKb_pup/(VKb_pup_W + 1E-7);
double CPTC_pup       = APTC_pup/(VPTC_pup_W + 1E-7);                                      
double CFil_pup       = AFil_pup/(VFil_pup_W + 1E-7);  
double CRest_pup      = ARest_pup/(VRest_pup_W + 1E-7);
double CVL_pup        = CL_pup/ PL_pup;
double CVK_pup        = CKb_pup;
double CVRest_pup     = CRest_pup /PRest_pup;   

// #+ Equation for estimation the rate of each compartment in the PBPK model for dam
// #+ RA_baso         : mg/h, Rate of basolateral transporters
// #+ RA_apical       : mg/h, Rate of apical transporter
// #+ Rdif            : mg/h, Rate of diffusion from blood into the PTC
// #+ RAefflux        : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI             : mg/h, Rate of clerance (CL) to via glomerular filtration (GFR)
// #+ RPTC            : mg/h, Rate of change in PTC
// #+ RFil            : mg/h, Rate of change in fil
// #+ RKb             : mg/h, Rate of change in Kidney-serum compartment
// #+ RST             : mg/h, Rate of change in stomach compartment
// #+ RSI             : mg/h, Rate of change in small intestines
// #+ RabsST          : mg/h, Rate of change of absorption in Stomach
// #+ RabsSI          : mg/h, Rate of change of absorption in small intestines
// #+ RL              : mg/h, Rate of change in liver compartment
// #+ RF              : mg/h, Rate of chnage in fat comaprtment
// #+ RM              : mg/h, Rate of change in mammary tissues
// #+ RRest           : mg/h, Rate of change in rest of body
// #+ RPlas_free      : mg/h, Rate of free PFOS change in the dam plasma

double RA_baso        = (Vmax_baso*CKb)/(Km_baso + CKb);                
double RA_apical      = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif           = Kdif*(CKb - CPTC);                               
double RAefflux       = Kefflux*APTC;                                
double RPTC           = Rdif + RA_apical + RA_baso - RAefflux;            
double RCI            = CPlas*GFR*Free;                              
double RFil           = RCI - RA_apical - AFil*Kurine;           
double RKb            = QK*(CPlas-CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RM             = QM*(CPlas - CVM)*Free - PAMilk*(CVM*Free - CMilk/PMilkM);
double Rtrans         = KMilk*CMilk; 
double RMilk          = PAMilk*(CVM*Free - CMilk/PMilkM) - Rtrans;
double RF             = QF*(CPlas - CVF)*Free;                            
double RRest          = QRest*(CPlas - CVRest)*Free;                         
double RST            = -K0*AST - GE*AST + AFil_pup*Kurine_pup_W;
double RabsST         = K0*AST;
double RSI            = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsSI         = Kabs*ASI;
double RL             = QL*(CPlas-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;   
double Rurine         = Kurine*AFil;                                   
double Abile          = Kbile*AL;
double Rfeces         = Kbile*AL + Kunabs*ASI;                         
double RPlas_free     = (QF*CVF*Free) + (QM*CVM*Free) + (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) - (QC*CPlas*Free) + RAefflux;

// #+ ODE equation for the dam compartments
dxdt_A_baso           = RA_baso;                                           
dxdt_A_apical         = RA_apical;                                       
dxdt_Adif             = Rdif;                                              
dxdt_Aefflux          = RAefflux;                                       
dxdt_APTC             = RPTC;                                            
dxdt_AFil             = RFil;                                            
dxdt_AKb              = RKb;                                                
dxdt_AM               = RM;
dxdt_Atrans           = Rtrans;
dxdt_AMilk            = RMilk;
dxdt_AF               = RF;                                            
dxdt_ARest            = RRest;                                           
dxdt_AST              = RST;
dxdt_AabsST           = RabsST;
dxdt_ASI              = RSI;
dxdt_AabsSI           = RabsSI;
dxdt_AL               = RL;                                                  
dxdt_APlas_free       = RPlas_free;                                  
dxdt_Aurine           = Rurine;                                          
dxdt_Afeces           = Rfeces;                                          
dxdt_AUC_CPlas        = CPlas;
dxdt_AUC_CPlas_pup    = CPlas_pup;

// #+ Virtual compartment for estmating input dose
dxdt_ADOSE            = 0;

// #+ Equation for pup compartment
// #+ RA_baso_pup     : mg/h, Rate of basolateral transporters
// #+ RA_apical_pup   : mg/h, Rate of apical transporter
// #+ Rdif_pup        : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux_pup    : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI_pup         : mg/h, Rate of clerance (CL) to via glomerular filtration (GFR)
// #+ RPTC_pup        : mg/h, Rate of change in PTC
// #+ RFil_pup        : mg/h, Rate of change in Fil
// #+ RKb_pup         : mg/h, Rate of change in Kidney-serum compartment
// #+ RGI_pup         : mg/h, Rate of change in GI tract compartment
// #+ RL_pup          : mg/h, Rate of change in liver compartment
// #+ RRest_pup       : mg/h, Rate of change in rest of body
// #+ RPlas_free_pup  : mg/h, Rate of free PFOS change in the pup plasma

double RAbaso_pup     = (Vmax_baso_pup*CKb_pup)/(Km_baso_pup + CKb_pup);             
double RAapical_pup   = (Vmax_apical_pup*CFil_pup)/(Km_apical_pup + CFil_pup);        
double Rdif_pup       = Kdif_pup*(CKb_pup - CPTC_pup);                            
double RAefflux_pup   = Kefflux_pup*APTC_pup;                                 
double RCI_pup        = CPlas_pup*GFR_pup_W*Free_pup;                             
double RPTC_pup       = Rdif_pup + RAapical_pup + RAbaso_pup - RAefflux_pup;            
double RFil_pup       = RCI_pup - RAapical_pup - AFil_pup*Kurine_pup_W;           
double RKb_pup        = QK_pup*(CPlas_pup - CVK_pup)*Free_pup - RCI_pup - Rdif_pup - RAbaso_pup; 
double RRest_pup      = QRest_pup*(CPlas_pup - CVRest_pup)*Free_pup;
double RGI_pup        = Rtrans - Kabs_pup_W*AGI_pup;
double RL_pup         = QL_pup*(CPlas_pup - CVL_pup)*Free_pup + Kabs_pup_W*AGI_pup - Kbile_pup_W*AL_pup*Free_pup; 
double Rfeces_pup     = Kbile_pup_W*AL_pup*Free_pup;
double Rurine_pup     = AFil_pup*Kurine_pup_W;
double RPlas_free_pup = (QRest_pup*CVRest_pup*Free_pup) + (QL_pup*CVL_pup*Free_pup) + (QK_pup*CVK_pup*Free_pup) - (QC_pup*CPlas_pup*Free_pup) + RAefflux_pup;

// #+ ODE equation for pup compartment
dxdt_Abaso_pup        = RAbaso_pup;                                     
dxdt_Aapical_pup      = RAapical_pup;                                   
dxdt_Adif_pup         = Rdif_pup;                                         
dxdt_Aefflux_pup      = RAefflux_pup;                                     
dxdt_APTC_pup         = RPTC_pup;                                            
dxdt_AFil_pup         = RFil_pup;     
dxdt_AKb_pup          = RKb_pup;
dxdt_Aurine_pup       = Rurine_pup;
dxdt_ARest_pup        = RRest_pup;
dxdt_AGI_pup          = RGI_pup;
dxdt_AL_pup           = RL_pup;
dxdt_Afeces_pup       = Rfeces_pup;
dxdt_APlas_free_pup   = RPlas_free_pup;

// #+ Mass Balance check (dam)
double ATissue        = AM + AF + ARest + AKb + AL + APlas_free + AFil + APTC + ASI + AST; 
double ALoss          = Aurine + Afeces + AMilk - Aurine_pup + Atrans;
double ATotal         = ATissue + ALoss;
double Mbal           = ADOSE - ATotal; 

//## Mass Balance check (pup)
double ATissue_pup    = ARest_pup + AKb_pup + AL_pup + APlas_free_pup + AGI_pup + AFil_pup + APTC_pup;
double ALoss_pup      = Aurine_pup + Afeces_pup; 
double Atotal_pup     = ATissue_pup + ALoss_pup;
double Mbal_pup       = Atrans - Atotal_pup;


$TABLE
capture Plasma        = CPlas;
capture Liver         = CL;
capture Plasma_pup    = CPlas_pup;
capture Liver_pup     = CL_pup;
capture Milk          = CMilk;
'

