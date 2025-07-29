PreGHumanPBPK.code <- '
$PROB
# Pre-pregnant PFOS PBPK model 
- Author    : Wei-Chun Chou
- Date      : March, 2020
- Strucutre : GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC
- Note1     : Pysiological parameters were obtained from Brown, 1997; Yoon et al., 2011; Loccisano et al., 2013; Corley, 2005 
- Note2     : Chemical-specific parameters were optimized from Chou and Lin, 2019

$PARAM @annotated
BW                  : 60        : kg,                  Body weight                                   ; Value from Yoon et al. (2011)
Htc                 : 0.44      : Unitless             Hematocrit for human                          ; ICRP Publication 89 (2003)
QCC                 : 16.4      : L/h/kg^0.75,         Cardiac output                                ; Value obtained from yoon et al. (2011); Brown 1997; Forsyth 1968
QLC                 : 0.25      : Unitless,            Fraction blood flow to liver (%QCC)           ; Value obtained from Brown 1997; Fisher 2000
QKC                 : 0.175	    : Unitless,            Fraction blood flow to kidney (%QCC)          ; Value obtained from Brown 1997, Forsyth 1968
QMC                 : 0.027     : Unitless,            Fraction blood flow to Mammary gland (%QCC)   ; Value obtained from yoon et al., 2011
QFC                 : 0.052     : Unitless,            Fraction blood flow to Fat (%QCC)             ; Value obtained from Brown 1997; yoon et al., 2011
VLC                 : 0.026     : Unitless,            Fraction of liver tissue volume (%BW)         ; Value obtained from Brown 1997
VKC                 : 0.004     : Unitless,            Fraction of kidney tissue volume (%BW)        ; Value obtained from Brown 1997
VMC                 : 0.0062    : Unitless,            Fraction of Mammary gland volume (%BW)        ; Value obtained from yoon et al., 2011
VFC                 : 0.214     : Unitless,            Fraction of fat tissue (%BW)                  ; Value was assumed to be 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0428    : L/kg BW,             Fraction of plasma volume (%BW)               ; Value obtained from Davies and Morries, 1993
VFilC               : 8.4e-4    : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obatedin from Worley et al., 2017
VPTCC               : 1.35e-4   : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6    : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtaedin from Addis et al., 1936
PL                  : 2.03      : Unitless,            Liver-to-plasma partition coefficient         ; Value optimized from Chou and Lin, 2019
PK                  : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value optimized from Chou and Lin, 2019
PM                  : 0.16      : Unitless,            Mammary gland-to-plasma partition coefficient ; Value from Loccisano et al., 2013
PF                  : 0.13      : Unitless,            Fat-to-plasma partition coefficient           ; Value from Loccisano et al., 2013
PRest               : 0.20      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value from Loccisano et al., 2013
MW                  : 500.126   : g/mol,               PFOS molecular mass                           ; Value from Worley et al. (2017)
Free                : 0.014     : Unitless,            Free fraction of PFOS in plasma               ; Value optimized from Chou and Lin, 2019
KbileC              : 1.3e-4    : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
KurineC             : 0.096     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
K0C                 : 1.000     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value from Worley et al. (2017); Assumed to be same as PFOA
KabsC               : 2.120     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value from Worley et al. (2017); Assumed to be same as PFOA 
KunabsC             : 7.05e-5   : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value from Worley et al. (2017); Assumed to be same as PFOA  
GFRC                : 27.28     : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value from Corley, 2005
GEC                 : 3.510     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value from Yang et al., 2014
Vmax_baso_invitro   : 479       : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value optimized from Chou and Lin, 2019
Km_baso             : 20.1      : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al. (2007)
Vmax_apical_invitro : 51803     : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimized from Chou and Lin, 2019 
Km_apical           : 64.4      : mg/L,                Km of apical transporter                       ; Value optimized from Chou and Lin, 2019
RAFbaso             : 1         : Unitless             Relative activity factor for baso             ; Value calculated from Worley et al. (2017); Assumed to be same as PFOA 
RAFapi              : 0.001     : Unitless             Relative activity factor for api              ; Value optimized from Chou and Lin, 2019 
Kdif                : 0.001     : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value from Worley et al. (2017); Assumed to be same as PFOA
KeffluxC            : 0.150     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimized from Chou and Lin, 2019

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ YEAR          

double YEAR           = TIME/24*365;                                     

// #+ Physiological parameters: Blood flows
// #+ QC            : L/h,       Cardiac output 
// #+ QK            : L/h,       Blood flows of kidney
// #+ QL            : L/h,       Blood flows of liver  
// #+ QM            : L/h,       Blood flows of mammary gland 
// #+ QF            : L/h,       Blood flows of fat tissue 
// #+ QRest         : L/h,       Blood flows of rest of body  
// #+ QBal          : L/h,       Mass balance check for total blood flows 

double QC           = QCC*pow(BW, 0.75)*(1-Htc);                  
double QK           = QKC*QC;                                     
double QL           = QLC*QC;                                     
double QM           = QMC*QC;                                     
double QF           = QFC*QC;                                     
double QRest        = QC - QK - QL - QM - QF;                          
double QBal         = QC - (QK + QL + QRest + QM + QF);

// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

// #+ Physiological parameters: tissue volume
// #+ VL            : L,         Volume of liver
// #+ VK            : L,         Volume of Kidney
// #+ VKb           : L,         Volume of blood in the kidney  
// #+ VFil          : L,         Volume of filtrate
// #+ VPTC          : L,         Volume of proximal tubule cells (PTCs)
// #+ MK            : g,         Kidney weight                                 ; based on density of kidney = 1.0 g/mL
// #+ ML            : g,         Liver weight                                  ; based on density of liver = 1.05 g/mL, Overmeyer 1987
// #+ VBal          : L,         Mass balance check for total tissue volume

double VL           = VLC*BW;                                     
double VK           = VKC*BW;                                     
double VKb          = VK*0.16;                                   
double VM           = VMC*BW;                                     
double VF           = VFC*BW;                                     
double VPlas        = VPlasC*BW;                               
double VFil         = VFilC*BW;                                  
double VPTC         = VK*VPTCC;                                  
double VRest        = (0.93*BW) - VL - VK - VM - VF - VPlas;  
double MK           = VKC*BW*1000;                                                                            
double ML           = VLC*1.05*BW*1000;                                       
double VBal         = (0.93*BW) - VRest - (VL + VK + VM + VF + VPlas);  

// #+ Kinetic parameters
// #+ Kbile         : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine        : 1/h, Urinary elimination
// #+ Kabs          : 1/h, Rate of absorption of PFOS from small intestine to liver
// #+ Kunabs        : 1/h, Rate of unabsorbed dose to appear in feces
// #+ PTC           : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC          : g, mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC    : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC  : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_baso     : mg/h,
// #+ Vmax_apical   : mg/h,
// #+ Kefflux       : 1/h, Efflux clearance rate from PTC to blood  
// #+ GFR           : L/h, Glomerular filtration rate
// #+ GE            : 1/h, Gasric emptying rate
// #+ K0            : 1/h, Rate of uptake from the stomach into the liver

double Kbile        = KbileC*pow(BW,(-0.25));                  
double Kurine       = KurineC*pow(BW,(-0.25));                 
double Kabs         = KabsC*pow(BW,(-0.25));                    
double Kunabs       = KunabsC*pow(BW,(-0.25));                
double PTC          = VKC*6e7*1000;                              
double MPTC         = VPTC*1000;                                
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);         
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);      
double Vmax_baso    = Vmax_basoC*pow(BW,0.75);             
double Vmax_apical  = Vmax_apicalC*pow(BW,0.75);         
double Kefflux      = KeffluxC*pow(BW,(-0.25));              
double GFR          = GFRC*(MK/1000);                            
double GE           = GEC*pow(BW,(-0.25));                        
double K0           = K0C*pow(BW,(-0.25));                        


$CMT ADOSE A_baso A_apical Adif Aefflux ACI APlas_free AUCCPlas_free APTC AUCCPTC AFil AUCFil Aurine AKb AUCKb ARest AUCCRest AST
AabsST ASI AabsSI Afeces AL AUCCL AF AUCCF AM AUCCM

$ODE
// #+ Concentrations in the tissues and in the venous palsma leaving each of the tissues (Unit: mg/L) 
// #+ CPlas_free  : mg/L, Free PFOS concentration in the plasma
// #+ CPlas       : mg/L, Concentration of total PFOS in the plasma
// #+ CL          : mg/L, Concentration of PFOS in the liver compartment
// #+ CKb         : mg/L, Concentration of PFOS in venous plasma leaving kidney
// #+ CK          : mg/L, Concentration of PFOS in Kidney compartment
// #+ CM          : mg/L, Concentration of PFOS in mammary gland compartment
// #+ CF          : mg/L, Concentration of PFOS in fat compartment
// #+ CRest       : mg/L, Concentration of PFOS in rest of the body
// #+ CPTC        : mg/L, Concentration of PFOS in proximal tubule cells (PTC)
// #+ CFil        : mg/L, Concentration of PFOS in Filtrate (fil)
// #+ CVL         : mg/L, Concentration in the venous blood leaving the liver
// #+ CVK         : mg/L, Concentration in the venous blood leaving the kidney
// #+ CVM         : mg/L, Concentration in the venous blood leaving the mammary gland
// #+ CVF         : mg/L, Concentration in the venous blood leaving the fat
// #+ CVRest      : mg/L, Concentration in the venous blood leaving the rest of body

double CPlas_free = APlas_free/VPlas;                      
double CPlas      = CPlas_free/Free;                               
double CL         = AL/VL;                                      
double CKb        = AKb/VKb;                                   
double CK         = CKb*PK;                                     
double CM         = AM/VM;                                      
double CF         = AF/VF;                                      
double CRest      = ARest/VRest;                             
double CPTC       = APTC/VPTC;                                       
double CFil       = AFil/VFil;                                       
double CVL        = CL/PL;                                     
double CVK        = CKb;                                       
double CVM        = CM/PM;                                     
double CVF        = CF/PF;                                     
double CVRest     = CRest/PRest;                    

// #+ Equation for estimating the rate of compartment in mother PBPK
// #+ Virtural kidney sub-compartment 
// #+ RA_baso     : mg/h, Rate of basolateral transporters
// #+ RA_apical   : mg/h, Rate of apical transporter
// #+ Rdif        : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux    : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI         : mg/h, Rate of clearance(CL) to via glomerular filtration (GFR)
// #+ RPTC        : mg/h, Rate of change in PTC
// #+ RFil        : mg/h, Rate of change in Fil
// #+ RKb         : mg/h, Rate of change in Kidney-serum compartment
// #+ RST         : mg/h, Rate of change in stomach compartment
// #+ RSI         : mg/h, Rate of change in small intestines
// #+ RabsST      : mg/h, Rate of change of absorption in Stomach
// #+ RabsSI      : mg/h, Rate of change of absorption in small intestines
// #+ RL          : mg/h, Rate of change in liver compartment
// #+ RF          : mg/h, Rate of change in fat compartment
// #+ RM          : mg/h, Rate of change in mammary tissues
// #+ RRest       : mg/h, Rate of change in rest of body
// #+ Rurine      : mg/h, Rate of change in urine
// #+ Rbile       : mg/h, Rate of change in bile compartment 
// #+ RPlas_free  : mg/h, Rate of change in the plasma

double RA_baso    = (Vmax_baso*CKb)/(Km_baso + CKb);                                            
double RA_apical  = (Vmax_apical*CFil)/(Km_apical + CFil);                                   
double Rdif       = Kdif*(CKb - CPTC);                                                           
double RAefflux   = Kefflux*APTC;                                                             
double RCI        = CPlas*GFR*Free;                                                                   
double RPTC       = Rdif + RA_apical + RA_baso - RAefflux;                                        
double RFil       = RCI - RA_apical - AFil*Kurine;                                        
double RKb        = QK*(CPlas - CVK)*Free - RCI - Rdif - RA_baso;                               
double RST        = - K0*AST - GE*AST;                                                             
double RSI        = GE*AST - Kabs*ASI - Kunabs*ASI;                                                
double RabsST     = K0*AST;                                                                     
double RabsSI     = Kabs*ASI;
double RL         = QL*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;                                
double RF         = QF*(CPlas - CVF)*Free;                                                               
double RM         = QM*(CPlas - CVM)*Free;   
double RRest      = QRest*(CPlas - CVRest)*Free;                                                      
double Rurine     = Kurine*AFil;                                                               
double Rbile      = Kbile*AL;                                                                    
double Rfeces     = Rbile + Kunabs*ASI; 
double RPlas_free = (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) + (QM*CVM*Free) + (QF*CVF*Free) - (QC*CPlas*Free) + RAefflux;  

// #+ ODE equation 
dxdt_A_baso       = RA_baso;                                                                      
dxdt_A_apical     = RA_apical;                                                                  
dxdt_Adif         = Rdif;                                                                           
dxdt_Aefflux      = RAefflux;                                                                    
dxdt_ACI          = RCI;                                                                             
dxdt_APTC         = RPTC;                                                                           
dxdt_AFil         = RFil;                                                                           
dxdt_AKb          = RKb;                                                                             
dxdt_AST          = RST; 
dxdt_AabsST       = RabsST;                                                                       
dxdt_ASI          = RSI;                                                                            
dxdt_AabsSI       = RabsSI;                                                                       
dxdt_AL           = RL;                                                                               
dxdt_AF           = RF;                                                                               
dxdt_AM           = RM; 
dxdt_ARest        = RRest;                                                                         
dxdt_Aurine       = Rurine;                                          
dxdt_Afeces       = Rfeces;                                          
dxdt_APlas_free   = RPlas_free;                                                               

// #+ Virtural compartment for estimating input dose
dxdt_ADOSE        = 0;

// #+ AUC equation 
dxdt_AUCCPlas_free = CPlas_free;                                                                  
dxdt_AUCCPTC      = CPTC;                                                                       
dxdt_AUCFil       = CFil;                                                                      
dxdt_AUCKb        = CK;                                                                            
dxdt_AUCCRest     = CRest;                                                                     
dxdt_AUCCL        = CL;                                                                            
dxdt_AUCCM        = CM;                                                                            
dxdt_AUCCF        = CF;                                                                            

// #+ Mass Balance check
//double ATissue    = APlas_free + ARest + AKb + AFil + APTC + AL + AM + AF + AST + ASI;
//double ALoss      = Aurine + Afeces;
//double ATotal     = ATissue + ALoss;
//double Mbal       = ADOSE - ATotal;

$TABLE
capture Plasma    = CPlas;
capture Liver     = CL;
capture Kidney    = CK;
'




#///////////////////////////////////////////////////////
GHumanPBPK.code <- '
$PROB
# Gestational PFOS PBPK model for human
- Author : Wei-Chun Chou
- Date   : March, 2020
- Structure: GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, rest of fetuses, fetal liver, amniotic fluid  
- Note1  : Initial physiological parameters and optimized parameters is matched with the values in non-pregnantPBPK models
- Note2  : Growth equation of physiological parameters values and fetus-specific parameters was taken from Loccisano et al. (2013); Yoon et al. (2011); Kapraun et al., 2019 

$PARAM @annotated
// #+  Maternal parameters
BW                  : 60        : kg,                  Body weight                                   ; Value from Yoon et al. (2011)
QCC                 : 16.4      : L/h/kg^0.75,         Cardiac output                                ; Value obtained from yoon et al. (2011); Brown 1997; Forsyth 1968
QLC                 : 0.25      : Unitless,            Fraction blood flow to liver (%QCC)           ; Value obtained from Brown 1997; Fisher 2000
QKC                 : 0.175	    : Unitless,            Fraction blood flow to kidney (%QCC)          ; Value obtained from Brown 1997, Forsyth 1968
QMC                 : 0.027     : Unitless,            Fraction blood flow to Mammary gland (%QCC)   ; Value obtained from yoon et al., 2011
QFC                 : 0.052     : Unitless,            Fraction blood flow to Fat (%QCC)             ; Value obtained from Brown 1997; yoon et al., 2011
VLC                 : 0.026     : Unitless,            Fraction of liver tissue volume (%BW)         ; Value obtained from Brown 1997
VKC                 : 0.004     : Unitless,            Fraction of kidney tissue volume (%BW)        ; Value obtained from Brown 1997
VMC                 : 0.0062    : Unitless,            Fraction of Mammary gland volume (%BW)        ; Value obtained from yoon et al., 2011
VFC                 : 0.214     : Unitless,            Fraction of fat tissue (%BW)                  ; Value was assumed to be 1.5 times of male parameter (0.0723) based on ICRP, 1975.
VPlasC              : 0.0428    : L/kg BW,             Fraction of plasma volume (%BW)               ; Value obtained from Davies 1993
VFilC               : 8.4e-4    : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obatedin from Worley et al., 2017
VPTCC               : 1.35e-4   : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6    : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtaedin from Addis et al., 1936
PL                  : 2.03      : Unitless,            Liver-to-plasma partition coefficient         ; Value optimized from Chou and Lin, 2019
PK                  : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value optimized from Chou and Lin, 2019
PM                  : 0.16      : Unitless,            Mammary gland-to-plasma partition coefficient ; Value from Loccisano et al., 2013
PF                  : 0.13      : Unitless,            Fat-to-plasma partition coefficient           ; Value from Loccisano et al., 2013
PRest               : 0.20      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value from Loccisano et al., 2013
PPla                : 0.41      : Unitless,            Placenta-to-plasma partition coefficient      ; Value from Loccisano et al., 2013
MW                  : 500.126   : g/mol,               PFOS molecular mass                           ; Value from Worley et al. (2017)
Free                : 0.014     : Unitless,            Free fraction of PFOS in plasma               ; Value optimized from Chou and Lin, 2019
KbileC              : 1.3e-4    : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
KurineC             : 0.096     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
K0C                 : 1.000     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value from Worley et al. (2017); Assumed to be same as PFOA
KabsC               : 2.120     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value from Worley et al. (2017); Assumed to be same as PFOA 
KunabsC             : 7.05e-5   : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value from Worley et al. (2017); Assumed to be same as PFOA  
GFRC                : 27.28     : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value from Corley, 2005
GEC                 : 3.510     : 1/(h*BW^0.25),       Gastric emptying constant                     ; Value from Yang et al., 2014
Vmax_baso_invitro   : 479       : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value optimized from Chou and Lin, 2019
Km_baso             : 20.1      : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al. (2007)
Vmax_apical_invitro : 51803     : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimized from Chou and Lin, 2019 
Km_apical           : 64.4      : mg/L,                Km of apical transporter                       ; Value optimized from Chou and Lin, 2019
RAFbaso             : 1         : Unitless             Relative activity factor for baso             ; Value calculated from Worley et al. (2017); Assumed to be same as PFOA 
RAFapi              : 0.001     : Unitless             Relative activity factor for api              ; Value optimized from Chou and Lin, 2019 
Kdif                : 0.001     : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value from Worley et al. (2017); Assumed to be same as PFOA
KeffluxC            : 0.150     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimized from Chou and Lin, 2019
Ktrans1C            : 0.46      : L/h/kg^0.75,         Mother-to-fetus placental transfer rate       ; Value from Loccisano et al., 2013
Ktrans2C            : 1.01      : L/h/kg^0.75,         Fetus-to-mother placental transfer rate       ; Value from Loccisano et al., 2013
Ktrans3C            : 0.006     : L/h/kg^0.75,         Fetus-to-amniotic fluid transfer rate         ; Value from Loccisano et al., 2013 
Ktrans4C            : 0.001     : L/h/kg^0.75,         Amniotic fluid-to-fetus transfer rate         ; Value from Loccisano et al., 2013 
VPlasC_Fet          : 0.0428    : Unitless,            Fraction of fetal plasma volume (%BW)         ; Value was assumed to be same with mother; Loccisano et al. (2013)
Free_Fet            : 0.014     : Unitless,            Free fraction of PFOS in plasma at t = 0      ; Value was assumed to be same with mother
PRest_Fet           : 0.20      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value was assumed to be same with mother
PL_Fet              : 2.03      : Unitless,            Liver-to-plasma partition coefficient in Fetus; Value was assumed to be same with mother

// #+ Sensitivity anlaysis (SA) parametes for growth equations
SA_VPlas            : 1 : Unitless, SA for maternal Volume of plasma (VPlas)
SA_Htc              : 1 : Unitless, SA for maternal Hematocrit (Htc)
SA_GFR              : 1 : Unitless, SA for maternal Glomerular filtration rate (GFR)
SA_VAm              : 1 : Unitless, SA for maternal volume of amniotic fluid volume (VAm)
SA_VPla             : 1 : Unitless, SA for maternal volume of placenta (VPla)
SA_VF_P             : 1 : Unitless, SA for maternal volume of fat during pregnant (VF_P)
SA_VM_P             : 1 : Unitless, SA for maternal volume of mammary during pregnant (VM_P)
SA_QC               : 1 : Unitless, SA for maternal Cardiac output during pregnancy (QC)
SA_QF_P             : 1 : Unitless, SA for maternal Blood flows of fat during prengnacy (QF_P)
SA_QK_P             : 1 : Unitless, SA for maternal Blood flows of kidney during prengnacy (QK_P)
SA_QL_P             : 1 : Unitless, SA for maternal Blood flows of liver during prengnacy (QL_P)
SA_QPla             : 1 : Unitless, SA for maternal Blood flows of placenta (QPla)
SA_Htc_Fet          : 1 : Unitless, SA for fetal Hematocrit (Htc_Fet)
SA_VFet             : 1 : Unitless, SA for volume of fetus (VFet)
SA_VL_Fet           : 1 : Unitless, SA for volume of fetal liver (VL_Fet)

$MAIN
// #+ Time varabiles: Gestational day and age
// #+ GD            : Day,       Gestational days
// #+ GA            : Week,      Gestational age

double GD           = TIME/24;                                     
double GA           = GD/7;                                        

// #+ Growth equations for maternal or fetal physiological parameters; Equation from Kapraun et al. (2019): https://doi.org/10.1371/journal.pone.0215906 
// #+ VPlas         : L,         Volume of plasma                              ; Equation from Kapraun et al., 2019                                         
// #+ Htc           : Unitless,  Hematocrit                                    ; Equation from Kapraun et al., 2019
// #+ GFR           : L/h,       Glomerular filtration rate                    ; Equation from Kapraun et al., 2019; orinral unit is mL/min; using (60/1000) convert unit (mL/min) to L/h
// #+ VAm           : L,         Growth equations for amniotic fluid volume (L); Equation from Kapraun et al., 2019 
// #+ VPla          : L,         Volume of placenta                            ; Equation from Kapraun et al., 2019
// #+ VF_0          : L,         Volume of fat tissue at GD0
// #+ VF_P          : L,         Volume of fat tissue during pregnnacy         ; Equation from Kapraun et al., 2019; 0.95 kg/L is the mean density (Martin et al., 1994)
// #+ VM_0          : L,         Volume of mammary gland at GD0
// #+ VM_P          : L,         Volume of mammary gland during pregnnacy      ; Equation from Loccisano et al. (2013)
// #+ VL            : L,         Volume of liver                               
// #+ VK            : L,         Volume of Kidney
// #+ VKb           : L,         Volume of blood in the kidney  
// #+ VFil          : L,         Volume of filtrate
// #+ VPTC          : L,         Volume of proximal tubule cells (PTC)
// #+ MK            : g,         Kidney weight                                 ; based on density of kidney = 1.0 g/mL
// #+ ML            : g,         Liver weight                                  ; based on density of liver = 1.05 g/mL, Overmeyer 1987
                                         
double VPlas        = SA_VPlas*((1.2406/(1 + exp(-0.31338*(GA - 17.813)))) + 2.4958);
double Htc          = SA_Htc*((39.192 - 0.10562*GA - (7.1045E-4)*pow(GA, 2))/100);
double GFR          = SA_GFR*((113.73 + 3.5784*GA - 0.067272*pow(GA, 2))*(0.06)); 
double VAm          = SA_VAm*((822.34/(1 + exp(-0.26988*(GA - 20.150))))/1000);
double VPla         = GA < 2? 0 : SA_VPla*((-1.7646*GA + 0.91775*pow(GA, 2) - 0.011543*pow(GA, 3))/1000);
double VF_P         = SA_VF_P*((1/0.95)*(17.067 + 0.14937*GA)); 
double VF_0         = BW*VFC;
double VM_0         = BW*VMC;
double VM_P         = SA_VM_P*BW*(((VMC + (0.0065*exp(-7.444868*exp(-0.000678*(GD*24)))))));
double VL           = VLC*BW;
double VK           = VKC*BW;
double MK           = VKC*BW*1000;                                                                            
double ML           = VLC*BW*1.05*1000;                                       
double VKb          = VK*0.16;                            
double VFil         = VFilC*BW;                                      
double VPTC         = MK*VPTCC;  


// #+ Changes of maternal body weight during pregnnacy
// #+ BW_P          : kg,        Maternal BW in pregnant women; Equation from Loccisano et al. (2013)
// #+ BWinc         : kg,        BW increased during prengncy
// #+ VRest         : L,         Volume of rest of body
// #+ VBal          : L,         Mass balance check for total tissue volume

double BW_P         = BW + (VF_P - VF_0) + (VM_P - VM_0) + VPla + VFet + VAm; 
double BWinc        = (VF_P - VF_0) + (VM_P - VM_0) + VPla + VFet + VAm;  
double VRest        = (0.93*BW_P) - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm); 
double VBal         = (0.93*BW_P) - VRest -(VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm);  
                   
// #+ Growth equations for blood flows; Euqations from  Loccisano et al. (2013) or Kapraun et al., 2019
// #+ QC            : L/h,       Cardiac output during pregnancy 
// #+ QC_0          : L/h,       Cardiac output at GD0           
// #+ QF_P          : L/h,       Blood flows of liver during pregnancy   
// #+ QF            : L/h,       Blood flows of fat tissue for non-pregnant women
// #+ QK_P          : L/h,       Blood flows of kidney during prengnacy
// #+ QK            : L/h,       Blood flows of kidney for non-pregnant women
// #+ QL_P          : L/h,       Blood flows of liver during pregnancy
// #+ QL            : L/h,       Blood flows of liver for non-pregnant women 
// #+ QPla          : L/h,       Blood flows of placenta; Equation from Kapraun et al., 2019
// #+ QM            : L/h,       Blood flows of mammary gland for non-pregnant women
// #+ QM_P          : L/h,       Blood flows of mammary gland during pregnnacy
// #+ QRest         : L/h,       Blood flows of rest of body for non-pregnant women 
// #+ QBal          : L/h,       Mass balance check for total blood flows 

double QC_0         = QCC*pow(BW,0.75)*(1-Htc);
double QC           = SA_QC*(QC_0 + 3.2512*GA + 0.15947*pow(GA, 2) - 0.0047059*pow(GA, 3));
double QF_P         = SA_QF_P*(((0.01)*(8.5 + (-0.0175)*GA))*QC);
double QF           = QFC*QC_0;
double QK_P         = SA_QK_P*((0.01)*(17 + (-0.01)*GA)*QC);
double QK           = QKC*QC_0;
double QL_P         = SA_QL_P*((0.01)*(27 + (-0.175)*GA)*QC);
double QL           = QLC*QC_0;
double QPla         = GA < 3.6 ? 0 : SA_QPla*((0.00022)*(GA - 3.6)*(0.4 + 0.29*GA)*QC);
double QM           = QMC*QC_0;
double QM_P         = QM*(VM_P/VM_0);
double QRest        = QC - (QK_P + QL_P + QM_P + QF_P + QPla);
double QBal         = QC - (QK_P + QL_P + QRest + QPla + QM_P + QF_P);

// #+ Kinetic parameters
// #+ Kbile         : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine        : 1/h, Urinary elimination
// #+ Kabs          : 1/h, Rate of absorption of PFOS from small intestine to liver
// #+ Kunabs        : 1/h, Rate of unabsorbed dose to appear in feces
// #+ PTC           : cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
// #+ MPTC          : g, mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC    : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC  : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_baso     : mg/h,
// #+ Vmax_apical   : mg/h,
// #+ Kefflux       : 1/h, Efflux clearance rate from PTC to blood  
// #+ Ktrans_1      : L/h, Rate constant for placental transfer; Mother to fetus
// #+ Ktrans_2      : L/h, Rate constant for placental transfer; fetus to Mother
// #+ Ktrans_3      : L/h, Amniotic fluid transfer rate; fetus to fluid
// #+ Ktrans_4      : L/h, Amniotic fluid transfer rate; fluid to fetus
// #+ GE            : 1/h, Gasric emptying time
// #+ K0            : 1/h, Rate of uptake from the stomach into the liver

double Kbile        = KbileC*pow(BW_P,(-0.25));                   
double Kurine       = KurineC*pow(BW_P,(-0.25));                
double Kabs         = KabsC*pow(BW_P,(-0.25));                     
double Kunabs       = KunabsC*pow(BW_P,(-0.25));                
double PTC          = VKC*6e7*1000;                                       
double MPTC         = VPTC*1000;                                         
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_baso    = Vmax_basoC*pow(BW_P,0.75);                          
double Vmax_apical  = Vmax_apicalC*pow(BW_P,0.75);                      
double Kefflux      = KeffluxC*pow(BW_P,(-0.25));                         
double Ktrans_1     = Ktrans1C*(pow(VFet,0.75));                      
double Ktrans_2     = Ktrans2C*(pow(VFet,0.75));                     
double Ktrans_3     = Ktrans3C*(pow(VFet,0.75));                      
double Ktrans_4     = Ktrans4C*(pow(VFet,0.75));                     
double GE           = GEC*pow(BW_P,(-0.25));                        
double K0           = K0C*pow(BW_P,(-0.25));                        

// #+ Fetus physiological parameters
// #+ Htc_Fet       : Unitness,  Fetal Hematocrit; Equation from Kapraun et al., 2019; (0.5) fetal hematocrit (same as newborn); Sisson, et al 1959
// #+ VFet          : L,         Volume or mass (kg) of a human fetus; Equation obtained from Kapraun et al., 2019
// #+ VPlas_Fet     : L,         Plasma volume of fetal 
// #+ VL_Fet        : L,         Volume of fetal liver; Equation obtained from Kapraun et al., 2019
// #+ VRest_Fet     : L,         Volume of fetal rest of body
// #+ VFetBal       : L,         Mass balance for fetal volume
// #+ QCC_Fet       : L/h/kg,    Cardiac output of Fetus; value obtained from Mendes et al. (2015); Table 3     
// #+ QC_Fet        : L/h,       Cardiac output of Fetus; equation obatained from Loccisano et al. (2013)
// #+ QL_Fet        : L/h,       Blood flows of fetal Liver; Equation modifed from Kapraun et al., 2019 (blood flow in the intra-abdominal umbilical vein "Qpla_f" was assumed as 0)
// #+ QRest_Fet     : L/h,       Blood flows of fetal rest of body
// #+ QBal_Fet      : L/h,       Mass balance check for fetal blood flows

double Htc_Fet      = SA_Htc_Fet*((4.5061*GA - 0.18487*pow(GA,2) + 0.0026766*pow(GA,3))/100);
double VFet         = SA_VFet*((0.0018282*exp((15.12691)*(1-exp(-0.077577*GA))))/1000);
double VPlas_Fet    = VPlasC_Fet*VFet;
double VL_Fet       = SA_VL_Fet*((0.0075*exp(10.68*(1-exp((-0.062)*GA))))/1050);
double VRest_Fet    = (0.93*VFet) - VPlas_Fet - VL_Fet;
double VFetBal      = (0.93*VFet) - (VRest_Fet + VPlas_Fet + VL_Fet);

double QCC_Fet      = 54;
double QC_Fet       = QCC_Fet*VPlas_Fet*(1 - Htc_Fet); 
double QL_Fet       = (6.5/54)*(1 - 26.5/75)*QC_Fet;
double QRest_Fet    = QC_Fet - QL_Fet;
double QBal_Fet     = QC_Fet - (QRest_Fet + QL_Fet);
            
// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

$INIT @annotated
// #+ Set up the initial concentration; the initialconcentration was obtained from the steady-state concentration from non-pregnant model
ADOSE             : 0   : mg, Amount of input dose; virtual compartment
APlas_free        : 0   : mg, Amount of PFOS in the plasma 
APTC              : 0   : mg, Amount of PFOS in the proximal tubule cells
AFil              : 0   : mg, Amount of PFOS in the filtrate
Aurine            : 0   : mg, Amount of PFOS in urine
AKb               : 0   : mg, Amount of PFOS in the kidney blood 
ARest             : 0   : mg, Amount of PFOS in rest of body
Afeces            : 0   : mg, Amount of PFOS in the feces  
AL                : 0   : mg, Amount of PFOS in the liver
AM                : 0   : mg, Amount of PFOS in the mammary gland
AF                : 0   : mg, Amount of PFOS in fat compartment
A_baso            : 0   : mg, Amount of PFOS in the apical subcompartment 
A_apical          : 0   : mg, Amount of PFOS in the apical subcompartment 
Adif              : 0   : mg, Amount of of PFOS in the diffused into the proximal tubule cells (PTC)
Aefflux           : 0   : mg, Amount of PFOS in the efflux virtual compartment; simulation of PFOS from PTC to blood
Atrans_1          : 0   : mg, Amount of PFOS in the placental transfer from mother to fetus  
Atrans_2          : 0   : mg, Amount of PFOS in the placental transfer from fetus to mother 
Atrans_3          : 0   : mg, Amount of PFOS in the Amniotic fluid transfer from Amniotic fluid to fetus  
Atrans_4          : 0   : mg, Amount of PFOS in the Amniotic fluid transfer from fetus to Amniotic fluid 
APla              : 0   : mg, Amount of PFOS in the placenta compartment
ASI               : 0   : mg, Amount of PFOS in the small intestine compartment
AST               : 0   : mg, Amount of PFOS in the stomach compartment
AabsST            : 0   : mg, Amount of absorbed PFOS in the stomach compartment
AabsSI            : 0   : mg, Amount of absorbed PFOS in the small intestine compartment
AAm               : 0   : mg, Amount of PFOS in the Amniotic fluid compartment
AL_Fet            : 0   : mg, Amount of PFOS in the Fetal liver compartment
ARest_Fet         : 0   : mg, Amount of PFOS in the rest of body compartment
APlas_Fet_free    : 0   : mg, Amount of PFOS in the plasma compartment
AUC_CPlas         : 0   : mg/L*hr, Area under curve of PFOS in plasma
AUC_CPlas_Fet     : 0   : mg/L*hr, Area under curve of PFOS in fetal plasma


$ODE
// #+ Concentrations in the tissues and in the venous palsma leaving each of the tissues (Unit: mg/L) 
// #+ Concentrations for mother; CX indicate the cocentration in the tissue (e.g., CL); CVX represent the concentration of chemical in tissue leaving tissue 
// #+ CPlas_free  : mg/L, Free PFOS concentration in the plasma
// #+ CPlas       : mg/L, Concentration of total PFOS in the plasma
// #+ CL          : mg/L, Concentration of PFOS in the liver compartment
// #+ CKb         : mg/L, Concentration of PFOS in venous plasma leaving kidney
// #+ CK          : mg/L, Concentration of PFOS in Kidney compartment
// #+ CM          : mg/L, Concentration of PFOS in the mammary gland compartment
// #+ CF          : mg/L, Concentration of PFOS in the fat compartment
// #+ CAm         : mg/L, Concentration of PFOS in the Amniotic fluid compartment
// #+ CRest       : mg/L, Concentration of PFOS in the rest of the body
// #+ CPTC        : mg/L, Concentration of PFOS in proximal tubule cells (PTC)
// #+ CFil        : mg/L, Concentration of PFOS in Filtrate (fil)
// #+ CPla        : mg/L, Concentration of PFOS in placenta
// #+ CVL         : mg/L, Concentration in the venous blood leaving liver
// #+ CVK         : mg/L, Concentration in the venous blood leaving kidney
// #+ CVM         : mg/L, Concentration in the venous blood leaving mammary gland
// #+ CVF         : mg/L, Concentration in the venous blood leaving fat
// #+ CVRest      : mg/L, Concentration in the venous blood leaving rest of body
// #+ CVPla       : mg/L, Concentration in the venous blood leaving placenta

double CPlas_free = APlas_free/VPlas;                             
double CPlas      = CPlas_free/Free;                                      
double CL         = AL/VL;                                             
double CKb        = AKb/VKb;                                          
double CK         = CKb*PK;                                            
double CM         = AM/VM_P;                                             
double CF         = AF/VF_P;                                            
double CAm        = AAm/(VAm + 1E-7);
double CRest      = ARest/VRest;                                    
double CPTC       = APTC/VPTC;                                       
double CFil       = AFil/VFil;                                       
double CPla       = APla/(VPla+ 1E-7);
double CVL        = CL/PL;                                            
double CVK        = CKb;                                              
double CVM        = CM/PM;                                            
double CVF        = CF/PF;                                            
double CVRest     = CRest/PRest;                                   
double CVPla      = CPla/PPla;

// #+ Concentrations for fetus; Fet represent the fetus
// #+ CPlas_Fet_free  : mg/L, Free PFOS concentration in the fetal plasma
// #+ CPlas_Fet       : mg/L, Total PFOS concentration in the fetal plasma
// #+ CRest_Fet       : mg/L, PFOS concentration in the fetal rest of body
// #+ CL_Fet          : mg/L, PFOS concentration in the fetal liver
// #+ CVL_Fet         : mg/L, PFOS concentration in the venous blood leaving fetal liver

double CPlas_Fet_free  = APlas_Fet_free/(VPlas_Fet + 1E-7);
double CPlas_Fet  = CPlas_Fet_free/Free_Fet;
double CRest_Fet  = ARest_Fet/(VRest_Fet + 1E-7);
double CVRest_Fet = ARest_Fet/((VRest_Fet + 1E-7)*PRest_Fet);
double CL_Fet     = AL_Fet/ (VL_Fet + 1E-7);
double CVL_Fet    = CL_Fet/ PL_Fet;

// #+ Equation for estimating the rate of compartment in mother PBPK
// #+ Virtural kidney sub-compartment 
// #+ RA_baso     : mg/h, Rate of basolateral transporters
// #+ RA_apical   : mg/h, Rate of apical transporter
// #+ Rdif        : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux    : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI         : mg/h, Rate of clearance(CL) to via glomerular filtration (GFR)
// #+ RPTC        : mg/h, Rate of change in PTC
// #+ RFil        : mg/h, Rate of change in Fil
// #+ RKb         : mg/h, Rate of change in Kidney-serum compartment
// #+ RST         : mg/h, Rate of change in stomach compartment
// #+ RSI         : mg/h, Rate of change in small intestines
// #+ RabsST      : mg/h, Rate of change of absorption in Stomach
// #+ RabsSI      : mg/h, Rate of change of absorption in small intestines
// #+ RL          : mg/h, Rate of change in liver compartment
// #+ RF          : mg/h, Rate of change in fat compartment
// #+ RM          : mg/h, Rate of change in mammary tissues
// #+ RRest       : mg/h, Rate of change in rest of body
// #+ RPla        : mg/h, Rate of change in placenta compartment
// #+ Rurine      : mg/h, Rate of change in urine
// #+ Rtrans_1    : mg/h, Rate of change in placenta tranfer from mother to fetus
// #+ Rtrans_2    : mg/h, Rate of change in placenta tranfer from fetus to mother
// #+ Rtrans_3    : mg/h, Rate of change in amniotic fluid tranfer from fetus to fluid
// #+ Rtrans_4    : mg/h, Rate of change in amniotic fluid tranfer from fluid to fetus
// #+ RAm         : mg/h, Rate of change in amniotic fluid compartment 
// #+ Rbile       : mg/h, Rate of change in bile compartment 
// #+ RL_Fet      : mg/h, Rate of change in fetal liver
// #+ RRest_Fet   : mg/h, Rate of change in fetal rest of body
// #+ RPlas_Fet_free : mg/h, Rate of free PFOS change in the fetal plasma

double RA_baso    = (Vmax_baso*CKb)/(Km_baso + CKb);                
double RA_apical  = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif       = Kdif*(CKb - CPTC);                               
double RAefflux   = Kefflux*APTC;                                
double RCI        = CPlas*GFR*Free;                              
double RPTC       = Rdif + RA_apical + RA_baso - RAefflux;           
double RFil       = RCI - RA_apical - AFil*Kurine;           
double RKb        = QK_P*(CPlas - CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RST        = -K0*AST - GE*AST;
double RSI        = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsST     = K0*AST;
double RabsSI     = Kabs*ASI;
double RL         = QL_P*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;   
double RF         = QF_P*(CPlas - CVF)*Free;                            
double RM         = QM_P*(CPlas - CVM)*Free;                            
double RRest      = QRest*(CPlas - CVRest)*Free;                         
double Rfeces     = Kbile*AL + Kunabs*ASI;                       
double Rurine     = Kurine*AFil;                                   
double Rtrans_1   = Ktrans_1*CVPla*Free;                          
double Rtrans_2   = Ktrans_2*CPlas_Fet*Free_Fet;                           
double Rtrans_3   = Ktrans_3*CVRest_Fet*Free_Fet;                              
double Rtrans_4   = Ktrans_4*CAm;
double RPla       = QPla*(CPlas - CVPla)*Free + Rtrans_2 - Rtrans_1;     
double RPlas_free = (QRest*CVRest*Free) + (QK_P*CVK*Free) + (QL_P*CVL*Free) + (QM_P*CVM*Free) + (QF_P*CVF*Free) + (QPla*CVPla*Free) - (QC*CPlas*Free) + RAefflux;  
double RAm        = Rtrans_3 - Rtrans_4;                               
double RL_Fet     = QL_Fet*(CPlas_Fet - CVL_Fet)*Free_Fet;
double RRest_Fet  = QRest_Fet*(CPlas_Fet - CVRest_Fet)*Free_Fet - Rtrans_3 + Rtrans_4;     
double RPlas_Fet  = (QRest_Fet*CVRest_Fet*Free_Fet) + (QL_Fet*CVL_Fet*Free_Fet) - (QC_Fet*CPlas_Fet*Free_Fet) + Rtrans_1 - Rtrans_2;
double Abile      = Kbile*AL;

// #+ ODE equation for mother compartment
dxdt_A_baso       = RA_baso;                                          
dxdt_A_apical     = RA_apical;                                       
dxdt_Adif         = Rdif;                                              
dxdt_Aefflux      = RAefflux;                                       
dxdt_APTC         = RPTC;                                           
dxdt_AFil         = RFil;                                            
dxdt_AKb          = RKb;                                                
dxdt_AST          = RST;
dxdt_AabsST       = RabsST;
dxdt_ASI          = RSI;
dxdt_AabsSI       = RabsSI;
dxdt_AL           = RL;                                                  
dxdt_AF           = RF;                                            
dxdt_AM           = RM;                                            
dxdt_ARest        = RRest;                                            
dxdt_Aurine       = Rurine;                                          
dxdt_Afeces       = Rfeces;                                          
dxdt_Atrans_1     = Rtrans_1;                                
dxdt_Atrans_2     = Rtrans_2;
dxdt_APla         = RPla;                                        
dxdt_APlas_free   = RPlas_free;                                  
dxdt_Atrans_3     = Rtrans_3;                                 
dxdt_Atrans_4     = Rtrans_4;                                  
dxdt_AAm          = RAm;                                         
dxdt_AL_Fet       = RL_Fet;
dxdt_ARest_Fet    = RRest_Fet;                                   
dxdt_APlas_Fet_free = RPlas_Fet;                         
dxdt_AUC_CPlas     = CPlas;
dxdt_AUC_CPlas_Fet = CPlas_Fet;

// #+ Virtural compartment for estmating input dose
dxdt_ADOSE        = 0;

// #+ Mass Balance check (Mother)
double ATissue    = AF + AM + ARest + APlas_free + AKb + AL + APla + AFil + APTC + AST + ASI;
double ALoss      = Aurine + Atrans_1 - Atrans_2 + Afeces;
double ATotal     = ATissue + ALoss;
double Mbal       = ADOSE - ATotal*KDOSE; 

// #+ Mass Balance check (Fetus)
double ATissueF   = APlas_Fet_free + ARest_Fet + AL_Fet;
double ALossF     = Atrans_2 + Atrans_3 - Atrans_4;
double ATotalF    = ATissueF + ALossF;
double DoseF      = Atrans_1;
double MbalF      = DoseF - ATotalF; 

// Embryo/fetus concentration
double CFtotal    = ATissueF/(VFet + 1.0e-7)/1000; 
            
$TABLE
capture Plasma    = CPlas;
capture Liver     = CL;
capture Kidney    = CK;
capture Placenta  = CPla;
capture CordB     = CPlas_Fet;
capture Fliver    = CL_Fet;
'
##///////////////////////////////////////////////////////////////////////////////////////////////////
LHumanPBPK.code <- '
$PROB
# Lactational PFOS PBPK model for human
- Author : Wei-Chun Chou
- Date   : Feb, 2020
- Strucutre: GI tract, Plasma, Liver, Fat, Mammary gland, Kidney, Filtrate, PTC, neonatal liver and neonatal compartments (Plasma, liver, fat, kidney, filtrate, PTC)  
- Note1  : Initial physiological parameters and optimized parameters is matchced with the values in non-pregnant PBPK models
- Note2  : Growth equation of physioloical parameters values was taken from Loccisano et al., 2013; yang et al., 2019 

$PARAM @annotated
// #+  Maternal parameters
BW0                 : 68.05     : kg,                  Body weight                                   ; Default value was the BW at the end of gestation excluding the fetus weight
Htc                 : 0.44      : Unitless             Hematocrit for human                          ; Value obtained Davies and Morris et al. (1993)
QCC                 : 16.4      : L/h/kg^0.75,         Cardiac blood output                          ; Value obtained from yoon et al. (2011)
QLC                 : 0.25      : Unitless,            Fraction blood flow to liver (%QCC)           ; Value obtained from Brown 1997; Fisher 2000
QKC                 : 0.175	    : Unitless,            Fraction blood flow to kidney (%QCC)          ; Value obtained from Brown 1997, Forsyth 1968
QMC                 : 0.027     : Unitless,            Fraction blood flow to Mammary gland (%QCC)   ; Value obtained from yoon et al., 2011
QFC                 : 0.052     : Unitless,            Fraction blood flow to Fat (%QCC)             ; Value obtained from Brown 1997; yoon et al., 2011
VLC                 : 0.026     : Unitless,            Fraction of liver tissue volume (%BW)         ; Value obtained from Brown 1997
VKC                 : 0.004     : Unitless,            Fraction of kidney tissue volume (%BW)        ; Value obtained from Brown 1997
VMC                 : 0.0062    : Unitless,            Fraction of Mammary gland volume (%BW)        ; Value obtained from yoon et al., 2011
VMilk               : 0.25      : L,                   Residual milk volume                          ; Value obtained from Loccisano et al., 2013; Gentry et al., 2003
VFilC               : 8.4e-4    : L/kg BW,             Fraction of filtrate (10% of Kidney volume)   ; Value obatedin from Worley et al., 2017
VPTCC               : 1.35e-4   : L/kg kidney,         Volume of proximal tubule cells               ; Value calculated from Hsu et al., 2014 (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)
protein             : 2.0e-6    : mg protein/PTCs,     Amount of protein in proximal tubule cells    ; Value obtaedin from Addis et al., 1936
PL                  : 2.03      : Unitless,            Liver-to-plasma partition coefficient         ; Value optimized from Chou and Lin, 2019
PK                  : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value optimized from Chou and Lin, 2019
PM                  : 0.16      : Unitless,            Mammary gland-to-plasma partition coefficient ; Value from Loccisano et al., 2013
PF                  : 0.13      : Unitless,            Fat-to-plasma partition coefficient           ; Value from Loccisano et al., 2013
PRest               : 0.20      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value from Loccisano et al., 2013
PMilkM              : 1.9       : Unitless,            Milk-to-mammary gland partition coefficient   ; Value from Loccisano et al., 2012; rat data
PMilkP              : 0.01      : Unitless,            Milk-to-plasma partition coefficient          ; Value from Loccisano et al., 2013 
PAMilkC             : 0.5       : Unitless,            Permeability area cross product (mammary to milk); L/h/kg
MW                  : 500.126   : g/mol,               PFOS molecular mass                           ; Value from Worley et al. (2017)
Free                : 0.014     : Unitless,            Free fraction of PFOS in plasma               ; Value optimized from Chou and Lin, 2019
KbileC              : 1.3e-4    : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
KurineC             : 0.096     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value optimized from Chou and Lin, 2019                        
K0C                 : 1.000     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in stomach         ; Value from Worley et al. (2017); Assumed to be same as PFOA
KabsC               : 2.120     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value from Worley et al. (2017); Assumed to be same as PFOA 
KunabsC             : 7.05e-5   : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces    ; Value from Worley et al. (2017); Assumed to be same as PFOA  
GFRC                : 27.28     : L/hr/kg kiney,       Glomerular filtration rate (female)           ; Value from Corley, 2005
GEC                 : 3.510     : 1/(h*BW^0.25),       Gastric emptying constant            ; Value from Yang et al., 2014
Vmax_baso_invitro   : 479       : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value optimized from Chou and Lin, 2019
Km_baso             : 20.1      : mg/L,                Km of basolateral transporter                  ; Value calculated from Worley et al. (2007)
Vmax_apical_invitro : 51803     : pmol/mg protein/min, Vmax of apical transporter                     ; Value optimized from Chou and Lin, 2019 
Km_apical           : 64.4      : mg/L,                Km of apical transporter                       ; Value optimized from Chou and Lin, 2019
RAFbaso             : 1         : Unitless             Relative activity factor for baso             ; Value calculated from Worley et al. (2017); Assumed to be same as PFOA 
RAFapi              : 0.001     : Unitless             Relative activity factor for api              ; Value optimized from Chou and Lin, 2019 
Kdif                : 0.001     : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value from Worley et al. (2017); Assumed to be same as PFOA
KeffluxC            : 0.150     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value optimized from Chou and Lin, 2019

// Neonate parameters
Vmax_baso_invitro_neo : 479     : pmol/mg Protein/min, Vmax of basolateral transporter               ; Value was assumed to be same with mother
Km_baso_neo         : 20.1      : mg/L,                Km of basolateral transporter                  ; Value was assumed to be same with mother
Vmax_apical_invitro_neo: 51803  : pmol/mg protein/min, Vmax of apical transporter                    ; Value was assumed to be same with mother
Km_apical_neo       : 64.4      : mg/L,                Km of apical transporter                       ; Value was assumed to be same with mother  
PL_neo              : 2.03      : Unitless,            Liver-to-plasma partition coefficient         ; Value was assumed to be same with mother  
PK_neo              : 1.26      : Unitless,            Kidney-to-plasma partition coefficient        ; Value was assumed to be same with mother 
PRest_neo           : 0.20      : Unitless,            Rest of body-to-plasma partition coefficient  ; Value was assumed to be same with mother 
Free_neo            : 0.014     : Unitless,            Free fraction of PFOS in plasma               ; Value was assumed to be same with mother 
KabsC_neo           : 2.120     : 1/(h*BW^-0.25),      Rate of absorption of PFOS in small intestines; Value was assumed to be same with mother 
KbileC_neo          : 1.3e-4    : 1/(h*BW^-0.25),      Biliary elimination rate                      ; Value was assumed to be same with mother 
KurineC_neo         : 0.001     : 1/(h*BW^-0.25),      Urinary elimination rate                      ; Value from Loccisano et al., 2013                    
Kdif_neo            : 0.001     : L/h,                 Diffusion rate from PTCs to kidney serum      ; Value was assumed to be same with mother 
KeffluxC_neo        : 0.150     : 1/(h*BW^-0.25),      Rate of clearance of PFOS from PTCs into blood; Value was assumed to be same with mother

// #+ Sensitivity anlaysis (SA) parametes for growth equations
SA_BW               : 1 : Unitless, SA for maternal body weight during lactation (BW)
SA_VPlasC           : 1 : Unitless, SA for maternal Fraction of volume of plasma (VPlasC)
SA_VFC              : 1 : Unitless, SA for maternal Fraction of volume of fat (VFC)
SA_KMilk            : 1 : Unitless, SA for maternal maternal milk production rate (KMilkC)
SA_VM               : 1 : Unitless, SA for maternal volume of mammary tissue (VM)
SA_GFR              : 1 : Unitless, SA for maternal Glomerular Filtration Rate (GFR)
SA_BW_neo           : 1 : Unitless, SA for neonatal body weight (BW_neo)
SA_Htc_neo          : 1 : Unitless, SA for neonatal Hematocrit (Htc_neo)
SA_VPlasC_neo       : 1 : Unitless, SA for neonatal Fraction of volume of plasma (VPlasC_neo)



$MAIN
// #+ Time varabiles: Gestational day and age
// Gestational day and age
double PND          = TIME/24;                                     
double WK           = PND/7;                                        
double Mon          = PND/30;
double Yr           = Mon/12;
double PMA          = (40 + WK);


// #+ Growth equations for maternal or fetal physiological parameters; Equation from Loccisano et al. (2013) and Yang et al. (2019) 
// #+ BW            : kg, lactating women body weight; equation fitting from Loccisano et al. (2013)
// #+ MK            : g,         Kidney weight                                 ; based on density of kidney = 1.0 g/mL
// #+ ML            : g,         Liver weight                                  ; based on density of liver = 1.05 g/mL, Overmeyer 1987


double BW           = SA_BW*(0.0014*pow(WK,2) - 0.1227*WK + BW0);
double VPlasC       = SA_VPlasC*((1E-4)*pow(WK,2) - 0.0024*WK + 0.0469);
double VPlasC_0     = 0.0469;
double VFC          = SA_VFC*((-7E-4)*WK + 0.3026);
double VFC_0        = 0.3026;
double VFC_F        = (-7E-4)*48 + 0.3026;
double KMilk_0      = (1E-9)*pow(PND, 3) - (1E-6)*pow(PND, 2) + 0.0002*PND + 0.0211;
double KMilk        = Mon <= 6? SA_KMilk*(KMilk_0): 0;
double VMT          = SA_VM*((-3E-6)*pow(PND, 2) + 0.0006*PND + 0.0059);
double GFR          = SA_GFR*(0.0259*pow(WK,2) - 0.4369*WK + 7.7972);

// #+ Tissue volumes
double VPlas        = VPlasC*BW;
double VL           = VLC*BW0;
double VK           = VKC*BW0;
double VF           = VFC*BW;
double VF_0         = VFC_0*BW;
double VM           = (1 + (VMT/VMC))*(VMC*BW0);
double VM_0         = (VMC + 0.05)*BW0;

// #+ Volumes of kidney and related compartmnet
double MK           = VKC*BW*1000;                                                                            
double ML           = VLC*BW*1.05*1000;                                       
double VKb          = VK*0.16;                            
double VFil         = VFilC*BW;                                      
double VPTC         = MK*VPTCC;  
double VRest        = (0.93*BW) - (VL + VK + VM + VF + VPlas + VMilk); 
double VBal         = (0.93*BW) - VRest - (VL + VK + VM + VF + VPlas + VMilk);  
                   
// #+ Maternal physiological parameters 
// #+ Blood flows
double QC_0         = QCC*pow(BW0,0.75)*(1 - Htc);
double QC_P         = (QC_0*(1 + (QFC*((0.4/VFC_F)-1)) + (QMC*((0.05/VMC)-1))));
double QF_0         = QFC*QC_0*(0.4/VFC_F);
double QF           = (QF_0*(VF/VF_0));
double QM_0         = QMC*QC_0*(0.05/VMC);
double QM           = (QM_0*(VM/VM_0));
double QL           = QLC*QC_0;
double QK           = QKC*QC_0;
double QC           = QC_P + (QF - QF_0) + (QM - QM_0);
double QRest        = QC - (QL + QK + QM + QF);
double QBal         = QC - (QRest + QL + QK + QM + QF);

// #+ Scaled rate constants 
// #+ Kabs          : 1/h, Rate of absorption of PFOS from small intestine to liver
// #+ Kunabs        : 1/h, Rate of absorption of PFOS from small intestine to liver
// #+ Kbile         : 1/h, Biliary elimination, liver to feces storage
// #+ Kurine        : 1/h, Urinary elimination
// #+ GE            : 1/h, Gasric emptying time
// #+ K0            : 1/h, Rate of uptake from the stomach into the liver
// #+ PTC           : n,   Proximial tubule cells
// #+ MPTC          : g,   mass of the proximal tubule cells (assuming density 1 kg/L)
// #+ Vmax_basoC    : mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3)
// #+ Vmax_apicalC  : mg/h/kg BW^0.75, Vmax of apical transporters in in vitro studies (Oatp1a1)
// #+ Vmax_baso     : mg/h,
// #+ Vmax_apical   : mg/h,
// #+ Kefflux       : 1/h, Efflux clearance rate from PTC to blood  

double Kabs         = KabsC*pow(BW,(-0.25));                    
double Kunabs       = KunabsC*pow(BW,(-0.25));                
double Kbile        = KbileC*pow(BW,(-0.25));                  
double Kurine       = KurineC*pow(BW,(-0.25));               
double GE           = GEC*pow(BW,(-0.25));
double K0           = K0C*pow(BW,(-0.25));
double PTC          = VKC*6e7*1000;                                       
double MPTC         = VPTC*1000;                                         
double Vmax_basoC   = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_baso    = Vmax_basoC*pow(BW,0.75);                          
double Vmax_apical  = Vmax_apicalC*pow(BW,0.75);                      
double Kefflux      = KeffluxC*pow(BW,(-0.25));                                   

// #+ Neonatal physiological parameters
// #+ BW_neo        : kg, Neonatal body weight; fitting data from Locisanno et al. (2013); 0 - 20.5 month
// #+ Htc_neo       : Unitness, Neonatal hematocrit; equation obtained from Yang et al. (2019) table 1
// #+ VPlasC_neo    : L, Neonatal volume of plasma (as fraction of BW_neo); fitting data from Yoon et al. (2011), table S5
// #+ Kabs_neo      : 1/h, Neonatal absorption rate of PFOS in small intestines
// #+ QC_neo        : L/h, Neonatal cardiac output; equation from Yang et al. (2019) table 2; <6 month QC_neo = 34.8 (Björkman, 2005)
// #+ QK_neo        : L/h, Neonatal Total Kidney blood flow; 
// #+ VLC_neo       : L, Neonatal volume of liver; equation from Yang et al. (2019) table 1
// #+ VK_neo        : L, Neonatal volume of kidney; equation from Yang et al. (2019) table 1
// #+ MK_neo        : kg, Neonatal kidney weight
// #+ ML_neo        : kg, Neonatal liver weight
// #+ VKb_neo       : L, Neonatal volume of kidney blood; fraction blood volume of kidney (0.16) from Brown, 1997
// #+ VFil_neo      : L, Neonatal filtrate volume; VFilC was assumed to be same as mother
// #+ VRest_neo     : L, Neonatal rest of body volume
// #+ VBal_neo      : L, Blance check for volume
// #+ QC_neo        : L/h, Neonatal cardiac output; equation from Yang et al. (2019) table 2
// #+ QK_neo        : L/h, Neonatal kidney blood flow; equation from Yang et al. (2019) table 2
// #+ QL_neo        : L/h, Neonatal liver blood flow; equation from Yang et al. (2019) table 2; Unit transfromed using (0.06) 


double BW_neo       = SA_BW_neo*((-0.0183)*pow(Mon, 2) + 0.7732*Mon + 3.7711);
double Htc_neo      = WK > 27 ? SA_Htc_neo*((0.744*WK + 24.656)/100): 0.37;
double VPlasC_neo   = SA_VPlasC_neo*((4E-6)*pow(Yr,3) - 0.0003*pow(Yr,2) + 0.0059*Yr + 0.05); 
double VPlas_neo    = VPlasC_neo*BW_neo;
double VL_neo       = (0.036*BW_neo*1000 + 17.851)/1000;
double VK_neo       = ((0.0034*BW_neo*1000 + 24.36) + (0.0036*BW_neo*1000 + 2.207))/1000; 
double MK_neo       = VKC*BW_neo*1000;
double ML_neo       = VLC*BW_neo*1000;
double VPTC_neo     = MK_neo*VPTCC;
double VKb_neo      = VK_neo*0.16;
double VFil_neo     = VFilC*BW_neo;
double VRest_neo    = (0.93*BW_neo) - (VL_neo + VPlas_neo + VK_neo + VPTC_neo + VFil_neo);
double VBal_neo     = (0.93*BW_neo) - (VL_neo + VPlas_neo + VK_neo + VRest_neo + VPTC_neo + VFil_neo);
double QC_neo       = (214.81*BW_neo + 76.57)*(0.06);
double QK_neo       = 0.221*VK_neo*1000 + 0.572;
double QL_neo       = 0.026*QC_neo/0.06;
double Kabs_neo     = KabsC_neo*pow(BW_neo,(-0.25));                    



// Neonatal scaled rat constant
double QRest_neo        = QC_neo - (QL_neo + QK_neo);
double QBal_neo         = QC_neo - (QL_neo + QK_neo + QRest_neo);

// Kinetics parameters of kidney for neonate
double Vmax_basoC_neo   = (Vmax_baso_invitro_neo*RAFbaso*PTC*protein*60*(MW/1e12)*1000);       
double Vmax_apicalC_neo = (Vmax_apical_invitro_neo*RAFapi*PTC*protein*60*(MW/1e12)*1000);  
double Vmax_baso_neo    = Vmax_basoC_neo*pow(BW_neo,0.75);                          
double Vmax_apical_neo  = Vmax_apicalC_neo*pow(BW_neo,0.75);                     
double Kefflux_neo      = KeffluxC_neo*pow(BW_neo,(-0.25));                         
double PAMilk           = PAMilkC*pow(BW_neo, 0.75);
double Kurine_neo       = KurineC_neo*pow(BW_neo, (-0.25));
double Kbile_neo        = KbileC_neo*pow(BW_neo, (-0.25)); 
double GFR_neo          = 0.00108*exp(0.1328*PMA);

// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE       = (TIME==0)?0:1;

$INIT @annotated
// #+ Set up the initial concentration; the initialconcentration was obtained from the steady-state concentration from non-pregnant model
// #+ Monther
ADOSE               : 0   : mg, Amount of input PFOS dose; virtual compartment
APTC                : 0   : mg, Amount of PFOS in the proximal tubule cells
A_baso              : 0   : mg, Amount of PFOS in the apical subcompartment 
A_apical            : 0   : mg, Amount of PFOS in the apical subcompartment 
Adif                : 0   : mg, Amount of PFOS diffused into the proximal tubule cells (PTC)
Aefflux             : 0   : mg, Amount of PFOS effluxed from PTC to blood
AFil                : 0   : mg, Amount of PFOS in the filtrate
AKb                 : 0   : mg, Amount of PFOS in the kidney blood 
AM                  : 0   : mg, Amount of PFOS in the mammary gland
AMilk               : 0   : mg, Amount of PFOS in the milk compartmnet
AF                  : 0   : mg, Amount of PFOS in fat compartment
ARest               : 0   : mg, Amount of PFOS in rest of body
ASI                 : 0   : mg, Amount of PFOS in the small intestine compartment
AST                 : 0   : mg, Amount of PFOS in the stomach compartment
AabsST              : 0   : mg, Amount of absorbed PFOS in the stomach compartment
AabsSI              : 0   : mg, Amount of absorbed PFOS in the small intestine compartment
AL                  : 0   : mg, Amount of PFOS in the liver compartmnet
Aurine              : 0   : mg, Amount of PFOS in urine
Afeces              : 0   : mg, Amount of PFOS in the feces  
APlas_free          : 0   : mg, Amount of PFOS in the plasma 
Atrans              : 0   : mg, Amount of PFOS trasfered from milk to neonate
AUC_CPlas           : 0   : mg/L*hr, Area under curve of PFOS in plasma

// #+  Fuetus
A_baso_neo          : 0   : mg, Amount of PFOS in the baso subcompartment 
A_apical_neo        : 0   : mg, Amount of PFOS in the apical subcompartment
Adif_neo            : 0   : mg, Amount of PFOS diffused into the proximal tubule cells (PTC)
Aefflux_neo         : 0   : mg, Amount of PFOS effluxed from PTC to blood
APTC_neo            : 0   : mg, Amount of PFOS in proximal tubule cells
AFil_neo            : 0   : mg, Amount of PFOS in filtrate
AKb_neo             : 0   : mg, Amount of PFOS in the kidney blood 
Aurine_neo          : 0   : mg, Amount of PFOS in urine
Afeces_neo          : 0   : mg, Amount of PFOS in the feces  
ARest_neo           : 0   : mg, Amount of PFOS in the rest of body compartment
AGI_neo             : 0   : mg, Amount of PFOS in the GI tract compartment (Neonate only)
AL_neo              : 0   : mg, Amount of PFOS in the liver compartment
APlas_free_neo      : 0   : mg, Amount of PFOS in the plasma compartment
AUC_CPlas_neo       : 0   : mg/L*hr, Area under curve of PFOS in neonatal plasma

$ODE
// Model equations for Monther
// #+ Concentrations in the tissues and in the venous palsma leaving each of the tissues (Unit: mg/L) 
// #+ CVL           : mg/L, Concentration of PFOS in plasma leaving liver
// #+ CVK           : mg/L, Concentration of PFOS in plasma leaving Kidney
// #+ CVM           : mg/L, Concentration of PFOS in plasma leaving mammary gland
// #+ CVF           : mg/L, Concentration of PFOS in plasma leaving fat
// #+ CVRest        : mg/L, Concentration of PFOS in plasma leaving rest of body
// #+ CPlas_free    : mg/L, Concentration of free PFOS in the plasma
// #+ CPlas         : mg/L, Concentration of total PFOS in the plasma
// #+ CKb           : mg/L, Concentration of PFOS in venous plasma leaving kidney
// #+ CPTC          : mg/L, Concentration of PFOS in proximal tubule cells (PTC)
// #+ CFil          : mg/L, Concentration of PFOS in Filtrate (fil)
// #+ CK            : mg/L, Concentration of PFOS in Kidney compartment
// #+ CM            : mg/L, Concentration of PFOS in the mammary gland compartment
// #+ CF            : mg/L, Concentration of PFOS in the fat compartment
// #+ CRest         : mg/L, Concentration of PFOS in the rest of the body
// #+ CL            : mg/L, Concentration of PFOS in the liver compartment

double CVL            = AL/(VL*PL);                                            
double CVK            = AKb/VKb;                                              
double CVM            = AM/(VM*PM);                                            
double CVF            = AF/(VF*PF);
double CVRest         = ARest/((VRest + 1E-7)*PRest);                                   
double CPlas_free     = APlas_free/VPlas;                             
double CPlas          = CPlas_free/Free;                                      
double CKb            = AKb/VKb;                                          
double CPTC           = APTC/VPTC;                                      
double CFil           = AFil/(VFil + 1E-7);                                     
double CMilk          = AMilk/VMilk;
double CK             = CVK*PK;                                           
double CM             = AM/VM;                                             
double CF             = AF/VF;                                             
double CRest          = ARest/VRest;                                   
double CL             = AL/VL;                                             

// Concentrations in the tissues and in the venous palsma leaving each of the tissues           
double CVL_neo        = AL_neo/(VL_neo*PL_neo);
double CVK_neo        = AKb_neo/(VKb_neo);
double CVRest_neo     = ARest_neo/((VRest_neo + 1E-7)*PRest_neo);
double CPlas_free_neo = APlas_free_neo/(VPlas_neo + 1E-7);
double CPlas_neo      = CPlas_free_neo/Free_neo;
double CKb_neo        = AKb_neo/VKb_neo;
double CPTC_neo       = APTC_neo/VPTC_neo;
double CFil_neo       = AFil_neo/VFil_neo;


// #+ Equation for estimating the rate of compartment in mother PBPK
// #+ RA_baso         : mg/h, Rate of basolateral transporters
// #+ RA_apical       : mg/h, Rate of apical transporter
// #+ Rdif            : mg/h, Rate of diffusion from into the PTC
// #+ RAefflux        : mg/h, Rate of efflux clearance rate from PTC to blood
// #+ RCI             : mg/h, Rate of clearance(CL) to via glomerular filtration (GFR)
// #+ RPTC            : mg/h, Rate of change in PTC
// #+ RFil            : mg/h, Rate of change in Fil
// #+ RKb             : mg/h, Rate of change in Kidney-serum compartment
// #+ RST             : mg/h, Rate of change in stomach compartment
// #+ RSI             : mg/h, Rate of change in small intestines
// #+ RabsST          : mg/h, Rate of change of absorption in stomach
// #+ RabsSI          : mg/h, Rate of change of absorption in small intestines
// #+ RL              : mg/h, Rate of change in liver compartment
// #+ RF              : mg/h, Rate of change in fat compartment
// #+ RM              : mg/h, Rate of change in mammary tissues
// #+ RRest           : mg/h, Rate of change in rest of body
// #+ RPla            : mg/h, Rate of change in placenta compartment
// #+ Rurine          : mg/h, Rate of change in urine
// #+ Rbile           : mg/h, Rate of change in bile compartment 
// #+ RL_Fet          : mg/h, Rate of change in fetal liver
// #+ RRest_Fet       : mg/h, Rate of change in fetal rest of body
// #+ RPlas_Fet_free  : mg/h, Rate of free PFOS change in the fetal plasma

double RA_baso             = (Vmax_baso*CKb)/(Km_baso + CKb);                
double RA_apical           = (Vmax_apical*CFil)/(Km_apical + CFil);        
double Rdif                = Kdif*(CKb - CPTC);                              
double RAefflux            = Kefflux*APTC;                               
double RCI                 = CPlas*GFR*Free;                            
double RPTC                = Rdif + RA_apical + RA_baso - RAefflux;                         
double RFil                = RCI - RA_apical - AFil*Kurine;           
double RKb                 = QK*(CPlas - CVK)*Free - CPlas*GFR*Free - Rdif - RA_baso;  
double RST                 = -K0*AST - GE*AST;
double RabsST              = K0*AST;
double RSI                 = GE*AST - Kabs*ASI - Kunabs*ASI;
double RabsSI              = Kabs*ASI;
double RL                  = QL*(CPlas - CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;  
double RF                  = QF*(CPlas - CVF)*Free;                            
double RM                  = QM*(CPlas - CVM)*Free - PAMilk*(CVM*Free - CMilk/PMilkM);                           
double Rtrans              = KMilk*CMilk;
double RMilk               = PAMilk*(CVM*Free - CMilk/PMilkM) - Rtrans;
double RRest               = QRest*(CPlas - CVRest)*Free;                         
double RPlas_free          = (QF*CVF*Free) + (QM*CVM*Free) + (QRest*CVRest*Free) + (QK*CVK*Free) + (QL*CVL*Free) - (QC*CPlas*Free) + RAefflux;  
double Rurine              = Kurine*AFil;                                   
double Rfeces              = Kbile*AL + Kunabs*ASI;
double Rbile               = Kbile*AL;

// #+ Equation for the rate of PFOS amount in the compartment of neonatal PBPK
double RA_baso_neo         = (Vmax_baso_neo*CKb_neo)/(Km_baso_neo + CKb_neo);                
double RA_apical_neo       = (Vmax_apical_neo*CFil_neo)/(Km_apical_neo + CFil_neo);        
double Rdif_neo            = Kdif_neo*(CKb_neo  - CPTC_neo);                               
double RAefflux_neo        = Kefflux_neo *APTC_neo ;                                
double RCI_neo             = CPlas_neo*GFR_neo*Free_neo;                             
double RPTC_neo            = Rdif_neo + RA_apical_neo + RA_baso_neo - RAefflux_neo;           
double RFil_neo            = RCI_neo - RA_apical_neo - AFil_neo*Kurine_neo;           
double RKb_neo             = QK_neo*(CPlas_neo - CVK_neo)*Free_neo - RCI_neo - Rdif_neo - RA_baso_neo;  
double RRest_neo           = QRest_neo*(CPlas_neo - CVRest_neo)*Free_neo;     
double RGI_neo             = Rtrans - Kabs_neo*AGI_neo;     
double RL_neo              = QL_neo*(CPlas_neo - CVL_neo)*Free_neo + Kabs_neo*AGI_neo - Kbile_neo*AL_neo*Free_neo;
double RPlas_free_neo      = (QRest_neo*CVRest_neo*Free_neo) + (QL_neo*CVL_neo*Free_neo) + (QK_neo*CVK_neo*Free_neo) - (QC_neo*CPlas_neo*Free_neo) + RAefflux_neo; ;
double Rurine_neo          = AFil_neo*Kurine_neo;
double Rfeces_neo          = Kbile_neo*AL_neo*Free_neo;

// #+ ODE equation for mother compartment
dxdt_A_baso                = RA_baso;                                           
dxdt_A_apical              = RA_apical;                                       
dxdt_Adif                  = Rdif;                                             
dxdt_Aefflux               = RAefflux;                                       
dxdt_APTC                  = RPTC;                                            
dxdt_AFil                  = RFil;                                            
dxdt_AKb                   = RKb;                                              
dxdt_AST                   = RST;
dxdt_AabsST                = RabsST;
dxdt_ASI                   = RSI;
dxdt_AabsSI                = RabsSI;
dxdt_AL                    = RL;                                                  
dxdt_AF                    = RF;                                           
dxdt_AM                    = RM;                                            
dxdt_Atrans                = Rtrans;
dxdt_AMilk                 = RMilk;
dxdt_ARest                 = RRest;                                            
dxdt_APlas_free            = RPlas_free;                                  
dxdt_Aurine                = Rurine;                                         
dxdt_Afeces                = Rfeces;                                          
dxdt_ADOSE                 = 0;
dxdt_AUC_CPlas             = CPlas;

// #+ ODE equation for neonatal compartment
dxdt_A_baso_neo            = RA_baso_neo ;                                           
dxdt_A_apical_neo          = RA_apical_neo ;                                       
dxdt_Adif_neo              = Rdif_neo ;                                              
dxdt_Aefflux_neo           = RAefflux_neo;                                       
dxdt_APTC_neo              = RPTC_neo;                                            
dxdt_AFil_neo              = RFil_neo;                                            
dxdt_AKb_neo               = RKb_neo;                                                
dxdt_AGI_neo               = RGI_neo;                                   
dxdt_AL_neo                = RL_neo;
dxdt_ARest_neo             = RRest_neo;                                   
dxdt_Aurine_neo            = Rurine_neo;
dxdt_APlas_free_neo        = RPlas_free_neo;                         
dxdt_Afeces_neo            = Rfeces_neo;  
dxdt_AUC_CPlas_neo         = CPlas_neo;

//## Mass Balance check (Monther)
double ATissue             = AF + AM + ARest + APlas_free + AKb + AL + AFil + APTC + AST + ASI;
double ALoss               = Aurine + Atrans + Afeces + AMilk;
double ATotal              = ATissue + ALoss;
double Mbal                = ADOSE - ATotal*KDOSE; 

//## Mass Balance check (Neonate)
double ATissue_neo         = ARest_neo + AKb_neo + AL_neo + APlas_free_neo + AGI_neo + AFil_neo + APTC_neo;
double ALoss_neo           = Aurine_neo + Afeces_neo;
double ATotal_neo          = ATissue_neo + ALoss_neo;
double MbalF               = Atrans - ATotal_neo; 

$TABLE
capture Plasma    = CPlas;
capture Liver     = CL;
capture Kidney    = CK;
capture CPneo     = CPlas_neo;
capture Milk      = CMilk;
'




