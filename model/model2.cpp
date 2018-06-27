//Voriconazole PBPK model for a typical male child


$PROB voriPBPK_Pediatric


$PARAM 
fup = 0.42 //fraction of unbound drug in plasma
fumic = 0.711 //fraction of unbound drug in microsomes
WEIGHT = 19 //(kg)
MPPGL = 26 //pediatric mg microsomal protein per g liver (mg/g); from Yanni, Souzan et al. (2009)
MPPGI = 26 //pediatric mg microsomal protein per g intestine (mg/g);
C_OUTPUT = 3.4 //cardiac output (l/min); from ICRP Publication 89
VmaxH = 120.5 //pediatric hepatic Vmax (pmol/min/mg); from Yanni, Souzan et al. (2009)
VmaxG = 120.5 //pediatric intestinal Vmax (pmol/min/mg)
KmH = 11 //pediatric hepatic Km (uM); from Yanni, Souzan et al. (2009)
KmG = 11 //pediatric intestinal Km (uM)
BP = 1 //blood to plasma ratio; initial estimate

//partition coefficients estimated by Poulin and Theil method (2001)
Kpad = 9.89, Kpbo = 7.91, Kpbr = 7.35, Kpgu = 5.82, Kphe = 1.95, Kpki = 2.9, Kpli = 4.66
Kplu = 0.83, Kpmu = 2.94, Kpsp = 2.96, Kpre = 4 //calculated as average of non adipose Kps

//absorption model parameters
MW = 349.317  //(g/mol)
logP = 2.56  //log10 octanol oil:water partition coefficient; will be used as proxy for membrane affinity; preferably we will have phospholipid bilayer:water partition instead
S_lumen = 0.39*1000  //(mg/L) voriconazole intestinal lumen solubility https://www.ncbi.nlm.nih.gov/pubmed/24557773
L = 170  //(cm) small intestine length; from ICRP Publication 89
d = 2.5  //(cm) diameter of small intestine lumen
PF = 1.57  //3  //2.29  //3  //1.57  //plicae circulare factor https://www.ncbi.nlm.nih.gov/pubmed/24694282
VF = 6.5  //villi factor
MF = 13  //microvilli factor
ITT = 3.32  //(h) small intestine transit time; https://www.ncbi.nlm.nih.gov/pubmed/25986421
A = 7440  //this and the rest of parameters are constants in the permeability calculation equation https://www.ncbi.nlm.nih.gov/pubmed/15267240
B = 1e7
alpha = 0.6
beta = 4.395
fabs = 1  //absorption factor to manipulate ka
fdis = 1  //disappearance from gut lumen factor to manipulate kd
fperm = 1  //permeability factor to manipulate Pm


$CMT 
GUTLUMEN GUTWALL GUT ADIPOSE BRAIN HEART BONE 
KIDNEY LIVER LUNG MUSCLE SPLEEN REST 
ART VEN


$MAIN

//Tissue volumes (L); from ICRP Publication 89
double Vad = 5.5; //adipose
double Vbo = 2.43; //bone
double Vbr = 1.31; //brain
double VguWall = 0.22; //just small intestines; 0.4; //gut wall
double VguLumen = 0.117; //just small intestines; 0.3; //gut lumen
double Vgu = VguWall + VguLumen;
double Vhe = 0.085; //heart
double Vki = 0.11; //kidneys
double Vli = 0.467; //liver
double Vlu = 0.125; //lungs
double Vmu = 5.6; //muscle
double Vsp = 0.05; //spleen
double Vbl = 1.5; //total blood
double Vve = 0.705*Vbl; //venous blood
double Var = 0.295*Vbl; //arterial blood
double Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+VguWall+Vbl); //volume of rest of the body compartment

//fractions of tissue blood flows; from ICRP Publication 89
double FQad = 0.05;
double FQbo = 0.05;  
double FQbr = 0.12;
double FQgu = 0.16;
double FQhe = 0.04;
double FQki = 0.19;
double FQli = 0.255;
double FQmu = 0.17;
double FQsp = 0.03;

//computing the blood flows for each tissue
double CO = C_OUTPUT*60; //scaled cardiac output (L/hr)
double Qad = FQad*CO;
double Qbo = FQbo*CO;
double Qbr = FQbr*CO;
double Qgu = FQgu*CO;
double Qhe = FQhe*CO;
double Qki = FQki*CO;
double Qli = FQli*CO;
double Qmu = FQmu*CO;
double Qsp = FQsp*CO;
double Qha = Qli - (Qgu + Qsp); //hepatic artery 
double Qtot = Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr;
double Qre = CO - Qtot;
double Qlu = CO;

//absorption model parameters derivation
double SA_abs = M_PI*L*d*PF*VF*MF*1e-4;  //(m^2) mucosal absorption surface area https://www.ncbi.nlm.nih.gov/pubmed/24694282
double SA_basal = M_PI*L*d*PF*VF*1e-4;  //(m^2) basal membrane surface area https://www.ncbi.nlm.nih.gov/pubmed/24694282
double MA = pow(10,logP);  //membrane affinity
double MW_eff = MW - (3*17);  //effective molecular weight; voriconazole has 3 F atoms so we subtract 17 mass units per atom https://www.ncbi.nlm.nih.gov/pubmed/15267240
double Peff = fperm*A*((pow(MW_eff,(-alpha-beta))*MA)/(pow(MW_eff,(-alpha)) + B*pow(MW_eff,(-beta))*MA) * 1e-2 * 3600);  //(m/h) intestinal permeability 
double kd = fdis*Peff*SA_abs*1000/VguLumen;  //(h-1) rate constant for drug disappearing from lumen and into enterocytes
double ka = fabs*Peff*SA_basal*1000/VguWall; //(h-1) rate constant for drug absorption from enterocytes to gut circulation
double kt = 1/ITT;  //(h-1) intestinal transit rate constatnt

$ODE
//Calculation of tissue drug concentrations (mg/L)
double Cadipose = ADIPOSE/Vad;
double Cbone = BONE/Vbo;
double Cbrain = BRAIN/Vbr; 
double Cheart = HEART/Vhe; 
double Ckidney = KIDNEY/Vki;
double Cliver = LIVER/Vli; 
double Clung = LUNG/Vlu; 
double Cmuscle = MUSCLE/Vmu;
double Cspleen = SPLEEN/Vsp;
double Crest = REST/Vre;
double Carterial = ART/Var;
double Cvenous = VEN/Vve;
double CgutLumen = GUTLUMEN/VguLumen;
double CgutWall = GUTWALL/VguWall;
double Cgut = GUT/VguWall;

//Free Concentration Calculations
double Cliverfree = Cliver*fup; 
double Ckidneyfree = Ckidney*fup;

//scaling factor for intrinsic hepatic and intestinal clearances
double scale_factor_H = MPPGL*Vli*1000; //hepatic (mg)
double scale_factor_G = MPPGI*VguWall*1000; //intestinal (mg)                                         

//intrinsic hepatic and intestinal clearances calculation
double CLintHep = (VmaxH/KmH)*scale_factor_H*60*pow(10,-6); //(L/hr)
double CLintGut = (VmaxG/KmG)*scale_factor_G*60*pow(10,-6); //(L/hr)
CLintHep = CLintHep/fumic; 
CLintGut = CLintGut/fumic;

//renal clearance
double CLrenal = 0.096; //(L/hr); from Zane and Thakker (2014)

//accounting for solubility
double f = 1;
if(CgutLumen > S_lumen){
  f = 0;
}

//ODEs
dxdt_GUTLUMEN = - kd*VguLumen*(f*CgutLumen + (1-f)*S_lumen) - kt*GUTLUMEN;
dxdt_GUTWALL = kd*VguLumen*(f*CgutLumen + (1-f)*S_lumen) - ka*GUTWALL - CLintGut*CgutWall;
dxdt_GUT = ka*GUTWALL + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP)) - CLrenal*(Ckidneyfree/(Kpki/BP));
dxdt_LIVER = Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - 
  CLintHep*(Cliverfree/(Kpli/BP)); 
dxdt_LUNG = Qlu*(Cvenous - Clung/(Kplu/BP));
dxdt_MUSCLE = Qmu*(Carterial - Cmuscle/(Kpmu/BP));
dxdt_SPLEEN = Qsp*(Carterial - Cspleen/(Kpsp/BP));
dxdt_BONE = Qbo*(Carterial - Cbone/(Kpbo/BP));
dxdt_REST = Qre*(Carterial - Crest/(Kpre/BP));
dxdt_VEN = Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
  Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
  Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qre*(Crest/(Kpre/BP)) - Qlu*Cvenous;
dxdt_ART = Qlu*(Clung/(Kplu/BP) - Carterial);


$CAPTURE Cvenous
  
  
  