//Model 2; Voriconazole PBPK model for a typical male child
$PROB model2


$PARAM 
//References for the following parameters are:
//http://journals.sagepub.com/doi/abs/10.1016/S0146-6453%2803%2900002-2
//http://dmd.aspetjournals.org/content/38/1/25.long
//https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
Ka = 0.849 //absorption rate constant(/hr) 
fup = 0.42 //fraction of unbound drug in plasma
fumic = 0.711 //fraction of unbound drug in microsomes
WEIGHT = 19 //(kg)
MPPGL = 26 //pediatric mg microsomal protein per g liver (mg/g)
MPPGI = 0 //pediatric mg microsomal protein per g intestine (mg/g); 0 indicates no intestinal clearance
C_OUTPUT = 3.4 //cardiac output (l/min)
VmaxH = 120.5 //pediatric hepatic Vmax (pmol/min/mg)
VmaxG = 120.5 //pediatric intestinal Vmax (pmol/min/mg)
KmH = 11 //pediatric hepatic Km (uM)
KmG = 11 //pediatric intestinal Km (uM)
BP = 1 //blood to plasma ratio

//partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
//calculated using the script calcKp_PT.R
Kpad = 9.89  //adipose
Kpbo = 7.91  //bone
Kpbr = 7.35  //brain
Kpgu = 5.82  //gut
Kphe = 1.95  //heart
Kpki = 2.9  //kidneys
Kpli = 4.66  //liver
Kplu = 0.83  //lungs
Kpmu = 2.94  //muscle
Kpsp = 2.96  //spleen
Kpre = 4 //Rest of body; calculated as average of non adipose Kps


$CMT 
D ADIPOSE BRAIN GUT HEART BONE 
KIDNEY LIVER LUNG MUSCLE SPLEEN REST 
ART VEN


$MAIN
//Tissue volumes (L); http://journals.sagepub.com/doi/abs/10.1016/S0146-6453%2803%2900002-2
double Vad = 5.5; //adipose
double Vbo = 2.43; //bone
double Vbr = 1.31; //brain
double VguWall = 0.22; //gut wall
double VguLumen = 0.117; //gut lumen
double Vgu = VguWall + VguLumen;  //total gut
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

//fractions of tissue blood flows; http://journals.sagepub.com/doi/abs/10.1016/S0146-6453%2803%2900002-2
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


$ODE
//Calculation of tissue drug concentrations (mg/L)
double Cadipose = ADIPOSE/Vad;
double Cbone = BONE/Vbo;
double Cbrain = BRAIN/Vbr; 
double Cgut = GUT/VguWall; 
double Cheart = HEART/Vhe; 
double Ckidney = KIDNEY/Vki;
double Cliver = LIVER/Vli; 
double Clung = LUNG/Vlu; 
double Cmuscle = MUSCLE/Vmu;
double Cspleen = SPLEEN/Vsp;
double Crest = REST/Vre;
double Carterial = ART/Var;
double Cvenous = VEN/Vve;
double Cgut_D = D/Vgu;

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

//ODEs
dxdt_D = -Ka*D - CLintGut*Cgut_D; 
dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP)) - CLrenal*(Ckidneyfree/(Kpki/BP));
dxdt_GUT = Ka*D + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
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
  
  
  