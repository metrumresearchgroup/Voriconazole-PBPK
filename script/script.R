#This file contains the scripts used for simulations, plots and tables for the voriconazole PBPK manuscript
#Plots and tables are shown in the order of appearance in the manuscript
#First thing is to set working directory to "script" directory using the function: setwd()
#Note: All observed data and Zane and Thakker predictions were digitized from Zane and Thakker (ZT) paper. Zane and Thakker. Clin Pharmacokinet (2014) 53:1171–1182
#Note: All digitization was done using WebPlotDigitizer version 3.12 developed by Ankit Rohatgi https://automeris.io/WebPlotDigitizer

#loading libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(magrittr)
library(PKPDmisc)
library(gridExtra)
library(magick)
library(kableExtra)
library(mrgsolve)  #https://github.com/metrumresearchgroup/mrgsolve
source("funs.R")

#compile models
model1 <- mread_cache("model1", "../model")  # model1 (Adult physiology)
model2 <- mread_cache("model2", "../model")  # model2 (Pediatric physiology)


###################################################################################################
######################################### Chunk 1: Figure 3 #######################################
###################################################################################################
## This chunk plots figure 1 that compares model predictions using different calculation methods
# Source the different calculation methods
source("CalcKp_P&T.R")
source("CalcKp_R&R.R")
source("CalcKp_Berez.R")
source("CalcKp_Schmitt.R")
source("CalcKp_pksim.R")

# load observed data
load("../data/Fig3a_obs.Rda")  #load observed data (digitized) from fig 3a in the ZT paper
load("../data/Fig3a_ZT.Rda")  #load ZT redicitons (digitized) from fig 3a in the ZT paper

## drug-related parameters
type <- 3
logP <- 2.56
pKa <- 1.76
fup <- .42
BP <- 1

BW <- 73

## calculate the different Kps
dat <- read.csv("../data/tissue_comp_P&T.csv")
Kp_PT <- calcKp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type)
dat <- read.csv("../data/PKSim_tissue_comp_PT_Berez.csv")
Kp_Berez <- calcKp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type)
dat <- read.csv("../data/tissue_comp_R&R.csv")
Kp_RR <- calcKp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type)
#dat <- read.csv("../data/tissue_comp_Schmitt.csv")
dat <- read.csv("../data/PKSim_tissue_comp_Schmitt.csv")
Kp_Schmitt <- calcKp_Schmitt(logP=logP, pKa=pKa, fup=fup, type=type)
dat <- read.csv("../data/PKSim_tissue_comp_pksim.csv")
Kp_pksim <- calcKp_pksim(logP=logP, fup=fup)

## get all predictions
pred_PT <- model1 %>% 
  param(Kp_PT) %>% 
  ev(cmt = "VEN", amt = 4*BW, ii = 12, addl =  13, rate = 4*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_Berez <- model1 %>% 
  param(Kp_Berez) %>% 
  ev(cmt = "VEN", amt = 4*BW, ii = 12, addl =  13, rate = 4*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_RR <- model1 %>% 
  param(Kp_RR) %>% 
  ev(cmt = "VEN", amt = 4*BW, ii = 12, addl =  13, rate = 4*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_Schmitt <- model1 %>% 
  param(Kp_Schmitt) %>% 
  ev(cmt = "VEN", amt = 4*BW, ii = 12, addl =  13, rate = 4*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_pksim <- model1 %>% 
  param(Kp_pksim) %>% 
  ev(cmt = "VEN", amt = 4*BW, ii = 12, addl =  13, rate = 4*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

## plot
gp <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = pred_PT, mapping = aes(y=Cvenous, x=time, col="pred_PT"), lwd=1) + 
  geom_line(data = pred_Berez, mapping = aes(y=Cvenous, x=time, col="pred_Berez"), lwd=1, lty=2) +
  geom_line(data = pred_RR, mapping = aes(y=Cvenous, x=time, col="pred_RR"), lwd=1, lty=3) +
  geom_line(data = pred_Schmitt, mapping = aes(y=Cvenous, x=time, col="pred_Schmitt"), lwd=1, lty=4) +
  geom_line(data = pred_pksim, mapping = aes(y=Cvenous, x=time, col="pred_pksim"), lwd=1, lty=6) +
  scale_colour_manual(name='', values=c('ZT'='grey', 
                                        'pred_PT'='black', 
                                        'observed'='black',
                                        'pred_RR'='black',
                                        'pred_pksim'='black',
                                        'pred_Schmitt'='black',
                                        'pred_Berez'='black'), 
                      breaks=c("observed","ZT","pred_PT","pred_Berez","pred_RR","pred_Schmitt","pred_pksim"), 
                      labels=c("observed","ZT","PT","Berez","RR","Schmitt","PKSim")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1,2,3,4,6), shape=c(16, NA, NA, NA, NA, NA, NA)))) +
  ggtitle("Adult 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,2)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))

# plot and save
g <- grid.arrange(gp, ncol=1, nrow=2)
ggsave(file="../deliv/fig3.pdf", g, width=6, height=8)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 2: Table 1 #########################################
###################################################################################################
## This chunk tables out the prediction erros of different calculation methods
## pad observed data till 12 hours
df_temp <- obs %>%
  dplyr::filter(time > 7)
mod <- lm(obs~time, data=df_temp)
obs2 <- data.frame(time=12, obs=as.numeric(as.character(predict.lm(mod, newdata=data.frame(time=12)))), sd=NA)
obs <- bind_rows(obs, obs2)

## calculate AUCs for all methods
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,12))
auc_PT <- auc_partial(pred_PT$time, pred_PT$Cvenous, range=c(0,12))
auc_Berez <- auc_partial(pred_Berez$time, pred_Berez$Cvenous, range=c(0,12))
auc_RR <- auc_partial(pred_RR$time, pred_RR$Cvenous, range=c(0,12))
auc_Schmitt <- auc_partial(pred_Schmitt$time, pred_Schmitt$Cvenous, range=c(0,12))
auc_pksim <- auc_partial(pred_pksim$time, pred_pksim$Cvenous, range=c(0,12))

cmax_obs <- max(obs$obs)
cmax_PT <- max(pred_PT$Cvenous)
cmax_Berez <- max(pred_Berez$Cvenous)
cmax_RR <- max(pred_RR$Cvenous)
cmax_Schmitt <- max(pred_Schmitt$Cvenous)
cmax_pksim <- max(pred_pksim$Cvenous)

table1 <- data.frame(Method=c("Observed","PT","Berez","RR","Schmitt","PKSim"),
                     AUC=c(auc_obs, auc_PT, auc_Berez, auc_RR, auc_Schmitt, auc_pksim),
                     Cmax=c(cmax_obs, cmax_PT, cmax_Berez, cmax_RR, cmax_Schmitt, cmax_pksim))
table1 <- table1 %>%
  dplyr::mutate("RE_AUC (%)"=ifelse(Method=="Observed", NA, abs((first(AUC)-AUC)/first(AUC)*100)),
                "RE_Cmax (%)"=ifelse(Method=="Observed", NA, abs((first(Cmax)-Cmax)/first(Cmax)*100))) %>%
  dplyr::mutate("RE_total (%)"=`RE_AUC (%)` + `RE_Cmax (%)`) %>%
  dplyr::select(Method, AUC, `RE_AUC (%)`, Cmax, `RE_Cmax (%)`, `RE_total (%)`)
table1[,-1] <- round(table1[,-1], 2)

## view table
table1 %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover"))

## save as .csv
write.csv(table1, file="../deliv/table1.csv", row.names=F, quote=F)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 3: Figure 4 #########################################
###################################################################################################
## This chunk reproduces Figure 3 plots; all observed and ZT data were digitized and saved as .Rda 
## files whose numbers follow the ZT publication numbering
## Update models 1 and 2 with the RR Kps
model1 <- param(model1, Kp_RR)
model2 <- param(model2, Kp_RR)

## Figure 4a; Model 1 with 4 mg/kg IV infusion ##
load("../data/Fig3a_obs.Rda")  #load observed data (digitized) from fig 3a in the ZT paper
load("../data/Fig3a_ZT.Rda")  #load ZT redicitons (digitized) from fig 3a in the ZT paper
wt <- 73  #adult body weight
dose <- 4*wt  
rate <- 4*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(model1 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, rate=rate, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>% 
                            dplyr::filter(row_number() != 1)  

gp1 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 1")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("a   Adult 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,8) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp1


## Figure 3b; Model 2 with 4 mg/kg IV infusion ##
load("../data/Fig4a_obs.Rda")  #load observed data (digitized) from fig 3a in the ZT paper
load("../data/Fig4a_ZT.Rda")  #load ZT redicitons (digitized) from fig 3a in the ZT paper
wt <- 19  #pediatric body weight
dose <- 4*wt  
rate <- 3*wt
cmt <- "VEN"  #intravenous infusion

# simulate
sim <- as.data.frame(model2 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, rate=rate, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

gp2 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 2")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("b  Pediatric 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,8) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp2


#Figure 3c; Model 1 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  #load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

# simulate
sim <- as.data.frame(model1 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
                       dplyr::filter(row_number() != 1)  

gp3 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 1")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("c   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,3.5) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp3


#Figure 3d; Model 2 with 4 mg/kg PO 
load("../data/Fig4d_obs.Rda")  ##load observed data (digitized) from fig 4d in the ZT paper
load("../data/Fig4d_ZT.Rda")  ##load ZT predictions (digitized) from fig 4d in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

# simulate
sim <- as.data.frame(model2 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>% 
  dplyr::filter(row_number() != 1)  

gp4 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 2")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("d  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,3.5) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp4

#plot and save
g <- grid.arrange(gp1, gp2, gp3, gp4, ncol=2, nrow=3)
ggsave(file="../deliv/fig4.pdf", g, width=10, height=12)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 4: Figure 5 #########################################
###################################################################################################
### sensitivity analyses for absorption parameters: permeability, intestinal transit time and solubility ###
gp1 <- model1 %>%  
  knobs(S_lumen=c(0.39*500, 0.39*1000, 0.39*2000), amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>% 
  plot(Cvenous~., 
       main=list(label=expression(S["int"]~(mg/mL)), cex=0.75), 
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), 
       scale=list(cex=0.5), 
       lty=1:3, 
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.195","0.39","0.78"), cex=0.5)
       ))
gp1

gp2 <- model1 %>% 
  knobs(fperm=c(0.5, 1, 2), amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>% 
  plot(Cvenous~., 
       main=list(label=expression(paste("P"[eff], " x10",phantom()^{-4}," (cm/s)")), cex=0.75), 
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), 
       scale=list(cex=0.5), 
       lty=1:3, 
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.073","0.145","0.29"), cex=0.5)
       ))
gp2

gp3 <- model1 %>%
  knobs(ITT=c(3.32/2, 3.32, 3.32*2), amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>%
  plot(Cvenous~.,
       main=list(label=expression(ITT~(h)), cex=0.75),
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75),
       scale=list(cex=0.5),
       lty=1:3,
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1.66","3.32","6.64"), cex=0.5)
       ))
gp3

gp4 <- model1 %>%
  knobs(MPPGI=c(1.212/2, 1.212, 1.212*2), amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>%
  plot(Cvenous~.,
       main=list(label=expression(Cl["Gu"]~(mL/min/kg)), cex=0.75),
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75),
       scale=list(cex=0.5),
       lty=1:3,
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.035","0.07","0.14"), cex=0.5)
       ))
gp4

#plot and save
g <- grid.arrange(gp1, gp2, gp3, gp4, ncol=2, nrow=3)
ggsave(file="../deliv/fig5.pdf", g, width=10, height=12)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 5: Figure 7 #########################################
###################################################################################################
## This chunk reproduces Figure 5 plots
## calculate fperm (permeability factor)
VguWall <- 0.65
VguLumen <- 0.35
MW = 349.317  #(g/mol)
logP = 2.56  #log10 octanol oil:water partition coefficient; will be used as proxy for membrane affinity; preferably we will have phospholipid bilayer:water partition instead
S_lumen = 0.39*1000  #(mg/L) voriconazole intestinal lumen solubility https://www.ncbi.nlm.nih.gov/pubmed/24557773
L = 280  #(cm) small intestine length; from ICRP Publication 89
d = 2.5  #(cm) diameter of small intestine lumen
PF = 1.57  #3  //2.29  //average of 1.57 and 3; 1.57  //plicae circulare factor https://www.ncbi.nlm.nih.gov/pubmed/24694282
VF = 6.5  #villi factor
MF = 13  #microvilli factor
A = 7440  #this and the rest of parameters are constants in the permeability calculation equation https://www.ncbi.nlm.nih.gov/pubmed/15267240
B = 1e7
alpha = 0.6
beta = 4.395
fabs = 1  
fperm = 1  
SA_basal = pi*L*d*PF*VF*1e-4
MA = 10^logP
MW_eff = MW - (3*17)
Peff = fperm*A*((MW_eff^(-alpha-beta))*MA/(MW_eff^(-alpha) + B*(MW_eff^(-beta))*MA) * 1e-2 * 3600)
ka = fabs*Peff*SA_basal*1000/VguWall  #(h-1)
fperm <- 0.849/ka  #0.849 is the reported absorption rate constant

## Generate models 3,4,5,6,7 and 8
model3 <- param(model1, list(MPPGI=30.3/25))
model4 <- param(model2, list(MPPGI=26/25))
model5 <- param(model1, list(fperm=fperm))
model6 <- param(model2, list(fperm=fperm))
model7 <- param(model1, list(MPPGI=30.3/25, fperm=fperm))
model8 <- param(model2, list(MPPGI=26/25, fperm=fperm))

#Figure 5a; Model 3 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

sim <- as.data.frame(model1 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

sim_mppgi <- as.data.frame(model3 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  

sim_fperm <- as.data.frame(model5 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

sim_mppgi_fperm <- as.data.frame(model7 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

gp1 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1, lty=6) +
  geom_line(data = sim_mppgi, mapping = aes(y=Cvenous, x=time, col="sim_mppgi"), lwd=1, lty=2) + 
  geom_line(data = sim_fperm, mapping = aes(y=Cvenous, x=time, col="sim_fperm"), lwd=1, lty=3) +
  geom_line(data = sim_mppgi_fperm, mapping = aes(y=Cvenous, x=time, col="sim_mppgi_fperm"), lwd=1) +
  scale_colour_manual(name='', values=c('observed'='black',
                                        'ZT'='grey',
                                        'sim'='black',
                                        'sim_mppgi'='black', 
                                        'sim_fperm'='black',
                                        'sim_mppgi_fperm'='black'), 
                      breaks=c("observed","ZT","sim","sim_mppgi","sim_fperm","sim_mppgi_fperm"),
                      labels=c("observed","ZT","model 1","model 3","model 5","model 7")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,6,2,3,1), shape=c(16, NA,NA, NA,NA,NA)))) +
  ggtitle("a   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,3.5) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp1


#Figure 5b; Model 4 with 4 mg/kg IV infusion 
load("../data/Fig4d_obs.Rda")  ##load observed data (digitized) from fig 4d in the ZT paper
load("../data/Fig4e_ZT.Rda")  ##load ZT predictions (digitized) from fig 4d in the ZT paper
ZT_Gu <- ZT
load("../data/Fig4d_ZT.Rda")  ##load ZT predictions (digitized) from fig 4d in the ZT paper

wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

sim <- as.data.frame(model2 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

sim_mppgi <- as.data.frame(model4 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  

sim_fperm <- as.data.frame(model6 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

sim_mppgi_fperm <- as.data.frame(model8 %>% 
                                   ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                                   mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1) 

gp2 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = ZT_Gu, mapping = aes(x=time, y=ZT, col="ZT_Gu"), lwd=1, lty=2) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1, lty=6) +
  geom_line(data = sim_mppgi, mapping = aes(y=Cvenous, x=time, col="sim_mppgi"), lwd=1, lty=2) + 
  geom_line(data = sim_fperm, mapping = aes(y=Cvenous, x=time, col="sim_fperm"), lwd=1, lty=3) +
  geom_line(data = sim_mppgi_fperm, mapping = aes(y=Cvenous, x=time, col="sim_mppgi_fperm"), lwd=1) +
  scale_colour_manual(name='', values=c('observed'='black',
                                        'ZT'='grey',
                                        'ZT_Gu'='grey',
                                        'sim'='black',
                                        'sim_mppgi'='black', 
                                        'sim_fperm'='black',
                                        'sim_mppgi_fperm'='black'), 
                      breaks=c("observed","ZT","ZT_Gu","sim","sim_mppgi","sim_fperm","sim_mppgi_fperm"),
                      labels=c("observed","ZT",expression(ZT[Gu]),"model 2","model 4","model 6","model 8")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,2,6,2,3,1), shape=c(16, NA,NA,NA, NA,NA,NA)))) +
  ggtitle("b   Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  ylim(-0.5,3.5) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp2


#plot and save
g <- grid.arrange(gp1, gp2, ncol=2, nrow=3)
ggsave(file="../deliv/fig7.pdf", g, width=10, height=12)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 6: Table 2 #########################################
###################################################################################################
## This chunk tables out the prediction erros for IV predcitions

## Adult
## load data
load("../data/Fig3a_obs.Rda")
load("../data/Fig3a_ZT.Rda")

## pad observed data till 12 hours
df_temp <- obs %>%
  dplyr::filter(time > 7)
mod <- lm(obs~time, data=df_temp)
obs2 <- data.frame(time=12, obs=as.numeric(as.character(predict.lm(mod, newdata=data.frame(time=12)))), sd=NA)
obs <- bind_rows(obs, obs2)

## make prediction
BW <- 73
dose <- 4*BW
pred <- model1 %>% 
  ev(cmt = "VEN", amt = dose, ii = 12, addl =  13, rate = dose, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

## calculate AUCs for observed, ZT and prediction
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,12))
auc_ZT <- auc_partial(ZT$time, ZT$ZT, range=c(0,12))
auc_pred <- auc_partial(pred$time, pred$Cvenous, range=c(0,12))

## get cmax for observed and prediction
cmax_obs <- max(obs$obs)
cmax_ZT <- max(ZT$ZT)
cmax_pred <- max(pred$Cvenous)

df_temp <- data.frame(Type=c("Observed","ZT","Model 1"),
                      AUC=c(auc_obs, auc_ZT, auc_pred),
                      Cmax=c(cmax_obs, cmax_ZT, cmax_pred))
t1 <- df_temp %>%
  dplyr::mutate("RE_AUC (%)"=ifelse(Type=="Observed", NA, abs((first(AUC)-AUC)/first(AUC)*100)),
                "RE_Cmax (%)"=ifelse(Type=="Observed", NA, abs((first(Cmax)-Cmax)/first(Cmax)*100))) %>%
  dplyr::mutate("RE_total (%)"=`RE_AUC (%)` + `RE_Cmax (%)`) %>%
  dplyr::select(Type, AUC, `RE_AUC (%)`, Cmax, `RE_Cmax (%)`, `RE_total (%)`)

## Pediatric
## load data
load("../data/Fig4a_obs.Rda")
load("../data/Fig4a_ZT.Rda")

## make prediction
BW <- 19
dose <- 4*BW
pred <- model2 %>% 
  ev(cmt = "VEN", amt = dose, ii = 12, addl =  13, rate = 3*BW, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

## calculate AUCs for observed, ZT and prediction
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,12))
auc_ZT <- auc_partial(ZT$time, ZT$ZT, range=c(0,12))
auc_pred <- auc_partial(pred$time, pred$Cvenous, range=c(0,12))

## get cmax for observed and prediction
cmax_obs <- max(obs$obs)
cmax_ZT <- max(ZT$ZT)
cmax_pred <- max(pred$Cvenous)

df_temp <- data.frame(Type=c("Observed","ZT","Model 2"),
                      AUC=c(auc_obs, auc_ZT, auc_pred),
                      Cmax=c(cmax_obs, cmax_ZT, cmax_pred))
t2 <- df_temp %>%
  dplyr::mutate("RE_AUC (%)"=ifelse(Type=="Observed", NA, abs((first(AUC)-AUC)/first(AUC)*100)),
                "RE_Cmax (%)"=ifelse(Type=="Observed", NA, abs((first(Cmax)-Cmax)/first(Cmax)*100))) %>%
  dplyr::mutate("RE_total (%)"=`RE_AUC (%)` + `RE_Cmax (%)`) %>%
  dplyr::select(Type, AUC, `RE_AUC (%)`, Cmax, `RE_Cmax (%)`, `RE_total (%)`)

table2 <- bind_rows(t1,t2)
table2[,-1] <- round(table2[,-1], 2)

## view
table2 %>%
  kable() %>%
  kable_styling() %>%
  group_rows("Adult 4 mg/kg IV", 1, 3) %>%
  group_rows("Pediatric 4 mg/kg IV", 4, 6)

## save as .csv
write.csv(table2, file="../deliv/table2.csv", row.names=F, quote=F)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 7: Table 3 #########################################
###################################################################################################
## This chunk tables out the prediction erros for IV predcitions

## Adult
## load data
load("../data/Fig3b_obs.Rda")
load("../data/Fig3b_ZT.Rda")

## pad observed data till 12 hours
df_temp <- obs %>%
  dplyr::filter(time > 7)
mod <- lm(obs~time, data=df_temp)
obs2 <- data.frame(time=12, obs=as.numeric(as.character(predict.lm(mod, newdata=data.frame(time=12)))), sd=NA)
obs <- bind_rows(obs, obs2)

## make predictions
BW <- 73
dose <- 200
cmt <- "GUTLUMEN"
pred_init <- model1 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_mppgi <- model3 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_fperm <- model5 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_mppgi_fperm <- model7 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

## calculate AUCs for observed, ZT and prediction
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,12))
auc_ZT <- auc_partial(ZT$time, ZT$ZT, range=c(0,12))
auc_pred_init <- auc_partial(pred_init$time, pred_init$Cvenous, range=c(0,12))
auc_pred_mppgi <- auc_partial(pred_mppgi$time, pred_mppgi$Cvenous, range=c(0,12))
auc_pred_fperm <- auc_partial(pred_fperm$time, pred_fperm$Cvenous, range=c(0,12))
auc_pred_mppgi_fperm <- auc_partial(pred_mppgi_fperm$time, pred_mppgi_fperm$Cvenous, range=c(0,12))

## get cmax for observed and prediction
cmax_obs <- max(obs$obs)
cmax_ZT <- max(ZT$ZT)
cmax_pred_init <- max(pred_init$Cvenous)
cmax_pred_mppgi <- max(pred_mppgi$Cvenous)
cmax_pred_fperm <- max(pred_fperm$Cvenous)
cmax_pred_mppgi_fperm <- max(pred_mppgi_fperm$Cvenous)

df_temp <- data.frame(Type=c("Observed","ZT","Model 1","Model 3","Model 5","Model 7"),
                      AUC=c(auc_obs, auc_ZT, auc_pred_init, auc_pred_mppgi, auc_pred_fperm, auc_pred_mppgi_fperm),
                      Cmax=c(cmax_obs, cmax_ZT, cmax_pred_init, cmax_pred_mppgi, cmax_pred_fperm, cmax_pred_mppgi_fperm))
t1 <- df_temp %>%
  dplyr::mutate("RE_AUC (%)"=ifelse(Type=="Observed", NA, abs((first(AUC)-AUC)/first(AUC)*100)),
                "RE_Cmax (%)"=ifelse(Type=="Observed", NA, abs((first(Cmax)-Cmax)/first(Cmax)*100))) %>%
  dplyr::mutate("RE_total (%)"=`RE_AUC (%)` + `RE_Cmax (%)`) %>%
  dplyr::select(Type, AUC, `RE_AUC (%)`, Cmax, `RE_Cmax (%)`, `RE_total (%)`)


## Pediatric
## load data
load("../data/Fig4d_obs.Rda")
load("../data/Fig4e_ZT.Rda")

## make prediction
BW <- 19
dose <- 4*BW
cmt <- "GUTLUMEN"
pred_init <- model2 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_mppgi <- model4 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_fperm <- model6 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

pred_mppgi_fperm <- model8 %>% 
  ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss = 1) %>% 
  mrgsim(delta = 0.1, end = 12) %>%
  as.data.frame() %>%
  dplyr::filter(row_number() > 1)

## calculate AUCs for observed, ZT and prediction
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,12))
auc_ZT <- auc_partial(ZT$time, ZT$ZT, range=c(0,12))
auc_pred_init <- auc_partial(pred_init$time, pred_init$Cvenous, range=c(0,12))
auc_pred_mppgi <- auc_partial(pred_mppgi$time, pred_mppgi$Cvenous, range=c(0,12))
auc_pred_fperm <- auc_partial(pred_fperm$time, pred_fperm$Cvenous, range=c(0,12))
auc_pred_mppgi_fperm <- auc_partial(pred_mppgi_fperm$time, pred_mppgi_fperm$Cvenous, range=c(0,12))

## get cmax for observed and prediction
cmax_obs <- max(obs$obs)
cmax_ZT <- max(ZT$ZT)
cmax_pred_init <- max(pred_init$Cvenous)
cmax_pred_mppgi <- max(pred_mppgi$Cvenous)
cmax_pred_fperm <- max(pred_fperm$Cvenous)
cmax_pred_mppgi_fperm <- max(pred_mppgi_fperm$Cvenous)

df_temp <- data.frame(Type=c("Observed","ZT_Gu","Model 2","Model 4","Model 6","Model 8"),
                      AUC=c(auc_obs, auc_ZT, auc_pred_init, auc_pred_mppgi, auc_pred_fperm, auc_pred_mppgi_fperm),
                      Cmax=c(cmax_obs, cmax_ZT, cmax_pred_init, cmax_pred_mppgi, cmax_pred_fperm, cmax_pred_mppgi_fperm))
t2 <- df_temp %>%
  dplyr::mutate("RE_AUC (%)"=ifelse(Type=="Observed", NA, abs((first(AUC)-AUC)/first(AUC)*100)),
                "RE_Cmax (%)"=ifelse(Type=="Observed", NA, abs((first(Cmax)-Cmax)/first(Cmax)*100))) %>%
  dplyr::mutate("RE_total (%)"=`RE_AUC (%)` + `RE_Cmax (%)`) %>%
  dplyr::select(Type, AUC, `RE_AUC (%)`, Cmax, `RE_Cmax (%)`, `RE_total (%)`)

table3 <- bind_rows(t1,t2)
table3[,-1] <- round(table3[,-1], 2)

## view
table3 %>%
  kable() %>%
  kable_styling() %>%
  group_rows("Adult 200 mg PO", 1, 6) %>%
  group_rows("Pediatric 4 mg/kg PO", 7, 12)

## save as .csv
write.csv(table3, file="../deliv/table3.csv", row.names=F, quote=F)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 8: Figure 8 #########################################
###################################################################################################
# ## This chunk reproduces Figure S5 plots with the simulations from IIV and uncertainty of CL_Gu and Peff
# ## Figure S5a
# 
# ## Adults
# load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
# 
# ## sample Cl_Gu and Peff by sampling from MPPGI and fperm
# set.seed(23142)
# nSampl <- 1000
# mppgiSampl <- exp(rnorm(nSampl, mean=log(30.3/25), sd=0.33))  #lognormal; CV = 33%   
# fpermSampl <- exp(rnorm(nSampl, mean=log(fperm), sd=0.4))  #lognormal; CV = 40%
# lpars <- purrr::map2(mppgiSampl, fpermSampl, .f=function(x,y){list(x,y)})  #make a list of parameter samples
# lpars <- lapply(lpars, function(x) setNames(x, c("MPPGI","fperm")))  #set parameter names
# 
# #simulate
# wt <- 73
# dose <- 200
# cmt <- "GUTLUMEN"
# 
# sim <- function(pars, mod){
#   mod %>%
#     param(pars) %>%
#     ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss=1) %>% 
#     mrgsim(end = 12, delta = 0.1) %>%
#     dplyr::filter(row_number() > 1) %>%
#     dplyr::select(time, Cvenous)
# }
# 
# sims <- lapply(lpars, sim, mod=model1) %>% 
#   bind_rows %>%
#   group_by(time) %>%
#   dplyr::mutate(avg = mean(Cvenous),
#                 low = quantile(Cvenous, probs=0.025), 
#                 high=quantile(Cvenous, probs=0.975))
# 
# 
# # plot
# gp1 <- ggplot() + 
#   geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
#   geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
#   geom_line(data = sims, mapping = aes(y=avg, x=time, col="avg"), lwd=1) + 
#   geom_line(data = sims, mapping = aes(y=low, x=time, col="low"), lwd=1, lty=2) +
#   geom_line(data = sims, mapping = aes(y=high, x=time, col="high"), lwd=1, lty=2) +
#   scale_colour_manual(name='', values=c('observed'='black', 
#                                         'avg'='black', 
#                                         'low'='black', 
#                                         'high'='black'), 
#                       breaks=c("observed","avg","low","high"),
#                       labels=c("observed","mean","","")) +
#   guides(colour = guide_legend(override.aes = list(linetype=c(0,1,0,0), shape=c(16, NA, NA, NA)))) +
#   ggtitle("a   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
#   scale_y_continuous(breaks = seq(0,100,1)) +
#   #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
#   ylim(-0.5,3.5) +
#   scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
#   theme_bw() + 
#   theme(axis.title=element_text(size=12)) +
#   theme(axis.text=element_text(size=10)) + 
#   theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
#   theme(plot.title=element_text(size=15, face="bold")) +
#   theme(legend.text=element_text(size=10))
# gp1
# 
# 
# ## Pediatrics
# load("../data/Fig4d_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
# 
# ## sample Cl_Gu and Peff by sampling from MPPGI and fperm
# set.seed(23142)
# nSampl <- 1000
# mppgiSampl <- exp(rnorm(nSampl, mean=log(26/25), sd=0.33))  #lognormal; CV = 33%   
# fpermSampl <- exp(rnorm(nSampl, mean=log(fperm), sd=0.4))  #lognormal; CV = 40%
# lpars <- purrr::map2(mppgiSampl, fpermSampl, .f=function(x,y){list(x,y)})  #make a list of parameter samples
# lpars <- lapply(lpars, function(x) setNames(x, c("MPPGI","fperm")))  #set parameter names
# 
# #simulate
# wt <- 19
# dose <- 4*wt
# cmt <- "GUTLUMEN"
# 
# sim <- function(pars, mod){
#   mod %>%
#     param(pars) %>%
#     ev(cmt = cmt, amt = dose, ii = 12, addl =  13, ss=1) %>% 
#     mrgsim(end = 12, delta = 0.1) %>%
#     dplyr::filter(row_number() > 1) %>%
#     dplyr::select(time, Cvenous)
# }
# 
# sims <- lapply(lpars, sim, mod=model2) %>% 
#   bind_rows %>%
#   group_by(time) %>%
#   dplyr::mutate(avg = mean(Cvenous),
#                 low = quantile(Cvenous, probs=0.025), 
#                 high=quantile(Cvenous, probs=0.975))
# 
# 
# gp2 <- ggplot() + 
#   geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
#   geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
#   geom_line(data = sims, mapping = aes(y=avg, x=time, col="avg"), lwd=1) + 
#   geom_line(data = sims, mapping = aes(y=low, x=time, col="low"), lwd=1, lty=2) +
#   geom_line(data = sims, mapping = aes(y=high, x=time, col="high"), lwd=1, lty=2) +
#   scale_colour_manual(name='', values=c('observed'='black', 
#                                         'avg'='black', 
#                                         'low'='black', 
#                                         'high'='black'), 
#                       breaks=c("observed","avg","low","high"),
#                       labels=c("observed","mean","","")) +
#   guides(colour = guide_legend(override.aes = list(linetype=c(0,1,0,0), shape=c(16, NA, NA, NA)))) +
#   ggtitle("b   Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
#   scale_y_continuous(breaks = seq(0,100,1)) +
#   #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
#   ylim(-0.5,3.5) +
#   scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
#   theme_bw() + 
#   theme(axis.title=element_text(size=12)) +
#   theme(axis.text=element_text(size=10)) + 
#   theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
#   theme(plot.title=element_text(size=15, face="bold")) +
#   theme(legend.text=element_text(size=10))
# gp2
# 
# #plot and save
# g <- grid.arrange(gp1, gp2, ncol=2, nrow=3)
# ggsave(file="../deliv/fig8.pdf", g, width=10, height=12)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 9: Figure 9 #########################################
###################################################################################################
## This chunk reproduces Figure 9 which compares exposures of adult standard oral dose vs several
## pediatric doses
## Figure 9
pred_adult <- run_sim(model7, cmt="GUTLUMEN", dose=200)
pred_ped_4 <- run_sim(model8, cmt="GUTLUMEN", dose=4*19)
pred_ped_5 <- run_sim(model8, cmt="GUTLUMEN", dose=5*19)
pred_ped_6 <- run_sim(model8, cmt="GUTLUMEN", dose=6*19)
pred_ped_7 <- run_sim(model8, cmt="GUTLUMEN", dose=7*19)
pred_ped_8 <- run_sim(model8, cmt="GUTLUMEN", dose=8*19)
  
## calculate AUCs for observed, ZT and prediction
auc_pred_adult <- auc_partial(pred_adult$time, pred_adult$Cvenous, range=c(0,12))
auc_pred_ped_4 <- auc_partial(pred_ped_4$time, pred_ped_4$Cvenous, range=c(0,12))
auc_pred_ped_5 <- auc_partial(pred_ped_5$time, pred_ped_5$Cvenous, range=c(0,12))
auc_pred_ped_6 <- auc_partial(pred_ped_6$time, pred_ped_6$Cvenous, range=c(0,12))
auc_pred_ped_7 <- auc_partial(pred_ped_7$time, pred_ped_7$Cvenous, range=c(0,12))
auc_pred_ped_8 <- auc_partial(pred_ped_8$time, pred_ped_8$Cvenous, range=c(0,12))

## setup dataframe
df_temp <- data.frame(Dose = c("Adult 200 mg", 
                               "Pediatric 4 mg/kg", 
                               "Pediatric 5 mg/kg", 
                               "Pediatric 6 mg/kg", 
                               "Pediatric 7 mg/kg",
                               "Pediatric 8 mg/kg"),
                      AUC = c(auc_pred_adult,
                              auc_pred_ped_4,
                              auc_pred_ped_5,
                              auc_pred_ped_6,
                              auc_pred_ped_7,
                              auc_pred_ped_8))

#plot and save
gp <- ggplot(data=df_temp, aes(Dose, AUC)) +
  geom_bar(stat="identity") +
  labs(x="", y=expression(paste("AUC"["0-12,ss"], " (mg.h/L)"))) +
  geom_hline(yintercept=auc_pred_adult, lty=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
gp

g <- grid.arrange(gp, ncol=1, nrow=2)
ggsave(file="../deliv/fig9.pdf", g, width=6, height=8)

###################################################################################################
###################################################################################################
