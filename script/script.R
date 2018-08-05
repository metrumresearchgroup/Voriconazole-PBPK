#This file contains the scripts used for simulations and plots for the voriconazole PBPK manuscript
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
library(kableExtra)
library(mrgsolve)  #https://github.com/metrumresearchgroup/mrgsolve

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
  theme_bw() + theme(axis.title=element_text(size=18)) +
  theme(axis.text=element_text(size=15)) +
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=20, face="bold")) +
  theme(legend.text=element_text(size=18))

# plot and save
g <- grid.arrange(gp, ncol=1, nrow=2)
ggsave(file="../deliv/fig3.pdf", g, width=6, height=8)

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 2: Table 1 #########################################
###################################################################################################
## This chink tables out the prediction erros of different calculation methods
auc_obs <- auc_partial(obs$time, obs$obs, range=c(0,10))
auc_PT <- auc_partial(pred_PT$time, pred_PT$Cvenous, range=c(0,10))
auc_Berez <- auc_partial(pred_Berez$time, pred_Berez$Cvenous, range=c(0,10))
auc_RR <- auc_partial(pred_RR$time, pred_RR$Cvenous, range=c(0,10))
auc_Schmitt <- auc_partial(pred_Schmitt$time, pred_Schmitt$Cvenous, range=c(0,10))
auc_pksim <- auc_partial(pred_pksim$time, pred_pksim$Cvenous, range=c(0,10))

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
  dplyr::mutate("AUC_error(%)"=ifelse(Method=="Observed", NA, round(abs((first(AUC)-AUC)/first(AUC)*100), digits=2)),
                "Cmax_error(%)"=ifelse(Method=="Observed", NA, round(abs((first(Cmax)-Cmax)/first(Cmax)*100), digits=2))) %>%
  dplyr::mutate("Total_error(%)"=`AUC_error(%)` + `Cmax_error(%)`)

table1 %>%
  kable() %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
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
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  ylim(0,3) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp4

#plot and save
gp <- grid.arrange(gp1, gp2, gp3, gp4, ncol=2, nrow=3)
ggsave(file="../deliv/fig4.pdf", gp, width=10, height=12)

###################################################################################################
###################################################################################################


###################################################################################################
########################### Chunk 3: Figure 4 (Sensitivity Plots) #################################
###################################################################################################
## This chunk reproduces Figure 4 sensitivity plots
## Model 1 with 4 mg/kg IV infusion dosing was used for all plots
bw <- 73
dose <- 4*bw
rate <- 4*bw

# Adipose
p1 <- model1 %>% knobs(Kpad = c(5,10,20), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpad", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("5","10", "20"), cex=0.5)
       ))
p1

# Bone
p2 <- model1 %>% knobs(Kpbo = c(4,8,16), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpbo", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("4","8", "16"), cex=0.5)
       ))
p2

# Brain
p3 <- model1 %>% knobs(Kpbr = c(3.5,7,14), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpbr", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("3.5","7", "14"), cex=0.5)
       ))
p3

# Gut
p4 <- model1 %>% knobs(Kpgu = c(3,6,12), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpgu", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("3","6", "12"), cex=0.5)
       ))
p4

# Heart
p5 <- model1 %>% knobs(Kpgu = c(1,2,4), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kphe", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1","2", "4"), cex=0.5)
       ))
p5

# Kidney
p6 <- model1 %>% knobs(Kpgu = c(1.5,3,6), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpki", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1.5","3", "6"), cex=0.5)
       ))
p6

# Liver
p7 <- model1 %>% knobs(Kpgu = c(2,4,8), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpli", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("2","4", "8"), cex=0.5)
       ))
p7

# Lungs
p8 <- model1 %>% knobs(Kplu = c(0.5,1,2), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kplu", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.5","1", "2"), cex=0.5)
       ))
p8

# Muscle
p9 <- model1 %>% knobs(Kpmu = c(1.5,3,6), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpmu", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1.5","3", "6"), cex=0.5)
       ))
p9

# Spleen
p10 <- model1 %>% knobs(Kpsp = c(1.5,3,6), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpsp", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("1.5","3", "6"), cex=0.5)
       ))
p10

# Rest of the body
p11 <- model1 %>% knobs(Kpre = c(2,4,8), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="Kpre", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("2","4", "8"), cex=0.5)
       ))
p11

# BP
p12 <- model1 %>% knobs(BP = c(0.5,1,2), amt = dose, cmt = "VEN", ii = 12, addl = 13, rate=rate, ss = 1, end = 12, delta = 0.1) %>% 
  plot(Cvenous~., main=list(label="BP", cex=0.75), xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), scale=list(cex=0.5), 
       lty=1:3, col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.5","1", "2"), cex=0.5)
       ))
p12

###################################################################################################
###################################################################################################


###################################################################################################
#################################### Chunk 4: Optimization ########################################
###################################################################################################
## This chunk runs parameter optimization for parameters: Kpmu, BP and sigma2 (residual error variance)
## Optimization was done using maximum likelihood (ML) with the adult 4 mg/kg IV infusion data
load("../data/Fig3a_obs.Rda") #load observed data (digitized) for figure 3a in ZT
#load("../data/Fig3a_ZT.Rda") #load Simcyp predictions (digitized) for figure 3a in ZT
sampl <- obs$time  #sampling times from observed data
wt <- 73
dose <- 4*wt
rate <- 4*wt

theta <- log(c(Kpmu=2.94, BP=1, sigma2=1))  #initial parameter estimates for MLE

#objective function calculator
OF <- function(pars,pred=FALSE){
  pars %<>% lapply(exp) #Get out of log domain for MLE
  names(pars) <- names(theta)
  
  ## Get a prediction
  out <- as.data.frame(model1 %>% 
                         param(pars) %>%
                         ev(amt = dose, cmt = "VEN", ii=12, addl=13, rate=rate, ss=1) %>%
                         mrgsim(end=-1, add=sampl) %>%
                         filter(row_number() != 1))
  
  if(pred) return(out)  #if pred=TRUE, just get out with simulation output (no optimization)
  
  ##maximum log likelihood
  return(-1*sum(dnorm(log(obs$obs),
                      log(out$Cvenous),
                      sqrt(pars$sigma2), log=TRUE)))
}

fit <- optim(theta,OF, method="Nelder", hessian=T)  #optimize the parameters using Nelder-Mead method
fpar <- exp(fit$par)  #final estimates of parameters from MLE
fpar
###################################################################################################
###################################################################################################


###################################################################################################
####################################### Chunk 5: Figure 5 #########################################
###################################################################################################
## This chunk reproduces Figure 5 plots
## Generate models 3 and 4; same as models 1 and 2 but using optimized Kpmu and BP from Chink 5
model3 <- param(model1, list(Kpmu=0.78, BP=1.21)) 
model4 <- param(model2, list(Kpmu=0.78, BP=1.21))  

## Figure 5a; Model 3 with 4 mg/kg IV infusion dosing
load("../data/Fig3a_obs.Rda")  ##load observed data (digitized) from fig 3a in the ZT paper
load("../data/Fig3a_ZT.Rda")  ##load ZT predictions (digitized) from fig 3a in the ZT paper
wt <- 73
dose <- 4*wt
rate <- 4*wt
cmt <- "VEN"

sim <- as.data.frame(model3 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, rate=rate, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp1 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 3")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("a   Adult 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp1

#Figure 5b;  Model 4 with 4 mg/kg IV infusion
load("../data/Fig4a_obs.Rda")  ##load observed data (digitized) from fig 4a in the ZT paper
load("../data/Fig4a_ZT.Rda")  ##load ZT predictions (digitized) from fig 4a in the ZT paper
wt <- 19
dose <- 4*wt
rate <- 3*wt
cmt <- "VEN"

sim <- as.data.frame(model4 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, rate=rate, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp2 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 4")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("b  Pediatric 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp2


#Figure 5c; Model 3 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

sim <- as.data.frame(model3 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp3 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 3")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("c   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp3


#Figure 5d; Model 4 with 4 mg/kg IV infusion 
load("../data/Fig4d_obs.Rda")  ##load observed data (digitized) from fig 4d in the ZT paper
load("../data/Fig4d_ZT.Rda")  ##load ZT predictions (digitized) from fig 4d in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

sim <- as.data.frame(model4 %>% 
                            ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                            mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp4 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 4")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("d  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  ylim(0,3) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp4


## For figures 5e and 5f we will first generate models 5 and 6 from models 3 and 4 by using optimized P_eff values
model5 <- param(model3, list(fperm=0.276))
model6 <- param(model4, list(fperm=0.276))

#Figure 5e; Model 5 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

sim <- as.data.frame(model5 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp5 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 5")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("e   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp5


#Figure 5f; Model 6 with 4 mg/kg PO 
load("../data/Fig4d_obs.Rda")  #load observed data (digitized) from fig 4d in the ZT paper
load("../data/Fig4d_ZT.Rda")  ##load ZT_Gu data (digitized) from fig 4e in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

sim <- as.data.frame(model6 %>% 
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp6 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 6")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("f  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  ylim(0,3) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp6

## For figures 5g and 5h we will first generate models 7 and 8 from models 5 and 6 by using MPPGI
## values calculated from corresponding MPPGL as: MPPGI = MPPGL/25 (see results in main text)
model7 <- param(model5, list(MPPGI=1.212))
model8 <- param(model6, list(MPPGI=1.04))

#Figure 5g; Model 7 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

sim <- as.data.frame(model7 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp7 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT","sim"),
                      labels=c("observed","ZT","model 7")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("g   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp7


#Figure 5f; Model 6 with 4 mg/kg PO 
load("../data/Fig4d_obs.Rda")  #load observed data (digitized) from fig 4d in the ZT paper
load("../data/Fig4d_ZT.Rda")  ##load ZT_Gu data (digitized) from fig 4e in the ZT paper
ZT1 <- ZT
load("../data/Fig4e_ZT.Rda")  ##load ZT_Gu data (digitized) from fig 4e in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

sim <- as.data.frame(model8 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp8 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = ZT, mapping = aes(x=time, y=ZT, col="ZT"), lwd=1, lty=2) +
  geom_line(data = ZT1, mapping = aes(x=time, y=ZT, col="ZT1"), lwd=1, lty=1) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  scale_colour_manual(name='', values=c('sim'='black', 'ZT'='grey', 'ZT1'='grey', 'observed'='black'), 
                      breaks=c("observed","ZT1","ZT","sim"),
                      labels=c("observed","ZT", expression(ZT[Gu]),"model 8")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,2,1), shape=c(16, NA, NA, NA)))) +
  ggtitle("h  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  ylim(0,3) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp8


###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 6: Figure S2 #########################################
###################################################################################################
### sensitivity analyses for absorption parameters: permeability, intestinal transit time and solubility ###
gp1 <- model3 %>%  
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

gp2 <- model3 %>% 
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

gp3 <- model3 %>%
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

###################################################################################################
###################################################################################################

###################################################################################################
####################################### Chunk 7: Figure S3 #########################################
###################################################################################################
### sensitivity analyses for intestinal clearance ###
p1 <- model5 %>%  
  knobs(MPPGI=c(1.212/2, 1.212, 1.212*2), amt = 200, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>% 
  plot(Cvenous~., 
       main=list(label=expression(Adult~Cl["Gu"]~(mL/min/kg)), cex=0.75), 
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
p1

p2 <- model6 %>%  
  knobs(MPPGI=c(1.04/2, 1.04, 1.04*2), amt = 4*19, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1, delta = 0.1, end = 12) %>% 
  plot(Cvenous~., 
       main=list(label=expression(Pediatric~Cl["Gu"]~(mL/min/kg)), cex=0.75), 
       xlab=list(label="time (h)", cex=0.75),
       ylab=list(label="Plasma concentration (mg/L)", cex=0.75), 
       scale=list(cex=0.5), 
       lty=1:3, 
       col="black",
       #auto.key=list(cex=2, lines=T, points=F),
       key=list(space="top",
                lines=list(col=c("black","black","black"), lty=1:3),
                text=list(c("0.095","0.19","0.38"), cex=0.5)
       ))
p2

###################################################################################################
###################################################################################################


###################################################################################################
####################################### Chunk 8: Figure S5 #########################################
###################################################################################################
## This chunk reproduces Figure S5 plots
## Figure S5a; Model 7 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "GUTLUMEN"

sim <- as.data.frame(model7 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

sim_lower <- as.data.frame(model7 %>%
                             param(MPPGI=0.61) %>%
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

sim_upper <- as.data.frame(model7 %>%
                             param(MPPGI=2.18) %>%
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp1 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  geom_line(data = sim_upper, mapping = aes(y=Cvenous, x=time, col="sim_upper"), lwd=1, lty=2) +
  geom_line(data = sim_lower, mapping = aes(y=Cvenous, x=time, col="sim_lower"), lwd=1, lty=2) +
  scale_colour_manual(name='', values=c('sim'='black', 'observed'='black', 'sim_upper'='black', 'sim_lower'='black'), 
                      breaks=c("observed","sim","sim_upper","sim_lower"),
                      labels=c("observed","model 7","","")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,0,0), shape=c(16, NA, NA, NA)))) +
  ggtitle("a   Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp1


#Figure S5b; Model 8 with 4 mg/kg PO 
load("../data/Fig4d_obs.Rda")  #load observed data (digitized) from fig 4d in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "GUTLUMEN"

sim <- as.data.frame(model8 %>% 
                       ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                       mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

sim_lower <- as.data.frame(model8 %>%
                             param(MPPGI=0.52) %>%
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

sim_upper <- as.data.frame(model8 %>%
                             param(MPPGI=1.89) %>%
                             ev(amt=dose, cmt=cmt, ii=12, addl=7, ss=1) %>% 
                             mrgsim(delta = 0.1, end = 12)) %>%
  dplyr::filter(row_number() != 1)  #sets the mrgsolve dataframe

gp2 <- ggplot() + geom_point(data = obs, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
  geom_errorbar(data = obs, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
  geom_line(data = sim, mapping = aes(y=Cvenous, x=time, col="sim"), lwd=1) + 
  geom_line(data = sim_upper, mapping = aes(y=Cvenous, x=time, col="sim_upper"), lwd=1, lty=2) +
  geom_line(data = sim_lower, mapping = aes(y=Cvenous, x=time, col="sim_lower"), lwd=1, lty=2) +
  scale_colour_manual(name='', values=c('sim'='black', 'observed'='black', 'sim_upper'='black','sim_lower'='black'), 
                      breaks=c("observed","sim","sim_upper","sim_lower"),
                      labels=c("observed","model 8","","")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,0,0), shape=c(16, NA, NA, NA)))) +
  ggtitle("b  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  ylim(0,3) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp2


###################################################################################################
###################################################################################################
