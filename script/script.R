#This file contains the scripts used for simulations and plots for the voriconazole PBPK manuscript
#First thing is to set working directory to "script" directory using the function: setwd()
#CAUTION: Chunk 4 takes several hours to run so you might want to comment this out before running the whole script.
#Note: All observed data and Zane and Thakker predictions were digitized from Zane and Thakker (ZT) paper. Zane and Thakker. Clin Pharmacokinet (2014) 53:1171–1182
#Note: All digitization was done using WebPlotDigitizer version 3.12 developed by Ankit Rohatgi https://automeris.io/WebPlotDigitizer

#loading libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(mrgsolve)  #https://github.com/metrumresearchgroup/mrgsolve
library(sensitivity)
source("calcKp_PT.R")
source("genSamples.R")

#compile models
model1 <- mread_cache("model1", "../model")  # model1 (Adult physiology)
model2 <- mread_cache("model2", "../model")  # model2 (Pediatric physiology)


###################################################################################################
######################## Chunk 1: Calculate partition coefficients ################################
###################################################################################################
## This chunk calculates tissue:plasma partition coefficients used in models 1&2 using Poulin & Theil
## method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
# First define Voriconazole physicochemical properties
logP <- 2.56  #oil:water partition
pKa <- 1.76
fup <-.42  #fraction unbound in plasma

# calculate partition coefficients
Kp <- calcKp(logP=logP, pKa=pKa, fup=fup, type=3)  #type=3 for monoprotic base
Kp
###################################################################################################
###################################################################################################


###################################################################################################
####################################### Chunk 2: Figure 3 #########################################
###################################################################################################
## This chunk reproduces Figure 3 plots; all observed and ZT data were digitized and saved as .Rda 
## files whose numbers follow the ZT publication numbering

## Figure 3a; Model 1 with 4 mg/kg IV infusion ##
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
cmt <- "D"

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
cmt <- "D"

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
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp4

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
##################################### Chunk 4: Figure S1 ##########################################
###################################################################################################
## This chunk reproduces Figure S1 Sobol global sensitivity analysis
## Analysis is done on model 1 with 4 mg/kg IV infusion dosing
## CAUTION: This chunk takes several hours to run
wt <- 73
dose <- 4*wt
rate <- 4*wt

oneRun <- function(pars){
  out <- model1 %>%
    param(pars) %>%
    ev(amt = dose, cmt = "VEN", ii = 12, addl = 13, rate = rate, ss = 1) %>%
    mrgsim(end=12, delta=0.1) %>%
    as_data_frame() %>%
    summarise(Cmax=max(Cvenous))
  return(out$Cmax)
}

#sobol function to get entire vector of outputs for every parameter combination
modRun <- function(X){
  l <- split(X, seq(nrow(X)))   #split dataframe into list of nested lists of each row of df
  op <- mapply(oneRun, l)  #apply the oneRun function to each parameter combination in the list
  op <- as.vector(op, mode="numeric")
  return(op)
}

#generate random samples using pseudo=random generator
set.seed(1751)
n <- 20000  #number of samples per parameter
l <- list(BP=c(1/3,1*3),
          Kpmu=c(2.82/3,2.82*3),
          Kpbo=c(7.9/3,7.9*3),
          Kpki=c(2.9/3,2.9*3),
          Kpli=c(4.66/3,4.66*3),
          Kplu=c(0.83/3,0.83*3),
          Kpsp=c(2.96/3,2.96*3),
          Kpad=c(9.89/3,9.89*3),
          Kpbr=c(7.35/3,7.35*3),
          Kpgu=c(5.82/3,5.82*3),
          Kphe=c(1.95/3,1.95*3)) #list of param ranges
X1 <- genSamples(n, l)  #generate first random samples set
X2 <- genSamples(n, l)  #generate second random samples set

#run sobol function for Cmax endpoint
globSens <- sobol2007(modRun, X1=X1, X2=X2, nboot=100)

#print and plot
print(globSens)

s <- globSens$S  #single indices
t <- globSens$T  #total indices

s$index <- rep(0, nrow(s))  #add flag to identify
t$index <- rep(1, nrow(t))

st <- bind_rows(s,t) %>%
  dplyr::mutate(index = as.factor(index)) %>%
  dplyr::rename(min=`min. c.i.`, max=`max. c.i.`)
levels(st$index) <- c("first-order", "total")
params <- as.character(rownames(s))

st <- st %>%
  group_by(index) %>%
  dplyr::mutate(params = params)

gp <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  labs(y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) 
gp

###################################################################################################
###################################################################################################


###################################################################################################
#################################### Chunk 5: Optimization ########################################
###################################################################################################
## This chunk runs parameter optimization for parameters: Kpmu, BP and sigma2 (residual error variance)
## Optimization was done using maximum likelihood (ML) with the adult 4 mg/kg IV infusion data
load("../data/Fig3a_obs.Rda") #load observed data (digitized) for figure 3a in ZT
load("../data/Fig3a_ZT.Rda") #load Simcyp predictions (digitized) for figure 3a in ZT
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
####################################### Chunk 6: Figure 5 #########################################
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
cmt <- "D"

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
cmt <- "D"

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
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp4


## For figures 5e and 5f we will first generate models 5 and 6 from models 3 and 4 by using MPPGI
## values calculated from corresponding MPPGL as: MPPGI = MPPGL/21 (see results in main text)
model5 <- param(model3, list(MPPGI=1.44))
model6 <- param(model4, list(MPPGI=1.24))

#Figure 5e; Model 5 with 200 mg PO 
load("../data/Fig3b_obs.Rda")  ##load observed data (digitized) from fig 3b in the ZT paper
load("../data/Fig3b_ZT.Rda")  ##load ZT predictions (digitized) from fig 3b in the ZT paper
wt <- 73
dose <- 200
cmt <- "D"

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
load("../data/Fig4e_ZT.Rda")  ##load ZT_Gu data (digitized) from fig 4e in the ZT paper
wt <- 19
dose <- 4*wt
cmt <- "D"

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
                      labels=c("observed","ZT_Gu","model 6")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(0,1,1), shape=c(16, NA, NA)))) +
  ggtitle("f  Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
  scale_y_continuous(breaks = seq(0,100,1)) +
  #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
  scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
  theme_bw() + 
  theme(axis.title=element_text(size=12)) +
  theme(axis.text=element_text(size=10)) + 
  theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
  theme(plot.title=element_text(size=15, face="bold")) +
  theme(legend.text=element_text(size=10))
gp6

###################################################################################################
###################################################################################################
