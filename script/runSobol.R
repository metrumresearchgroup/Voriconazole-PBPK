## This script runs the Sobol global sensitivity analysis to reproduce Figure S1
## CAUTION: This script takes several hours to run
## load libraries
#loading libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(mrgsolve)  #https://github.com/metrumresearchgroup/mrgsolve
library(sensitivity)
source("genSamples.R")

#compile model1
model1 <- mread_cache("model1", "../model")  # model1 (Adult physiology)


###################################################################################################
##################################### Figure S1 ##########################################
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
    ev(amt = dose, cmt = "VEN", ii = 12, rate = rate, ss = 1) %>%
    mrgsim(end=2, delta=0.1) %>%
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

#generate random samples
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
pt1 <- proc.time()
globSens <- sobol2007(modRun, X1=X1, X2=X2, nboot=100)
pt2 <- proc.time() - pt1

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

gp_sobol <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  labs(y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) 
gp_sobol

#save.image(file = "sobol_all_objects.RData")  #can be uncommented if you want to save workspace objects

###################################################################################################
###################################################################################################

###################################################################################################
##################################### Figure S4 ##########################################
###################################################################################################
## This chunk reproduces Figure S4 Sobol global sensitivity analysis
## Analysis is done on model 1 with 4 mg/kg IV infusion dosing
## CAUTION: This chunk takes several hours to run
wt <- 73
dose <- 200

oneRun <- function(pars){
  out <- model1 %>%
    param(pars) %>%
    ev(amt = dose, cmt = "GUTLUMEN", ii = 12, ss = 1) %>%
    mrgsim(end=2, delta=0.1) %>%
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

#generate random samples
set.seed(1751)
n <- 20000  #number of samples per parameter
l <- list(fperm=c(1/3,1*3),
          S_lumen=c(390/3,390*3),
          ITT=c(3.32/3,3.32*3),
          MPPGI=c(1.44/3,1.44*3))
X1 <- genSamples(n, l)  #generate first random samples set
X2 <- genSamples(n, l)  #generate second random samples set

#run sobol function for Cmax endpoint
pt1 <- proc.time()
globSens <- sobol2007(modRun, X1=X1, X2=X2, nboot=100)
pt2 <- proc.time() - pt1

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
  dplyr::mutate(params = params) %>%
  dplyr::mutate(params = case_when(params == "fperm" ~ "P_eff",
                                   params == "S_lumen" ~ "S_int",
                                   params == "MPPGI" ~ "Cl_Gu",
                                   TRUE ~ "ITT"))

gp_sobol <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  labs(y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) 
gp_sobol

#save.image(file = "sobol_all_objects.RData")  #can be uncommented if you want to save workspace objects

###################################################################################################
###################################################################################################
