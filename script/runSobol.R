## This script runs the Sobol global sensitivity analysis to reproduce Figure S1
## CAUTION: This script takes several hours to run
## load libraries
#loading libraries
.libPaths("lib")
library(dplyr)
library(ggplot2)
library(PKPDmisc)
library(gridExtra)
library(mrgsolve)  #https://github.com/metrumresearchgroup/mrgsolve
library(sensitivity)
source("genSamples.R")
source("CalcKp_R&R.R")

#compile model1
model1 <- mread_cache("model1", "../model")  # model1 (Adult physiology)


###################################################################################################
##################################### Figure S1 ##########################################
###################################################################################################
## This chunk reproduces Figure S1 Sobol global sensitivity analysis
## Analysis is done on model 1 with 4 mg/kg IV infusion dosing
## CAUTION: This chunk takes several hours to run

## drug-related parameters
type <- 3
logP <- 2.56
pKa <- 1.76
fup <- .42
BP <- 1

## calculate Kps and update model1
dat <- read.csv("../data/tissue_comp_R&R.csv")
Kp_RR <- calcKp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type)
model1 <- param(model1, Kp_RR)

## run analysis
wt <- 73
dose <- 200

#########################################################################################
### sobol run for Cmax ; if you uncomment this block you need to comment the AUC block ##
# oneRun <- function(pars){
#   out <- model1 %>%
#     param(pars) %>%
#     ev(amt = dose, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1) %>%
#     mrgsim(end=2, delta=0.1) %>%
#     as_data_frame() %>%
#     summarise(Cmax=max(Cvenous))
#   return(out$Cmax)
# }
##########################################################################################

#########################################################################################
### sobol run for AUC ; if you uncomment this block you need to comment the Cmax block###
oneRun <- function(pars){
  out <- model1 %>%
    param(pars) %>%
    ev(amt = dose, cmt = "GUTLUMEN", ii = 12, addl = 13, ss = 1) %>%
    mrgsim(end=12, delta=0.1) %>%
    as_data_frame() %>%
    dplyr::filter(row_number() != 1) #%>%
    #dplyr::summarise(AUC=auc_partial(time, Cvenous))
  AUC <- as.numeric(auc_partial(out$time, out$Cvenous))
  return(AUC)
}
#########################################################################################

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

#print and assemble results
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

save.image(file = "sobol_all_objects_auc.RData")  #can be uncommented if you want to save workspace objects

# plot
gp_sobol_auc <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  labs(y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) 
gp_sobol


## Uncomment the following block to reproduce the paper FigS2
## FigS2
#plot and save
load("sobol_all_objects_auc_previous.RData")  #load results from previous run
gp_sobol_auc <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  ylim(-0.01,1) +
  labs(title="AUC", y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5)) 
gp_sobol_auc

load("sobol_all_objects_cmax_previous.RData")  #load results from previous run
gp_sobol_cmax <- ggplot(data=st, aes(x=params, y=original, group=index, shape=index)) +
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(x=params, ymin=min, ymax=max), width=0, position=position_dodge(width=0.5)) +
  ylim(-0.01,1) +
  labs(title="Cmax", y="Index value", x="parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14),
        plot.title = element_text(hjust = 0.5)) 
gp_sobol_cmax

# save
g <- grid.arrange(gp_sobol_auc, gp_sobol_cmax, ncol=2, nrow=3)
ggsave(file="../deliv/fig6.pdf", g, width=10, height=12)
###################################################################################################
###################################################################################################
