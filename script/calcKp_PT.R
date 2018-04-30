#This script calculates tissue:plasma partition coefficients based on: Poulin and Theil http://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext

calcKp <- function(logP, pKa, fup, BP=1, type=1){
  
  dat <- read.csv("../data/tissue_comp_PT.csv")
  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)
  
  Vwp <- dat$FVwt[dat$tissue == "Plasma"]
  Vnlp <- dat$FVnlt[dat$tissue == "Plasma"]
  Vphp <- dat$FVpht[dat$tissue == "Plasma"]
  
  dat2 <- dat %>% filter(!tissue %in% c("Plasma"))
  
  Vwt <- dat2$FVwt[dat2$tissue != "Adipose"]
  Vwad <- dat2$FVwt[dat2$tissue == "Adipose"]
  Vnlt <- dat2$FVnlt[dat2$tissue != "Adipose"]
  Vnlad <- dat2$FVnlt[dat2$tissue == "Adipose"]
  Vpht <- dat2$FVpht[dat2$tissue != "Adipose"]
  Vphad <- dat2$FVpht[dat2$tissue == "Adipose"]
    
  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species
  
  logD_star <- switch(type,
                 #1-neutral
                 logD,   
                 #2-monoprotic acid
                 logD-log10(1+10^(pH-pKa)),
                 #3-monoprotic base
                 logD-log10(1+10^(pKa-pH)), 
                 #4-diprotic acid
                 logD-log10(1+10^(2*pH-pKa[1]-pKa[2])),
                 #5-diprotic base
                 logD-log10(1+10^(pKa[1]+pKa[2]-2*pH)), 
                 #6-monoprotic acid monoprotic base (acid comes first)
                 logD-log10(1+10^(pKa[2]-pKa[1])),  
                 #7-triprotic acid
                 logD-log10(1+10^(3*pH-pKa[1]-pKa[2]-pKa[3])),  
                 #8-triprotic base
                 logD-log10(1+10^(pKa[1]+pKa[2]+pKa[3]-3*pH)),  
                 #9-diprotic acid monoprotic base (first two are acid)
                 logD-log10(1+10^(pH-pKa[1]-pKa[2]+pKa[3])), 
                 #10-diprotic base monoprotic acid (first one is acid)
                 logD-log10(1+10^(pKa[2]+pKa[3]-pKa[1]-pH)))       

  D_star <- 10^logD_star   
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+(1*(Vwad+0.7*Vphad)))/(D_star*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * fup
      
  P <- 10^logP
  fut <- 1/(1+((1-fup)/fup)*0.5)
  Kpt <- ((P*(Vnlt+0.3*Vpht)+(1*(Vwt+0.7*Vpht)))/(P*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * (fup/fut)
  
  Kp <- c(Kpad, Kpt)
  name <- dat2$tissue %>% substr(1,2) %>% tolower()
  name <- paste("Kp", name, sep="")
  uParam <- split(Kp, name)
  
  return(uParam)
}
    
    
    
    



