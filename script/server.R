.libPaths("lib")
library(Rcpp)
library(dplyr)
library(mrgsolve)
library(shiny)
library(gridExtra)
shinyServer(
  function(input,output){
    load("../data/Fig3a_obs.Rda")  #load observed data (digitized) from fig 3a in the ZT paper
    obs1 <- obs
    load("../data/Fig4a_obs.Rda")  #load observed data (digitized) from fig 3a in the ZT paper
    obs2 <- obs
    load("../data/Fig3b_obs.Rda")  #load observed data (digitized) from fig 3b in the ZT paper
    obs3 <- obs
    load("../data/Fig4d_obs.Rda")  ##load observed data (digitized) from fig 4d in the ZT paper
    obs4 <- obs
    modA <- mread("../model/model1")
    modP <- mread("../model/model2")
    Adult_MPPGI <- reactive({
      ifelse(input$optimAdult_ClGu, 0.065*73*9.3*0.711/(40*0.65), as.numeric(input$Adult_ClGu)*73*9.3*0.711/(40*0.65))
    })
    Pediatric_MPPGI <- reactive({
      ifelse(input$optimPediatric_ClGu, 0.19*19*11*0.711/(120.5*0.22), as.numeric(input$Adult_ClGu)*19*11*0.711/(40*0.22))
    })
    Kpmu <- reactive({
      ifelse(input$optimKpmu, 0.78, as.numeric(input$Kpmu))
    })
    BP <- reactive({
      ifelse(input$optimBP, 1.21, as.numeric(input$BP))
    })
    fperm <- reactive({
      ifelse(input$optimPeff, 0.276, as.numeric(input$Peff)/0.145)
    })
    
    output$myplot <- renderPlot({
      parsA <- list(Kpmu=Kpmu(), BP=BP(), fperm=fperm(), MPPGI=Adult_MPPGI()) 
      parsP <- list(Kpmu=Kpmu(), BP=BP(), fperm=fperm(), MPPGI=Pediatric_MPPGI())
      outAIV <- as.data.frame(modA %>%
                            param(parsA) %>%
                            ev(cmt = "VEN", amt = 4*73, ii = 12, addl =  13, rate = 4*73, ss=1) %>%
                            mrgsim(end = 12, delta = 0.1) %>%
                              dplyr::filter(row_number() > 1))
      
      outPIV <- as.data.frame(modP %>%
                                param(parsP) %>%
                                ev(cmt = "VEN", amt = 4*19, ii = 12, addl =  13, rate = 4*19, ss=1) %>%
                                mrgsim(end = 12, delta = 0.1) %>%
                                dplyr::filter(row_number() > 1))
      
      outAPO <- as.data.frame(modA %>%
                                param(parsA) %>%
                                ev(cmt = "GUTLUMEN", amt = 200, ii = 12, addl =  13, ss=1) %>%
                                mrgsim(end = 12, delta = 0.1) %>%
                                dplyr::filter(row_number() > 1))
      
      outPPO <- as.data.frame(modP %>%
                                param(parsP) %>%
                                ev(cmt = "GUTLUMEN", amt = 4*19, ii = 12, addl =  13, ss=1) %>%
                                mrgsim(end = 12, delta = 0.1) %>%
                                dplyr::filter(row_number() > 1))
      
      gp1 <- ggplot() + geom_point(data = obs1, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
        geom_errorbar(data = obs1, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
        geom_line(data = outAIV, mapping = aes(y=Cvenous, x=time, col="outAIV"), lwd=1) + 
        scale_colour_manual(name='', values=c('outAIV'='black', 'observed'='black'), 
                            breaks=c("observed","outAIV"),
                            labels=c("observed","Predicted")) +
        guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
        ggtitle("Adult 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
        scale_y_continuous(breaks = seq(0,100,1)) +
        #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
        scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
        theme_bw() + 
        theme(axis.title=element_text(size=12)) +
        theme(axis.text=element_text(size=10)) + 
        theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
        theme(plot.title=element_text(size=15, face="bold")) +
        theme(legend.text=element_text(size=10))
 
      gp2 <- ggplot() + geom_point(data = obs2, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
        geom_errorbar(data = obs2, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
        geom_line(data = outPIV, mapping = aes(y=Cvenous, x=time, col="outPIV"), lwd=1) + 
        scale_colour_manual(name='', values=c('outPIV'='black', 'observed'='black'), 
                            breaks=c("observed","outPIV"),
                            labels=c("observed","Predicted")) +
        guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
        ggtitle("Pediatric 4 mg/kg IV") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
        scale_y_continuous(breaks = seq(0,100,1)) +
        #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
        scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
        theme_bw() + 
        theme(axis.title=element_text(size=12)) +
        theme(axis.text=element_text(size=10)) + 
        theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
        theme(plot.title=element_text(size=15, face="bold")) +
        theme(legend.text=element_text(size=10))
      
      gp3 <- ggplot() + geom_point(data = obs3, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
        geom_errorbar(data = obs3, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
        geom_line(data = outAPO, mapping = aes(y=Cvenous, x=time, col="outAPO"), lwd=1) + 
        scale_colour_manual(name='', values=c('outAPO'='black', 'observed'='black'), 
                            breaks=c("observed","outAPO"),
                            labels=c("observed","Predicted")) +
        guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
        ggtitle("Adult 200 mg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
        scale_y_continuous(breaks = seq(0,100,1)) +
        #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
        scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
        theme_bw() + 
        theme(axis.title=element_text(size=12)) +
        theme(axis.text=element_text(size=10)) + 
        theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
        theme(plot.title=element_text(size=15, face="bold")) +
        theme(legend.text=element_text(size=10))
      
      gp4 <- ggplot() + geom_point(data = obs4, mapping = aes(x=time, y=obs, col="observed"), size=2.5) + 
        geom_errorbar(data = obs4, mapping = aes(x = time, y = obs, ymin=obs-sd, ymax=obs+sd), width=0) +
        geom_line(data = outPPO, mapping = aes(y=Cvenous, x=time, col="outPPO"), lwd=1) + 
        scale_colour_manual(name='', values=c('outPPO'='black', 'observed'='black'), 
                            breaks=c("observed","outPPO"),
                            labels=c("observed","Predicted")) +
        guides(colour = guide_legend(override.aes = list(linetype=c(0,1), shape=c(16, NA)))) +
        ggtitle("Pediatric 4 mg/kg PO") + xlab("time (h)") + ylab("Plasma concentration (mg/L)") +
        scale_y_continuous(breaks = seq(0,100,1)) +
        #scale_y_continuous(trans="log", breaks = c(0.1,1,10), limits = c(0.1,10)) +
        scale_x_continuous(breaks = seq(0,12,2), limits = c(0,12)) +
        theme_bw() + 
        theme(axis.title=element_text(size=12)) +
        theme(axis.text=element_text(size=10)) + 
        theme(axis.title.y=element_text(margin=margin(0,15,0,0))) +
        theme(plot.title=element_text(size=15, face="bold")) +
        theme(legend.text=element_text(size=10))
      
      gp <- grid.arrange(gp1, gp2, gp3, gp4, ncol=2, nrow=2)
    })
  }
)




# .libPaths("lib")
# library(dplyr)
# library(ggplot2)
# library(mrgsolve)
# library(magrittr)
# library(shiny)

# shinyServer(
#   function(input, output){
#     
#     pcoeff <- reactive({
#       pred <- ifelse(input$prediction == "Poulin & Theil", "P&T", 
#                      ifelse(input$prediction == "Berezhkovskiy", "Berez", 
#                             ifelse(input$prediction == "Schmitt", "Schmitt",
#                                    ifelse(input$prediction == "Rodgers & Rowland", "R&R", "pksim"))))
#       type <- ifelse(input$compType == "Acid", 3, 2)
#       pcoeffs(logP=as.numeric(input$logP), pKa=as.numeric(input$pKa), fup=as.numeric(input$fup), pred=pred, type=type)  #get partition coefficients
#     })
#     
#     pars <- reactive({
#       if(input$simType == "Individual"){
#         if(input$Gender == "Male"){
#           is.male <- TRUE
#         }else{
#           is.male <- FALSE
#         }
#         pars <- genInd_shiny(age=as.numeric(input$Age), is.male=is.male, bw_targ=as.numeric(input$Weight), ht_targ=as.numeric(input$Height))  #get the individual parameters
#         pars
#       }else{
#         lpars <- genPop_shiny(nSubj=as.numeric(input$nSubj),
#                               minAge=as.numeric(input$AgePop[1]), maxAge=as.numeric(input$AgePop[2]), 
#                               femPerc=as.numeric(input$femalePerc),
#                               minBW=as.numeric(input$WeightPop[1]), maxBW=as.numeric(input$WeightPop[2]), 
#                               minHT=as.numeric(input$HeightPop[1]), maxHT=as.numeric(input$HeightPop[2]))
#         ##single core simulations
#         for(i in seq_along(lpars)) lpars[[i]]["rep"] <- i  #add indices to each subject
#         lpars
#       }
#     })
#     
#     ##set the steady state
#     ss <- reactive({
#       ifelse(input$SS, 1, 0)
#     })
#     
#     ##set route of administration
#     cmt <- reactive({
#       ifelse(input$ROI == "IV", "VEN", "GUT")
#     })
#     
#     ##simulation plot
#     output$myplot1 <- renderPlot({
#       if(input$simType == "Individual"){
#         pars <- c(pars(), pcoeff())  #get full indivual parameters
#         df <- as.data.frame(mod %>%
#                               param(pars) %>%
#                               ev(cmt = cmt(), amt = as.numeric(input$Dose)*as.numeric(input$Weight), ii = as.numeric(input$II), addl =  as.numeric(input$ADDL),
#                                  rate = as.numeric(input$Rate)*as.numeric(input$Weight), ss=ss()) %>%
#                               mrgsim(end = as.numeric(input$end), delta = as.numeric(input$delta)))
#         df <- df[-1,]
#         ggplot(data=df, aes(x=time, y=Cvenous)) +
#           geom_line() +
#           labs(x="time (h)", y="conc (mg/L)")
#       }else{
#         pars <- lapply(pars(), function(x) x <- c(x,pcoeff()))  #get full indivual parameters by adding pcoeffs to all nested lists
#         sim <- function(pars, mod){
#           mod %>%
#             param(pars) %>%
#             ev(cmt = cmt(), amt = as.numeric(input$Dose)*pars$BW, rate = as.numeric(input$Dose)*pars$BW, 
#                ii = as.numeric(input$II), addl =  as.numeric(input$ADDL), ss = ss()) %>% 
#             mrgsim(end = as.numeric(input$end), delta = as.numeric(input$delta)) %>%
#             mutate(rep=pars$rep, sex=pars$SEX, bw=pars$BW, ht=pars$HT)
#         }
#         
#         ##bit of data manipulation
#         sims <- lapply(pars, sim, mod=mod) %>% bind_rows %>% group_by(rep) %>% dplyr::filter(row_number() != 1)
#         sims$sex <- as.factor(sims$sex)
#         levels(sims$sex)
#         levels(sims$sex) <- c("male","female")
#         sims <- sims %>% group_by(time) %>% mutate(avg=mean(Cvenous), low=quantile(Cvenous, probs=0.025), high=quantile(Cvenous, probs=0.975))  #get mean and 95% interval for each timepoint
#         #sims2 <- sims %>% group_by(time, sex) %>% mutate(avg=mean(Cvenous), low=quantile(Cvenous, probs=0.025), high=quantile(Cvenous, probs=0.975))  #get mean and 95% interval for each timepoint stratified by sex
#         
#         #plot
#         ggplot(data=sims, aes(x=time)) + 
#           geom_line(aes(y=avg)) + 
#           geom_ribbon(aes(ymin=low, ymax=high), alpha=0.2) +
#           labs(x="time (h)", y="conc (mg/L)")
#       }
#     })
#     
#     ##administration schedule output
#     output$myplot2 <- renderPlot({
#       dur <- as.numeric(input$II)*(as.numeric(input$ADDL) + 1)  #total duration of therpay
#       x <- seq(0, dur, 1)  #set the time axis
#       ivals <- cumsum(rep(as.numeric(input$II), dur/as.numeric(input$II)))  #intervals
#       y <- rep(0, length(x))
#       y[c(1,ivals)] <- as.numeric(input$Dose)
#       df3 <- data.frame(time=x, Dose=y)
#       ggplot(data=df3, aes(x=time, y=Dose)) + 
#         geom_line() +
#         labs(x="time(h)", y="Dose (mg/kg)")
#     })
#     
#     ##Individual output
#     output$mytable1 <- renderTable({
#       if(input$simType == "Individual"){
#         df1 <- data.frame(parameter=names(pars()), value=as.numeric(pars()))
#         df1
#       }
#     })  
#     
#     ##Population output
#     output$mysummary <- renderPrint({
#       if(input$simType == "Population"){
#         ltmp <- unlist(pars())
#         df2 <- as.data.frame(split(ltmp, names(ltmp)))
#         summary(df2)
#       }
#     })
#     
#     ##Compound output
#     output$mytable2 <- renderTable({
#       df3 <- data.frame(parameter=names(pcoeff()), value=as.numeric(pcoeff()))
#       df3
#     })
#     
#   }
# )
