## clear workspace
rm(list=ls())

## set working directory
setwd("~/github/JMdiscfrail")

## load data
library(data.table)
load("data/fake_dataRD.Rdata")

## load functions
source('functions/data_format.R')
source('functions/JMdiscfrail.R')
source('functions/stratified_base_surv.R')


# Formatting data
dataR = formatting.data(dataR)
dataD = formatting.data(dataD)


###########################
# Gaussian Initialization #
###########################
SigmaG =  matrix(c(2*0.12,0.0,0.0,2*1.4),nrow = 2, ncol = 2)
muG = c(0,0)

result = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init = 'gauss',
                     distance = "euclidean",
                     Sigma = SigmaG,
                     mu = muG,
                     M = 1000, L = 1,
                     max.it = 10,  toll = 1e-3)


#---------#
# RESULTS #
#---------#

# Fixed effects - Recurrent (betas)
result$modelR
# Fixed effects - Terminal (gammas)
result$modelD

# Baseline Cumulative Hazard for Recurrent event
result$cumhazR
# Baseline Cumulative Hazard for Terminal event
result$cumhazD

# Random effects
result$K
result$w
result$P
# Subjects' subgroups
head(result$id.subgroup)
table(result$id.subgroup$subgroup)

# Estimated mass-points
masses = paste0('(',sprintf("%.3f", round(result$P[,1],3)),', ',
                sprintf("%.3f", round(result$P[,2],3)),')')

dev.new()
plot(result$P[,1], result$P[,2], cex=result$w*10,
     xlim = range(result$P[,1]), ylim = range(result$P[,2]),
     col=1:result$K, pch=16,
     xlab='Recurrent (u)', ylab='Terminal (v)', 
     main='Distribution of masses',
     cex.axis=1.3, cex.lab = 1.3, cex.main=1.3)
abline(h=0, lty=2)
abline(v=0, lty=2)
legend(0.1,-0.1, col=1:(result$K), 
       pch = rep(16,result$K), y.intersp=1,
       legend = paste0('P',seq(1,result$K),': ',masses,' - w',seq(1,result$K),': ',round(result$w,3)))


#----------------------------#
# STRATIFIED SURVIVAL CURVES #
#----------------------------#
survival = stratified.base.surv(result)

dev.new()
par(mfrow=c(1,2))
survR = as.data.frame(survival$recurrent)
plot(survR[,1], survR[,2], type='l', lty=2,  ylim =c(0,1),
     xlab = 'Time', ylab='Survival Probability',
     main = 'Stratified baseline survival curve for recurrent events')
for(k in 1:(result$K)){
  points(survR[,1], survR[,(k+2)], type='l', col=k)
}
legend('topright', col=c('black',1:(result$K)), 
       lty = c(2,rep(1,result$K)),
       legend = colnames(survR[,2:(result$K+2)]))

survD = as.data.frame(survival$terminal)
plot(survD[,1], survD[,2], type='l', lty=2, ylim =c(0,1),
     xlab = 'Time', ylab='Survival Probability',
     main = 'Stratified baseline survival curve for terminal event')
for(k in 1:(result$K)){
  points(survD[,1], survD[,(k+2)], type='l', col=k)
}
legend('topright', col=c('black',1:(result$K)), 
       lty = c(2,rep(1,result$K)),
       legend = colnames(survD[,2:(result$K+2)]))
