########################################################################
# Joint models with parametric frailty
########################################################################
# Main script for the application based on the fake recurrent-event
# and terminal-event data.
#
# The script fits:
#   (i) the Rondeau et al. joint frailty model with log-normal frailty;
#   (ii) the Ng et al. joint frailty model.
#
# It reproduces the application results reported in the manuscript,
# including Table 1 and Figure 6.
#
# All numerical results and manuscript figures/tables are saved in the
# 'results/' directory.
########################################################################

## clear workspace
rm(list=ls())

## set working directory
setwd("~/github/JMDF")

## load data
library(data.table)
library(frailtypack)

## load functions
source('functions/data_format.R')
source('functions/JointFrailtyNg.R')
source('functions/jmdf_plots.R')

###############################
# Rondeau et al's joint model #
###############################
load("data/fake_dataRD.Rdata")

# Formatting data
dataR = formatting.data(dataR)
dataD = formatting.data(dataD)


## Disjoint models via frailtypack (to identify kappa values)
##-----------------------------------------------------------
## Cox Model with random effect -- Recurrent
coxR <- frailtyPenal(Surv(gaptime,deltaR) ~ cluster(id) + sex + age + 
                       ncom + adherent, n.knots=14, kappa=1, 
                     data=dataR, cross.validation = TRUE)
## Cox Model with random effect -- Death
coxD <- frailtyPenal(Surv(gaptime,deltaD) ~ cluster(id) + sex + age + 
                       ncom + adherent, n.knots=14, kappa=1, 
                     data=dataD, cross.validation = TRUE)
## Selection of smoothing parameters via disjoint models
kappa1 <- coxR$kappa
kappa2 <- coxD$kappa


## Joint frailty model -- LogNormal
##----------------------------------------------------------------------------
jm.rondeau <- frailtyPenal(Surv(gaptime,deltaR) ~ cluster(id) + sex + age + 
                             ncom + adherent + terminal(deltaD),
                           formula.terminalEvent = ~ sex + age + ncom + adherent,
                           recurrentAG = FALSE,
                           data=dataR, n.knots=5, 
                           kappa=c(kappa1,kappa2), 
                           RandDist = "LogN" )

print(jm.rondeau)
summary(jm.rondeau)

jm.rondeau$coef

rand.eff.rondeau = rbind.data.frame(
  cbind.data.frame(
    Model = "Rondeau et al",
    Parameter = "sigma2",
    Estimate = jm.rondeau$sigma2,
    StdDev = 2 * sqrt(jm.rondeau$sigma2) * sqrt(diag(jm.rondeau$varH))[1],
    pvalue = jm.rondeau$sigma2_p.value),
  cbind.data.frame(
    Model = "Rondeau et al",
    Parameter = "alpha",
    Estimate = jm.rondeau$alpha,
    StdDev = sqrt(diag(jm.rondeau$varH))[2],
    pvalue = jm.rondeau$alpha_p.value)
)
rand.eff.rondeau

##########################
# Ng et al's joint model #
##########################
# Reload the original data because the data objects were reformatted
# above for use with frailtypack.
rm(dataR)
rm(dataD)
load("data/fake_dataRD.Rdata")

# The data are converted to the format required by the
# implementation of Ng et al.'s model.
dat.ng = dataR[, .(event,id,gaptime,deltaR,deltaD,sex,age,ncom,adherent)]
colnames(dat.ng) = c("event","id","time","delta1","delta2","sex","age","ncom","adherent")
dat.ng = formatting.data.ng(dat.ng)
dat.ng = as.matrix(dat.ng)

jm.ng = joint.frailty.Ng(dat.ng, patient=dat.ng[,2], 
                         theta01=0.1, theta02=0.1, rho0=0.5, itmax=15)

# Fixed effects - Recurrent (betas)
jm.ng$Recurrent
# Fixed effects - Terminal (gammas)
jm.ng$Death

# Random effects
jm.ng$Frailty
# Mass plot
plot(jm.ng$u, jm.ng$v, xlab='u', ylab='v', main='Random Effects Estimates', pch=19)
abline(h=0, lty=2)
abline(v=0, lty=2)


# Save fitted model objects for reproducibility and subsequent analyses
save(jm.rondeau, jm.ng, file = "results/JM_parametric.Rdata")

#-----------------------------------------------------------------------
# Manuscript outputs
#-----------------------------------------------------------------------

#---------#
# Table 1 #
#---------#
rand.eff.ng = data.frame(
  Model = "Ng et al",
  Parameter = c("theta_u2","theta_v2","rho"),
  Estimate = jm.ng$Frailty[,1],
  StdDev = jm.ng$Frailty[,2],
  pvalue = jm.ng$Frailty[,3]
)
rownames(rand.eff.ng) <- NULL
table1 <- rbind(rand.eff.rondeau,rand.eff.ng)

write.csv2(table1, "results/tables/table_1.csv", row.names = FALSE)


#----------#
# Figure 6 #
#----------#
pdf("results/figures/figure_6.pdf", width = 12, height = 5)
frailty.comparisons(jm.rondeau, jm.ng, result = gauss3, 
                    colors = c('#CC3300','#FF6600','#FF9933'))
dev.off() 

