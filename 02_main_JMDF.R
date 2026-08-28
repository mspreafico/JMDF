########################################################################
# Joint models with discrete non-parametric frailty
########################################################################
# This script fits the JMDF model to the example dataset using Gaussian
# and Uniform initializations and different values of the distance
# threshold L. The fitted models are used to reproduce the results, Table 7,
# Figures 3 and 4 presented in the manuscript.
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
load("data/fake_dataRD.Rdata")
# The example dataset contains recurrent-event and terminal-event data.

## load functions
source('functions/data_format.R')
source('functions/JMdiscfrail.R')
source('functions/jmdf_plots.R')


# Format recurrent-event and terminal-event data for JMDF estimation
dataR = formatting.data(dataR)
dataD = formatting.data(dataD)

###########################
# Gaussian Initialization #
###########################
# Fit JMDF using Gaussian initialization for two different values of
# the distance threshold L.

SigmaG =  matrix(c(2*0.12,0.0,0.0,2*1.4),nrow = 2, ncol = 2)
muG = c(0,0)

# L = 1
gauss3 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = FALSE,
                     distance = "euclidean",
                     Sigma = SigmaG,
                     mu = muG,
                     M = 1000, L = 1,
                     max.it = 10,  toll = 1e-3)

# Fixed effects - Recurrent (betas)
gauss3$modelR
# Fixed effects - Terminal (gammas)
gauss3$modelD
# Mass points
gauss3$K
print(plot.masses(gauss3, colors=c('#CC3300','#FF6600','#FF9933')))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(gauss3, colors=c('#CC3300','#FF6600','#FF9933')))

# L = 2
gauss2 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = FALSE,
                     distance = "euclidean",
                     Sigma = SigmaG,
                     mu = muG,
                     M = 1000, L = 2,
                     max.it = 10,  toll = 1e-3)
# Fixed effects - Recurrent (betas)
gauss2$modelR
# Fixed effects - Terminal (gammas)
gauss2$modelD
# Mass points
gauss2$K
print(plot.masses(gauss2, colors=c('#3399FF','#0066CC')))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(gauss2, colors=c('#3399FF','#0066CC')))


##########################
# Uniform Initialization #
##########################
# Fit JMDF using Uniform initialization for two different values of
# the distance threshold L.

Sigma = matrix(c(0.12,0.0,0.0,1.4), nrow = 2, ncol = 2)
mu = c(0,0)
ulim = c(mu[1]-6*sqrt(Sigma[1,1]), mu[1]+6*sqrt(Sigma[1,1]))
vlim = c(mu[2]-6*sqrt(Sigma[2,2]), mu[2]+6*sqrt(Sigma[2,2]))

# L = 1.5
unif3 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = TRUE,
                     distance = "euclidean",
                     ulim.unif = ulim,
                     vlim.unif = vlim,
                     M = 1000, L = 1.5,
                     max.it = 10,  toll = 1e-3)
# Fixed effects - Recurrent (betas)
unif3$modelR
# Fixed effects - Terminal (gammas)
unif3$modelD
# Mass points
unif3$K
print(plot.masses(unif3, colors=c('#FF6699','#FF33CC','#CC0066')))
# --> Mass point P3 is not visible is its weight is very small
# Frailty-stratified baseline survival curves
print(figure.survival.curves(unif3, colors=c('#FF6699','#FF33CC','#CC0066')))

# L = 2
unif2 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                    dataD, formulaD = '~ sex + age + ncom + adherent',
                    init.unif = TRUE,
                    distance = "euclidean",
                    ulim.unif = ulim,
                    vlim.unif = vlim,
                    M = 1000, L = 2,
                    max.it = 10,  toll = 1e-3)
# Fixed effects - Recurrent (betas)
unif2$modelR
# Fixed effects - Terminal (gammas)
unif2$modelD
# Mass points
unif2$K
print(plot.masses(unif2, colors=c('#0066CC','#3399FF')))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(unif2, colors=c('#0066CC','#3399FF')))


# Save fitted model objects for reproducibility and subsequent analyses
save(gauss2, gauss3, unif2, unif3, file = "results/JMDF_gauss_unif.Rdata")


#-----------------------------------------------------------------------
# Manuscript outputs
#-----------------------------------------------------------------------

#---------#
# Table 7 #
#---------#
data_masses = create.data.masses(gauss3, gauss2, unif3, unif2)
data_masses[, P := seq_len(.N), by = .(init, K)]
appendixB = data_masses[,.(init,K,P,u,SE_u,v,SE_v,w)]

write.csv2(appendixB, "results/tables/table_7.csv", row.names = FALSE)


#----------#
# Figure 3 #
#----------#
figure3 = figure.masses.combined(gauss3, gauss2, unif3, unif2)
ggsave("results/figures/figure_3.pdf", plot = figure3, width = 10, height = 5)


#----------#
# Figure 4 #
#----------#
pdf("results/figures/figure_4_unif.pdf", width = 12, height = 5)
figure.survival.curves(unif3, colors=c('#CC0066', '#FF33CC', '#FF6699'))
dev.off()

pdf("results/figures/figure_4_gauss.pdf", width = 12, height = 5)
figure.survival.curves(gauss3, colors=c('#CC3300','#FF6600','#FF9933'))
dev.off()


