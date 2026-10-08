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

# L = 1.25
gauss3 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = FALSE,
                     distance = "euclidean",
                     Sigma = SigmaG,
                     mu = muG,
                     M = 1000, L = 1.25,
                     max.it = 10,  toll = 1e-3, seed = 4)

# Fixed effects - Recurrent (betas)
gauss3$modelR
# Fixed effects - Terminal (gammas)
gauss3$modelD
# Mass points
gauss3$K
col_gauss3 = c('#FF9933','#FF6600','#CC3300')[order(gauss3$P[,1])]
print(plot.masses(gauss3, colors = col_gauss3))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(gauss3, colors = col_gauss3))

# L = 2
gauss2 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = FALSE,
                     distance = "euclidean",
                     Sigma = SigmaG,
                     mu = muG,
                     M = 1000, L = 2,
                     max.it = 10,  toll = 1e-3, seed = 3)
# Fixed effects - Recurrent (betas)
gauss2$modelR
# Fixed effects - Terminal (gammas)
gauss2$modelD
# Mass points
gauss2$K
col_gauss2 = c('#0066CC','#3399FF')[order(gauss2$P[,1])]
print(plot.masses(gauss2, colors = col_gauss2))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(gauss2, colors = col_gauss2))


##########################
# Uniform Initialization #
##########################
# Fit JMDF using Uniform initialization for two different values of
# the distance threshold L.

Sigma = matrix(c(0.12,0.0,0.0,1.4), nrow = 2, ncol = 2)
mu = c(0,0)
ulim = c(mu[1]-6*sqrt(Sigma[1,1]), mu[1]+6*sqrt(Sigma[1,1]))
vlim = c(mu[2]-6*sqrt(Sigma[2,2]), mu[2]+6*sqrt(Sigma[2,2]))

# L = 1.25
unif3 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                     dataD, formulaD = '~ sex + age + ncom + adherent',
                     init.unif = TRUE,
                     distance = "euclidean",
                     ulim.unif = ulim,
                     vlim.unif = vlim,
                     M = 1000, L = 1.25,
                     max.it = 10,  toll = 1e-3, seed = 2)
# Fixed effects - Recurrent (betas)
unif3$modelR
# Fixed effects - Terminal (gammas)
unif3$modelD
# Mass points
unif3$K
col_unif3 = c('#FF6699','#FF33CC','#CC0066')[order(unif3$P[,1])]
print(plot.masses(unif3, colors = col_unif3))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(unif3, colors = col_unif3))

# L = 2
unif2 = JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                    dataD, formulaD = '~ sex + age + ncom + adherent',
                    init.unif = TRUE,
                    distance = "euclidean",
                    ulim.unif = ulim,
                    vlim.unif = vlim,
                    M = 1000, L = 2,
                    max.it = 10,  toll = 1e-3, seed = 1)
# Fixed effects - Recurrent (betas)
unif2$modelR
# Fixed effects - Terminal (gammas)
unif2$modelD
# Mass points
unif2$K
col_unif2 = c('#0066CC','#3399FF')[order(unif2$P[,1])]
print(plot.masses(unif2, colors = col_unif2))
# Frailty-stratified baseline survival curves
print(figure.survival.curves(unif2, colors = col_unif2))


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

appendixB[, c("u", "v", "w") := lapply(.SD, formatC, format="f", digits = 3),
          .SDcols = c("u", "v", "w")]
appendixB[, c("SE_u", "SE_v") := lapply(.SD, formatC, format="f", digits = 4),
          .SDcols = c("SE_u", "SE_v")]

write.csv2(appendixB, "results/tables/table_7.csv", row.names = FALSE)


#----------#
# Figure 3 #
#----------#
figure3 = figure.masses.combined(gauss3, gauss2, unif3, unif2,
                                 col.g3 = col_gauss3, 
                                 col.g2 = col_gauss2,
                                 col.u3 = col_unif3, 
                                 col.u2 = col_unif2)
ggsave("results/figures/figure_3.pdf", plot = figure3, width = 10, height = 5)


#----------#
# Figure 4 #
#----------#
pdf("results/figures/figure_4_unif.pdf", width = 12, height = 5)
figure.survival.curves(unif3, colors = col_unif3)
dev.off()

pdf("results/figures/figure_4_gauss.pdf", width = 12, height = 5)
figure.survival.curves(gauss3, colors = col_gauss3)
dev.off()


