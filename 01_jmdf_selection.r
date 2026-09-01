########################################################################
# Joint models with discrete non-parametric frailty - Model selection ##
########################################################################
# This script is used to perform JMDF model selection in terms of the threshold L.
# The script fits the JMDF model to the example dataset using Gaussian
# and Uniform initializations and different values of the distance
# threshold L, ranging from 0.1 to 3, across 12 different runs. 
# The AIC and the estimated number of mass points of JMDF models are computed 
# for each run and results are saved and plotted to reproduce  
# Figure 2 presented in the manuscript.
#
# All numerical results and manuscript figures/tables are saved in the
# 'results/' directory.
########################################################################

## clear workspace
rm(list=ls())

## load libraries
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)

## set working directory
setwd("~/github/JMDF")

## load data
load("data/fake_dataRD.Rdata")
# The example dataset contains recurrent-event and terminal-event data.

## load functions
source('functions/data_format.R')
source('functions/JMdiscfrail.R')


# Format recurrent-event and terminal-event data for JMDF estimation
dataR = formatting.data(dataR)
dataD = formatting.data(dataD)

# Auxiliary function for running JMDF with different L and initializations 
fit_jmdf <- function(init, L,seed_jmdf) {
  if (init == 'gauss') {
    SigmaG <- matrix(c(2*.12,0,0,2*1.4),2,2)
    muG = c(0,0)
    JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                dataD, formulaD = '~ sex + age + ncom + adherent',
                init.unif = FALSE,
                distance = "euclidean",
                Sigma = SigmaG,
                mu = muG,
                M = 1000, L = L,
                max.it = 10,  toll = 1e-3, seed= seed_jmdf)
  } else {
    Sigma = matrix(c(0.12,0.0,0.0,1.4), nrow = 2, ncol = 2)
    mu = c(0,0)
    ulim = c(mu[1]-6*sqrt(Sigma[1,1]), mu[1]+6*sqrt(Sigma[1,1]))
    vlim = c(mu[2]-6*sqrt(Sigma[2,2]), mu[2]+6*sqrt(Sigma[2,2]))
    JMdiscfrail(dataR, formulaR = '~ sex + age + ncom + adherent',
                dataD, formulaD = '~ sex + age + ncom + adherent',
                init.unif = TRUE,
                distance = "euclidean",
                ulim.unif = ulim,
                vlim.unif = vlim,
                M = 1000, L = L,
                max.it = 10,  toll = 1e-3, seed = seed_jmdf)
  }
}


## Model-selection route across different L thresholds and 12 runs

# Set the 12 seeds
seeds <- c(210197:210208)

# Set the candidates for L
candidates <- seq(0.1,3.0,by=.1)

colgauss = c(brewer.pal(9, name="PuBu")[6:9],brewer.pal(9, name="YlGn")[-1])
colunif = c(brewer.pal(9, name="PuBuGn")[6:9],brewer.pal(9, name="YlGnBu")[-1])


#fits_un <- list() 
#fits_gauss <- list()
model_selection_un <- data.frame()
model_selection_gauss <- data.frame()

###########################
# Gaussian Initialization #
###########################
init = 'gauss'
for(seed in seeds){
  for (L in candidates) {
    message('Pseudo-data JMDF ',init,', L=',round(L,3))
    seed_jmdf <- round(seed+L*10)
    z <- fit_jmdf(init,L, seed_jmdf)
    #key <- paste(seed,format(L,digits=8),sep='_')
    #fits_gauss[[key]] <- z
    model_selection_gauss <- rbind(model_selection_gauss, data.frame(seed=as.factor(seed),L=L,
                                                                     AIC=z$AIC, 
                                                                     LogL = z$LogL,
                                                                     CLogL = z$classLogL,
                                                                     K=z$K))
  }}

save(model_selection_gauss, file = "results/JMDF_model_selection_gauss.Rdata")

# load(results/JMDF_model_selection_gauss.Rdata)
p1 <- ggplot(model_selection_gauss, aes(x=L, y=AIC, color=seed)) +  geom_line()+
  labs(title = "AIC vs Distance L across runs",
       x = "Distance L",
       y = "AIC") +
  scale_y_continuous(limits = c(18600,19400)) +
  scale_colour_manual(values = colgauss) +
  theme_light() +
  theme(legend.position='none', 
        axis.text=element_text(size=rel(1.1)),
        axis.title=element_text(size=rel(1.2)),
        plot.title = element_text(face="bold", size=rel(1.3)), 
        legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))

p1
p2 <- ggplot(model_selection_gauss, aes(x=L, y=K, color=seed))+
  geom_line()+
  scale_y_continuous(limits = c(1,17), breaks=seq(1,17,2)) +
  scale_colour_manual(values = colgauss) +
  labs(title = "Number of masses vs Distance L across runs",
       x = "Distance L",
       y = "Cardinality") +
  theme_light() +
  theme(legend.position='none', 
        axis.text=element_text(size=rel(1.1)),
        axis.title=element_text(size=rel(1.2)),
        plot.title = element_text(face="bold", size=rel(1.3)), 
        legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))

p2

##########################
# Uniform Initialization #
##########################
init = 'unif'
for(seed in seeds){
for (L in candidates) {
  message('Pseudo-data JMDF ',init,', L=',round(L,3))
  seed_jmdf <- round(seed+L*10)
  z <- fit_jmdf(init,L, seed_jmdf)
  #key <- paste(seed,format(L,digits=8),sep='_')
  #fits_un[[key]] <- z
  model_selection_un <- rbind(model_selection_un, data.frame(seed=as.factor(seed),L=L,
                                                       AIC=z$AIC, 
                                                       LogL = z$LogL,
                                                       CLogL = z$classLogL,
                                                       K=z$K))
}}

save(model_selection_un, file = "results/JMDF_model_selection_unif.Rdata")

#load("results/JMDF_model_selection_unif.Rdata")

p3 <- ggplot(model_selection_un, aes(x=L, y=AIC, color=seed)) +  geom_line()+
  labs(title = "AIC vs Distance L across runs",
       x = "Distance L",
       y = "AIC") +
  scale_y_continuous(limits = c(18600,19400)) +
  scale_colour_manual(values = colunif) +
  theme_light() +
  theme(legend.position='none', 
        axis.text=element_text(size=rel(1.1)),
        axis.title=element_text(size=rel(1.2)),
        plot.title = element_text(face="bold", size=rel(1.3)), 
        legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))

p3

p4 <- ggplot(model_selection_un, aes(x=L, y=K, color=seed))+
  geom_line()+
  scale_y_continuous(limits = c(1,17), breaks=seq(1,17,2)) +
  scale_colour_manual(values = colunif) +
  labs(title = "Number of masses vs Distance L across runs",
       x = "Distance L",
       y = "Cardinality") +
  theme_light() +
  theme(legend.position='none', 
        axis.text=element_text(size=rel(1.1)),
        axis.title=element_text(size=rel(1.2)),
        plot.title = element_text(face="bold", size=rel(1.3)), 
        legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))

p4


## Code to reproduce Figure 2
col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Gaussian Initialization", 
                            angle = 0, size=6, fontface = "bold") + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Uniform Initialization", 
                            angle = 0, size=6, fontface = "bold") + theme_void() 

figure2<- ggarrange(col1, col2,
          ggarrange(p1, p2, align = "v", ncol=1),
          ggarrange(p3, p4, align = "v", ncol=1),
          nrow=2, ncol=2, heights = c(0.06,0.94))

# Save Figure 2
ggsave("results/figures/figure_2.pdf", plot = figure2, width = 10, height = 5)

