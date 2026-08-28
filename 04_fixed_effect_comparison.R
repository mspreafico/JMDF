########################################################################
# Fixed effects comparisons across different models
########################################################################
# Main script for comparing the estimated fixed effects across the
# joint models considered in the application.
#
# The script loads the fitted model objects produced by the main
# application scripts and computes hazard ratios (HRs) with 95% CIs
# for the recurrent and terminal events.
#
# The resulting comparison is reported in Figure 5 of the manuscript.
# The figure is saved as a PDF file in 'results/figures/' with a
# manuscript-matching filename.
########################################################################


## clear workspace
rm(list=ls())

## set working directory
setwd("~/github/JMDF")

## load packages
library(ggplot2)
library(ggpubr)
library(data.table)

# Load results
load("results/JM_parametric.Rdata")
load("results/JMDF_gauss_unif.Rdata")

########################################################################
# Functions
########################################################################
# Extract hazard ratios and 95% confidence intervals for the fixed
# effects from the supported model classes (JMDF, Ng et al., and
# Rondeau et al.).
HR.fixed.effect <- function(model, namesR=NULL, namesT=NULL, model.type='JMDF'){
  
  if(is.null(model.type)){
    stop('Please specify model.type')
  }
  
  if(model.type=='JMDF'){
    betas = model[[1]]$coefficients
    betas_se = sqrt(diag(model[[1]]$var))
    gammas = model[[2]]$coefficients
    gammas_se = sqrt(diag(model[[2]]$var))
    if(is.null(namesR) | length(namesR)!=length(model[[1]]$coefficients)){
      namesR = names(model[[1]]$coefficients)
    }
    if(is.null(namesT) | length(namesT)!=length(model[[2]]$coefficients)){
      namesT = names(model[[2]]$coefficients)
    }
  }else if(model.type=='ng'){
    betas = model$Recurrent[,1]
    betas_se = model$Recurrent[,2]
    gammas = model$Death[,1]
    gammas_se = model$Death[,2]
    if(is.null(namesR) | length(namesR)!=dim(model$Recurrent)[1]){
      namesR = row.names(model$Recurrent)
    }
    if(is.null(namesT) | length(namesT)!=dim(model$Death)[1]){
      namesT = row.names(model$Death)
    }
  }else if(model.type=='rondeau'){
    betas = model$coef[1:model$nvar[1]]
    betas_se = sqrt(diag(model$varH))[2+(1:model$nvar[1])]
    gammas = model$coef[(model$nvar[1]+1):(model$nvar[1]+model$nvar[2])]
    gammas_se = sqrt(diag(model$varH))[2+(model$nvar[1]+1):(model$nvar[1]+model$nvar[2])]
    if(is.null(namesR) | length(namesR)!=length(names(model$coef[1:model$nvar[1]]))){
      namesR = names(model$coef[1:model$nvar[1]])
    }
    if(is.null(namesT) | length(namesT)!=length(names(model$coef[(model$nvar[1]+1):(model$nvar[1]+model$nvar[2])]))){
      namesT = names(model$coef[(model$nvar[1]+1):(model$nvar[1]+model$nvar[2])])
    }
  }
  
  HR_fixed_R = cbind.data.frame('Variable' = namesR,
                                'HR' = exp(betas),
                                'L95' = exp(betas-1.96*betas_se),
                                'U95' = exp(betas+1.96*betas_se))
  row.names(HR_fixed_R) = c(1:dim(HR_fixed_R)[1])
  
  HR_fixed_T = cbind.data.frame('Variable' = namesT,
                                'HR' = exp(gammas),
                                'L95' = exp(gammas-1.96*gammas_se),
                                'U95' = exp(gammas+1.96*gammas_se))
  row.names(HR_fixed_T) = c(1:dim(HR_fixed_T)[1])
  
  return(list('recurrent'=HR_fixed_R,'terminal'=HR_fixed_T))
}

# Create the model comparison plot for a selected fixed-effect variable
HRplot = function(data, limits, breaks, order = c('G2', 'G3', 'U2', 'U3','R', 'N')){
  
  data$Models <- factor(data$Models, levels=order)
  
  p = ggplot(data, aes(x=Models, y=HR, colour=Models)) + 
    geom_hline(yintercept = 1, linetype = "dashed", colour='gray20') +
    geom_errorbar(aes(ymin=L95, ymax=U95), width=.1) +
    geom_line() +
    geom_point(size=2) +
    #scale_colour_brewer(palette='Dark2',
    scale_colour_manual(values = c('#33CC99','#006633','#3366FF','#000099',
                                   '#996633','#FF6666'),
                        labels=c('G2: Gaussian, K=2', 'G3: Gaussian, K=3',
                                 'U2: Uniform, K=2', 'U3: Uniform, K=3',
                                 'R: Rondeau and others (2007)', 'N: Ng and others (2023)')) +
    scale_y_continuous(limits = limits, breaks = breaks) +
    xlab(NULL) +
    ylab('HRs with 95% CIs') +
    theme_bw() +
    theme(legend.title = element_text(face="bold",size=rel(1.1)), legend.text = element_text(size=rel(1.1)))
  return(p)
}


########################################################################
# Fixed effects comparisons
########################################################################

# Extract and combine the fixed-effect estimates from all fitted models
# considered in the application
var.names = c('Male', 'Age/10', 'Comorbidity', 'Adherent')
hr.rond = HR.fixed.effect(jm.rondeau, namesR = var.names, namesT = var.names,
                          model.type = 'rondeau')
hr.ng = HR.fixed.effect(jm.ng, namesR = var.names, namesT = var.names, 
                        model.type = 'ng')
hr.G3 = HR.fixed.effect(gauss3, namesR = var.names, namesT = var.names)
hr.G2 = HR.fixed.effect(gauss2, namesR = var.names, namesT = var.names)
hr.U3 = HR.fixed.effect(unif3, namesR = var.names, namesT = var.names)
hr.U2 = HR.fixed.effect(unif2, namesR = var.names, namesT = var.names)

# Create data
ncovR = 4
ncovT = 4
data.fixed = list()
# Recurrent fixed effects (HRs)
data.fixed$recurrent = data.table(rbind.data.frame(
  cbind.data.frame('Models' = rep('G2',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.G2$recurrent),
  cbind.data.frame('Models' = rep('G3',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.G3$recurrent),
  cbind.data.frame('Models' = rep('U2',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.U2$recurrent),
  cbind.data.frame('Models' = rep('U3',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.U3$recurrent),
  cbind.data.frame('Models' = rep('R',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.rond$recurrent),
  cbind.data.frame('Models' = rep('N',ncovR),
                   'Event' = rep('Recurrent',ncovR),
                   hr.ng$recurrent)
))
row.names(data.fixed$recurrent) = 1:dim(data.fixed$recurrent)[1]

# Terminal fixed effects (HRs)
data.fixed$terminal = data.table(rbind.data.frame(
  cbind.data.frame('Models' = rep('G2',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.G2$terminal),
  cbind.data.frame('Models' = rep('G3',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.G3$terminal),
  cbind.data.frame('Models' = rep('U2',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.U2$terminal),
  cbind.data.frame('Models' = rep('U3',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.U3$terminal),
  cbind.data.frame('Models' = rep('R',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.rond$terminal),
  cbind.data.frame('Models' = rep('N',ncovT),
                   'Event' = rep('Terminal',ncovT),
                   hr.ng$terminal)
))
row.names(data.fixed$terminal) = 1:dim(data.fixed$terminal)[1]

#----------#
# Figure 5 #
#----------#

# Panels
limT=c(0.1,2.7)
breakT = seq(0.1,2.7,by=0.2)
limR=c(0.5,1.3)
breakR = seq(0.5,1.3,by=0.2)
pR1 = HRplot(data = data.fixed$recurrent[Variable=='Male'],
             limits = limR, breaks = breakR)
pT1 = HRplot(data = data.fixed$terminal[Variable=='Male'],
             limits = limT, breaks = breakT)
pR2 = HRplot(data = data.fixed$recurrent[Variable=='Age/10'],
             limits = limR, breaks = breakR)
pT2 = HRplot(data = data.fixed$terminal[Variable=='Age/10'],
             limits = limT, breaks = breakT)
pR3 = HRplot(data = data.fixed$recurrent[Variable=='Comorbidity'],
             limits = limR, breaks = breakR)
pT3 = HRplot(data = data.fixed$terminal[Variable=='Comorbidity'],
             limits = limT, breaks = breakT)
pR4 = HRplot(data = data.fixed$recurrent[Variable=='Adherent'],
             limits = limR, breaks = breakR)
pT4 = HRplot(data = data.fixed$terminal[Variable=='Adherent'],
             limits = limT, breaks = breakT)

# Labels
row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Recurrent Event", angle = -90, size=6) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Terminal Event", angle = -90, size=6) + theme_void() 
col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Sex (Male)", angle = 0, size=5) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Age (10-year)", angle = 0, size=5) + theme_void() 
col3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Comorbidity", angle = 0, size=5) + theme_void() 
col4 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Adherent (Yes)", angle = 0, size=5) + theme_void() 

figure5 = ggarrange(col1, col2, col3, col4, NULL,
                   pR1 + rremove("ylab"), pR2 + rremove("ylab"), pR3 + rremove("ylab"), pR4 + rremove("ylab"), row1,
                   pT1 + rremove("ylab"), pT2 + rremove("ylab"), pT3 + rremove("ylab"), pT4 + rremove("ylab"), row2,
                   nrow=3, ncol=5, 
                   heights = c(0.04,0.48,0.48),
                   widths = c(0.24,0.24,0.24,0.24,0.04),
                   common.legend = TRUE, legend = "bottom")
figure5 = annotate_figure(figure5,
                left = text_grob("HRs with 95% CIs", size = 12, rot = 90))
ggsave("results/figures/figure_5.pdf", plot = figure5, width = 13, height = 7)


