################################################################################
# S2. Simulation study: additional results under setting B
################################################################################
# This file reproduces Tables and Figures in Supplementary Material S2.
#-------------------------------------------------------------------------------
# Before running this file:
# 1.  First modify the value of `setting` in "01_data_generation.R" to 'B'. 
#     After making this change, run the file to generate data under setting B 
#     (yielding approximately 50% death and 50% administrative censoring).  
# 2.  Then, modify the value of `setting` in "02_fit_jm_bivnorm_ng.R", 
#     "02_fit_jmdf_gauss.R", and "02_fit_jmdf_unif.R" to 'B' and run those files 
#     to obtain the results under setting B.  
#     For JMDF, select the desired initialization (i or ii) and the threshold L 
#     by modifying the value of `folder` ('I_L1', 'II_L1', 'II_L15', or 'II_L2'). 
#     You must run the JMDF files separately for each initialization to get results 
#     corresponding to all choices.
################################################################################


# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
file.path(getwd(),'simstudy')

library(dplyr)
library(gt)
gt.tables <- function(df, title){
  df %>% gt() %>% fmt_number(columns = -1, decimals = 3) %>%
    tab_header(title)
}

source('sim_functions/perf_measures_ng.R')
source('sim_functions/perf_measures_jmdf.R')

# TRUE VALUES
beta.true=c(-0.6, 0.8)
gamma.true=c(-0.8, 0.5)
true.ng.param = matrix(c(0.8,0.8,0.8,
                         0.8,0.8,0.2,
                         0.03,0.02,0.2,
                         rep(NA,3*6),
                         1,2,0.8,
                         NA,NA,NA), nrow=11, ncol=3, byrow=T)
true.K.masses = c('-','-','-',2,2,3,3,4,3,'-',3)
# True masses
Pu.true = list(); Pv.true= list(); weights.true = list()
Pu.true[[4]] = c(-1.2,1.2); Pv.true[[4]] = c(-1,1); weights.true[[4]] = c(0.5,0.5)
Pu.true[[5]] = c(-1.2,0.51); Pv.true[[5]] = c(-1,0.43); weights.true[[5]] = c(0.3,0.7)
Pu.true[[6]] = c(-1.2,0.68,1.2); Pv.true[[6]] = c(-1,0.29,1.2); weights.true[[6]] = c(0.45,0.23,0.32)
Pu.true[[7]] = c(-1.2,2.2,-0.71); Pv.true[[7]] = c(-1,0.7,0.98); weights.true[[7]] = c(0.45,0.32,0.23)
Pu.true[[8]] = c(-2,2,-2,2); Pv.true[[8]] = c(-2,-2,2,2); weights.true[[8]] = c(0.25,0.25,0.25,0.25)
Pu.true[[9]] = c(-1.5,0.75,2.25); Pv.true[[9]] = c(0,0,0); weights.true[[9]] = c(0.4,0.5,0.1)
Pu.true[[11]] = c(-1,0,1); Pv.true[[11]] = c(-2,0,2); weights.true[[11]] = c(1/3,1/3,1/3)

# Data Setting
setting = 'B'


# S2.1 Results from joint model by Ng et al
#-------------------------------------------------------------------------------
jfm.ng = perf.measures.ng.all.scenarios(setting, scenarios = 1:9,
                                        beta = beta.true, gamma = gamma.true, 
                                        coef.ng.true = true.ng.param, N = 500)


gt.tables(jfm.ng, title = "Table S2.1: Ng et al.")


# S2.2 Results from JMDF with Gaussian initializations
#-------------------------------------------------------------------------------
gauss.i.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'I_L15', init = 'gauss', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)

gauss.ii.1 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L1', init = 'gauss', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)

gauss.ii.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L15', init = 'gauss', 
                                               scenarios = 1:9, 
                                               beta = beta.true, gamma = gamma.true,
                                               true.K = true.K.masses)

gauss.ii.2 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L2', init = 'gauss', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)
# Tables
gt.tables(gauss.i.15$perf.coefs, title = "Table S2.2: Gaussian (i), L = 1.5")
gt.tables(gauss.ii.1$perf.coefs, title = "Table S2.3: Gaussian (ii), L = 1")
gt.tables(gauss.ii.1$perf.coefs, title = "Table S2.4: Gaussian (ii), L = 1.5")
gt.tables(gauss.ii.2$perf.coefs, title = "Table S2.5: Gaussian (ii), L = 2")
gt.tables(gauss.i.15$K.est, title = "Table S2.6: Gaussian (i), L = 1.5")
gt.tables(gauss.ii.1$K.est, title = "Table S2.6: Gaussian (ii), L = 1")
gt.tables(gauss.ii.15$K.est, title = "Table S2.6: Gaussian (ii), L = 1.5")
gt.tables(gauss.ii.2$K.est, title = "Table S2.6: Gaussian (ii), L = 2")


# Figure S2.1
correct.k.boxplots(setting, folder = 'I_L15', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.i.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure S2.2
correct.k.boxplots(setting, folder = 'II_L1', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.1,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure S2.3
correct.k.boxplots(setting, folder = 'II_L15', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.1,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure S2.4
correct.k.boxplots(setting, folder = 'II_L2', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.2,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)



# S2.3 Results from JMDF with Uniform initializations
#-------------------------------------------------------------------------------
unif.i.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'I_L15', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

unif.ii.1 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L1', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

unif.ii.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L15', init = 'unif', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)

unif.ii.2 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L2', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

gt.tables(unif.i.15$perf.coefs, title = "Table S2.7: Uniform (i), L = 1.5")
gt.tables(unif.ii.1$perf.coefs, title = "Table S2.8: Uniform (ii), L = 1")
gt.tables(unif.ii.15$perf.coefs, title = "Table S2.8: Uniform (ii), L = 1.5")
gt.tables(unif.ii.2$perf.coefs, title = "Table S2.10: Uniform (ii), L = 2")
gt.tables(unif.i.15$K.est, title = "Table S2.11: Uniform (i), L = 1.5")
gt.tables(unif.ii.1$K.est, title = "Table S2.11: Uniform (ii), L = 1")
gt.tables(unif.ii.15$K.est, title = "Table S2.11: Uniform (ii), L = 1.5")
gt.tables(unif.ii.2$K.est, title = "Table S2.11: Uniform (ii), L = 2")


# Figure S2.5
correct.k.boxplots(setting, folder = 'I_L15', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.i.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure S2.6
correct.k.boxplots(setting, folder = 'II_L1', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.ii.1,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure S2.7
correct.k.boxplots(setting, folder = 'II_L15', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.ii.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)

# Figure S2.8
correct.k.boxplots(setting, folder = 'II_L2', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.ii.2,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)





