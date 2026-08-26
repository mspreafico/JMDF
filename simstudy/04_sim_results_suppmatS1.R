################################################################################
# S1. Simulation study: additional results under setting A
################################################################################
# This file reproduces Tables and Figures in Supplementary Material S1.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

library(dplyr)
library(gt)
gt.tables <- function(df, title){
  df %>% gt() %>% fmt_number(columns = -1, decimals = 3) %>%
    tab_header(title)
}

source('sim_functions/perf_measures_jmdf.R')

# TRUE VALUES
beta.true=c(-0.6, 0.8)
gamma.true=c(-0.8, 0.5)
true.ng.param = matrix(c(0.8,0.8,0.8,
                         0.8,0.8,0.2,
                         0.03,0.02,0.2,
                         rep(NA,3*6)), nrow=9, ncol=3, byrow=T)
true.K.masses = c('-','-','-',2,2,3,3,4,3)
# True masses
Pu.true = list(); Pv.true= list(); weights.true = list()
Pu.true[[4]] = c(-1.2,1.2); Pv.true[[4]] = c(-1,1); weights.true[[4]] = c(0.5,0.5)
Pu.true[[5]] = c(-1.2,0.51); Pv.true[[5]] = c(-1,0.43); weights.true[[5]] = c(0.3,0.7)
Pu.true[[6]] = c(-1.2,0.68,1.2); Pv.true[[6]] = c(-1,0.29,1.2); weights.true[[6]] = c(0.45,0.23,0.32)
Pu.true[[7]] = c(-1.2,2.2,-0.71); Pv.true[[7]] = c(-1,0.7,0.98); weights.true[[7]] = c(0.45,0.32,0.23)
Pu.true[[8]] = c(-2,2,-2,2); Pv.true[[8]] = c(-2,-2,2,2); weights.true[[8]] = c(0.25,0.25,0.25,0.25)
Pu.true[[9]] = c(-1.5,0.75,2.25); Pv.true[[9]] = c(0,0,0); weights.true[[9]] = c(0.4,0.5,0.1)

# Data Setting
setting = 'A'
res.tab = paste0('sim_results_',setting,'/tables')


# S1.1 Results from JMDF with Gaussian initializations
#-------------------------------------------------------------------------------
gauss.i.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'I_L15', init = 'gauss', 
                                                scenarios = 1:9, 
                                                beta = beta.true, gamma = gamma.true,
                                                true.K = true.K.masses)

gauss.ii.1 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L1', init = 'gauss', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)

gauss.ii.2 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L2', init = 'gauss', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)
# Tables
gt.tables(gauss.i.15$perf.coefs, title = "Table S1.1: Gaussian (i), L = 1.5")
write.csv2(gauss.i.15$perf.coefs,file = paste0(res.tab,'/table_S1_1.csv'), row.names = FALSE)

gt.tables(gauss.ii.1$perf.coefs, title = "Table S1.2: Gaussian (ii), L = 1")
write.csv2(gauss.ii.1$perf.coefs,file = paste0(res.tab,'/table_S1_2.csv'), row.names = FALSE)

gt.tables(gauss.ii.2$perf.coefs, title = "Table S1.3: Gaussian (ii), L = 2")
write.csv2(gauss.ii.2$perf.coefs,file = paste0(res.tab,'/table_S1_3.csv'), row.names = FALSE)

gt.tables(gauss.i.15$K.est, title = "Table S1.4: Gaussian (i), L = 1.5")
write.csv2(gauss.i.15$K.est,file = paste0(res.tab,'/table_S1_4_I_L15.csv'), row.names = FALSE)

gt.tables(gauss.ii.1$K.est, title = "Table S1.4: Gaussian (ii), L = 1")
write.csv2(gauss.ii.1$K.est,file = paste0(res.tab,'/table_S1_4_II_L1.csv'), row.names = FALSE)

gt.tables(gauss.ii.2$K.est, title = "Table S1.4: Gaussian (ii), L = 2")
write.csv2(gauss.ii.2$K.est,file = paste0(res.tab,'/table_S1_4_II_L2.csv'), row.names = FALSE)



# Figure S1.1
correct.k.boxplots(setting, folder = 'I_L15', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.i.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_1.pdf")
# Figure S1.2
correct.k.boxplots(setting, folder = 'II_L1', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.1,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_2.pdf")
# Figure S1.3
correct.k.boxplots(setting, folder = 'II_L2', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.2,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_3.pdf")


# S1.2 Results from JMDF with Uniform initializations
#-------------------------------------------------------------------------------
unif.i.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'I_L15', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

unif.ii.1 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L1', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

unif.ii.2 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L2', init = 'unif', 
                                             scenarios = 1:9, 
                                             beta = beta.true, gamma = gamma.true,
                                             true.K = true.K.masses)

gt.tables(unif.i.15$perf.coefs, title = "Table S1.5: Uniform (i), L = 1.5")
write.csv2(unif.i.15$perf.coefs, file = paste0(res.tab,'/table_S1_5.csv'), row.names = FALSE)

gt.tables(unif.ii.1$perf.coefs, title = "Table S1.6: Uniform (ii), L = 1")
write.csv2(unif.ii.1$perf.coefs, file = paste0(res.tab,'/table_S1_6.csv'), row.names = FALSE)

gt.tables(unif.ii.2$perf.coefs, title = "Table S1.7: Uniform (ii), L = 2")
write.csv2(unif.ii.2$perf.coefs, file = paste0(res.tab,'/table_S1_7.csv'), row.names = FALSE)

gt.tables(unif.i.15$K.est, title = "Table S1.8: Uniform (i), L = 1.5")
write.csv2(unif.i.15$K.est, file = paste0(res.tab,'/table_S1_8_I_L15.csv'), row.names = FALSE)

gt.tables(unif.ii.1$K.est, title = "Table S1.8: Uniform (ii), L = 1")
write.csv2(unif.ii.1$K.est, file = paste0(res.tab,'/table_S1_8_II_L1.csv'), row.names = FALSE)

gt.tables(unif.ii.2$K.est, title = "Table S1.8: Uniform (ii), L = 2")
write.csv2(unif.ii.2$K.est, file = paste0(res.tab,'/table_S1_8_II_L2.csv'), row.names = FALSE)


# Figure S1.4
correct.k.boxplots(setting, folder = 'I_L15', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.i.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_4.pdf")
# Figure S1.5
correct.k.boxplots(setting, folder = 'II_L1', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.ii.1,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_5.pdf")
# Figure S1.6
correct.k.boxplots(setting, folder = 'II_L2', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = unif.ii.2,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true,
                   fig.name = "figure_S1_6.pdf")





