################################################################################
# 5.2 Simulation results (under setting A)
################################################################################
# This file reproduces Tables and Figures in Section 5.2 of the Manuscript.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

library(dplyr)
library(gt)

source('sim_functions/perf_measures_ng.R')
source('sim_functions/perf_measures_jmdf.R')

gt.tables <- function(df, title){
  df %>% gt() %>% fmt_number(columns = -1, decimals = 3) %>%
    tab_header(title)
}


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



# Data setting A 
setting = 'A'

# Joinf frailty model by NG et al
#----------------------------------------------------------------------------
tables.ng = perf.measures.ng.all.scenarios(setting, scenarios = 1:9,
                                           beta = beta.true, gamma = gamma.true, 
                                           coef.ng.true = true.ng.param, N = 500)

# JMDF with Gaussian initialization (ii), L = 1.5
#-------------------------------------------------------------------------------
gauss.ii.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L15', init = 'gauss', 
                                               scenarios = 1:9, 
                                               beta = beta.true, gamma = gamma.true,
                                               true.K = true.K.masses)



# JMDF with Gaussian initialization (ii), L = 1.5
#-------------------------------------------------------------------------------
unif.ii.15 = perf.measures.jmdf.all.scenarios(setting, folder = 'II_L15', init = 'unif', 
                                              scenarios = 1:9, 
                                              beta = beta.true, gamma = gamma.true,
                                              true.K = true.K.masses)


# 5.2.1 Simulation results: first three scenarios
#---------------------------------------------------------------------------
# Table 2
gt.tables(tables.ng[tables.ng$scenario %in% 1:3,], title = "Table 2: Ng et al.")
gt.tables(gauss.ii.15$perf.coefs[gauss.ii.15$perf.coefs$scenario %in% 1:3,], 
          title = "Table 2: Gaussian (ii), L = 1.5")
gt.tables(unif.ii.15$perf.coefs[unif.ii.15$perf.coefs$scenario %in% 1:3,], 
          title = "Table 2: Uniform (ii), L = 1.5")

# Table 3
gt.tables(gauss.ii.15$K.est, title = "Table 3: Gaussian (ii), L = 1.5")
gt.tables(unif.ii.15$K.est, title = "Table 3: Uniform (ii), L = 1.5")


# 5.2.2 Simulation results: last six scenarios
#---------------------------------------------------------------------------
# Table 4 
gt.tables(tables.ng[tables.ng$scenario %in% 4:9,], title = "Table 4: Ng et al.")

# Table 5 --> Simulation results: first three scenarios (Section 5.2.2)
gt.tables(gauss.ii.15$perf.coefs[gauss.ii.15$perf.coefs$scenario %in% 4:9,], 
          title = "Table 5: Gaussian (ii), L = 1.5")
gt.tables(unif.ii.15$perf.coefs[unif.ii.15$perf.coefs$scenario %in% 4:9,], 
          title = "Table 5: Uniform (ii), L = 1.5")

# Figure 8
correct.k.boxplots(setting, folder = 'II_L15', init = 'gauss', 
                   scenarios.K = 4:9, true.K = true.K.masses,
                   performances = gauss.ii.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)
# Figure 9
correct.k.boxplots(setting, folder = 'II_L15', init = 'unif', 
                   scenarios.K = 4:9, true.K = true.K.masses, 
                   performances = unif.ii.15,
                   Pu = Pu.true, Pv = Pv.true, weights = weights.true)


