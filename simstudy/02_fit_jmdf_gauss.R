################################################################################
# 5.2. Simulation results: Fitting JMDF using Gaussian initialization
################################################################################
# This file estimates JMDF using Gaussian initialization for data/scenarios 
# generated under setting A (yielding approximately 80% death and 20% administrative censoring),
# initialization (ii) and threshold L = 1.5.
#
# To select different initialization (i or ii) and threshold L, modify the value 
# of `folder` (lines 62-65: 'I_L1', 'II_L1', 'II_L15', or 'II_L2'). 
# You must run the JMDF files separately for each initialization to get results 
# corresponding to all choices.
#
# To estimate JMDFs for data generated under setting B 
# (yielding approximately 50% death and 50% administrative censoring):
# 1. First, modify the value of `setting` in "01_data_generation.R" to 'B' and run it 
#    to generate data under setting B.
# 2. Then, modify lines 43-44 in this file to set `setting` to 'B';
# 3. Select the desired initialization (i or ii) and threshold L (`folder`) and run it.
#
# When possible, this code runs in parallel using either the number of available 
# cores minus 2, or up to a maximum of 20 cores.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
file.path(getwd(),'simstudy')

# Number of cores
available = parallel::detectCores()
n.cores = if (available >= 20) 20 else max(1, available - 2)

# Load functions
source('sim_functions/sim_data_joint_frail.R')
source('../functions/JMdiscfrail.R')
source('sim_functions/init_jmdf.R')
source('sim_functions/extract_sim_results.R')

library(mcprogress)


# Data setting
#-----------------------------------------------------------------
# SELECT generated datasets:
#   - Setting A yielding approximately 80% death and 20% administrative censoring
#   - Setting B yielding approximately 50% death and 50% administrative censoring
setting = 'A' 
#setting = 'B'

data.dir = paste0('sim_data_',setting)
res.dir = paste0('sim_results_',setting)
if (!dir.exists(res.dir)) {
  dir.create(res.dir)
}

# Gaussian simulation results
#-------------------------------------------------------------------------
# Implemented initializations:
#   - Init I: mu = c(0,0) & Sigma = matrix(c(2,0.2,0.2,2), nrow = 2, ncol = 2)
#   - Init II: mu = c(0,0) & Sigma = matrix(c(0.2,0,0,2), nrow = 2, ncol = 2)
#   - Init III: mu = c(0,0) & Sigma = matrix(c(0.6,0,0,0.2), nrow = 2, ncol = 2)
# and thresholds:
#   - L = 1, 1.5, 2
#---------------------------
# SELECT Init and L
#folder = 'I_L15'
#folder = 'II_L1'
folder = 'II_L15'
#folder = 'II_L2'
#folder = 'III_L15'

# Create folder
res.dir = paste0(res.dir,'/',folder)
if (!dir.exists(res.dir)) {
  dir.create(res.dir)
}

mu.0 = init.min.gauss(folder)[[1]]
Sigma.0 = init.min.gauss(folder)[[2]]
L.mindist = init.min.gauss(folder)[[3]]
cat(paste0('Case ',folder,' as follows:
  Initialization: mu vector = (',mu.0[1],',',mu.0[2],'),
                  Sigma matrix = [',Sigma.0[1,1],',',Sigma.0[1,2],', ',Sigma.0[2,1],',',Sigma.0[2,2],']
  Threshold mindist: L  = ', L.mindist))


set.seed(20241106)
start_time = Sys.time()
for(s in 1:9){
  start_s = Sys.time()
  print(paste0('Processing Scenario ',s,' - JMDF Gaussian Initialization'))
  # Load data
  df.path = paste0(data.dir,'/scenario',s,'.Rdata')
  load(df.path)
  # Simulations
  dat.sim.df <- pmclapply(dat.sim, function(x) dataformat.JMDF(x), 
                          mc.cores = n.cores, mc.preschedule = FALSE) 
  sim.results <- pmclapply(dat.sim.df , function(x) JMdiscfrail(x$dataR, formulaR = '~ X1 + X2',
                                                                x$dataD, formulaD = '~ X1 + X2',
                                                                init.unif = FALSE,
                                                                distance = "euclidean",
                                                                Sigma = Sigma.0,
                                                                mu = mu.0,
                                                                M = 500, L = L.mindist,
                                                                max.it = 300,  toll = 1e-3), 
                           mc.cores = n.cores, mc.preschedule = FALSE) 
  coef.results = jmdf.coef.res(sim.results)
  frail.results = jmdf.frail.res(sim.results)
  # Save
  res.path = paste0(res.dir,'/s',s,'_JMDF_gauss',folder,'.Rdata')
  save(coef.results, frail.results, file = res.path)
  # Processing time
  end_s = Sys.time()
  cat('Processing ended after:')
  print(end_s - start_s)
  cat('\n')
}
end_time = Sys.time()
cat('Total processing time:')
print(end_time-start_time)




