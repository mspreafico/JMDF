################################################################################
# 5.2. Simulation results: Fitting joint frailty model by Ng et al.'s
################################################################################
# This file estimates joint frailty models according to Ng et al. for data/scenarios 
# generated under setting A (yielding approximately 80% death and 20% administrative censoring).
#
# To estimate Ng et al.'s models for data generated under setting B 
# (yielding approximately 50% death and 50% administrative censoring):
# 1. First, modify the value of `setting` in "01_data_generation.R" to 'B' and run it 
#    to generate data under setting B.
# 2. Then, modify lines 35-36 in this file to set `setting` to 'B' and run it.
#
# When possible, this code runs in parallel using either the number of available 
# cores minus 2, or up to a maximum of 40 cores.
################################################################################


# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
file.path(getwd(),'simstudy')

# Number of cores
available = parallel::detectCores()
n.cores = if (available >= 40) 40 else max(1, available - 2)

library(mcprogress)

# Load functions
source('sim_functions/sim_data_joint_frail.R')
source('../functions/JointFrailtyNg.R')
source('sim_functions/extract_sim_results.R')


# Data setting
#-----------------------------------------------------------------
# SELECT generated datasets:
#   - Setting A yielding approximately 80% death and 20% administrative censoring
#   - Setting B yielding approximately 50% death and 50% administrative censoring
setting = 'A' 
#setting = 'B'

data.dir = paste0('sim_data_',setting)
# Create folder
res.dir = paste0('sim_results_ng_',setting)
if (!dir.exists(res.dir)) {
  dir.create(res.dir)
}

# Simulation study for Ng et al.'s model
start_time = Sys.time()
for(s in 1:9){
  start_s = Sys.time()
  print(paste0('Processing Scenario ',s,' - Ng et al Model'))
  # Load data
  df.path = paste0(data.dir,'/scenario',s,'.Rdata')
  load(df.path)
  # Estimate Ng et al.'s models
  dat.sim.df <- pmclapply(dat.sim, function(x) dataformat.Ng(x), 
                          mc.cores = n.cores, mc.preschedule = FALSE) 
  ng.results <- pmclapply(dat.sim.df , function(x) joint.frailty.Ng(x, patient=x[,2], 
                                                                     theta01=0.1, theta02=0.1, rho0=0.1, 
                                                                     itmax=300), 
                           mc.cores = n.cores, mc.preschedule = FALSE) 
  sim.results = ng.jm.results(ng.results)
  # Save
  res.path = paste0(res.dir,'/s',s,'_JMNg.Rdata')
  save(sim.results, file = res.path)
  # Processing time
  end_s = Sys.time()
  cat('Processing ended after:')
  print(end_s - start_s)
  cat('\n')
}
end_time = Sys.time()
cat('Total processing time:')
print(end_time-start_time)




