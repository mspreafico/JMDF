################################################################################
# Simulation results: Fitting JMDF using Uniform initialization
################################################################################
# This file estimates JMDF using Uniform initialization for data/scenarios 
# generated under different settings (A, B), initializations (i, ii), and 
# thresholds L = 1, 1.5, 2.
#
# When possible, this code runs in parallel using either the number of available 
# cores minus 2, or up to a maximum of 20 cores.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

# Number of cores
available = parallel::detectCores()
n.cores = if (available >= 20) 20 else max(1, available - 2)

# Load functions
source('sim_functions/sim_data_joint_frail.R')
source('../functions/JMdiscfrail.R')
source('sim_functions/init_jmdf.R')
source('sim_functions/extract_sim_results.R')

library(mcprogress)

# Function to run simulations under a given setting and initialization/threshold
#-------------------------------------------------------------------------------
# parameter 'setting' = A or B:
#   - Setting A: yielding approximately 80% death and 20% administrative censoring
#   - Setting B: yielding approximately 50% death and 50% administrative censoring
#-------------------------------------------------------------------------------
# parameter 'folder' = 'I_L15', 'II_L1', 'II_L15', or 'II_L2', 
# which is a combination of
#   - implemented initialization:
#       - Init (i): ulim.unif = c(-4,4) & vlim.unif = c(-4,4)
#       - Init (ii): ulim.unif = c(-1.5,1.5) & vlim.unif = c(-4,4)
#   - and threshold:
#       - L = 1, 1.5, 2
#-------------------------------------------------------------------------------
run.simulations.jmdf.unif <- function(setting, folder){
  
  # Create directory for saving results
  data.dir = paste0('sim_data_',setting)
  res.dir = paste0('sim_results_',setting)
  if (!dir.exists(res.dir)) {
    dir.create(res.dir)
  }
  
  # Create folder
  res.dir = paste0(res.dir,'/',folder)
  if (!dir.exists(res.dir)) {
    dir.create(res.dir)
  }
  
  ulim.0 = init.min.unif(folder)[[1]]
  vlim.0 = init.min.unif(folder)[[2]]
  L.mindist = init.min.unif(folder)[[3]]
  cat(paste0('Case ',folder,' as follows:
  Initialization: Pu limits = [',ulim.0[1],',',ulim.0[2],'], Pv limits = [',vlim.0[1],',',vlim.0[2],']
  Threshold mindist: L  = ', L.mindist))
  
  set.seed(20241106)
  start_time = Sys.time()
  for(s in 1:9){
    start_s = Sys.time()
    print(paste0('Processing Scenario ',s,' - JMDF Uniform Initialization'))
    # Load data
    df.path = paste0(data.dir,'/scenario',s,'.Rdata')
    load(df.path)
    # Simulations
    dat.sim.df <- pmclapply(dat.sim, function(x) dataformat.JMDF(x), 
                            mc.cores = n.cores, mc.preschedule = FALSE) 
    sim.results <- pmclapply(dat.sim.df , function(x) JMdiscfrail(x$dataR, formulaR = '~ X1 + X2',
                                                                  x$dataD, formulaD = '~ X1 + X2',
                                                                  init.unif = TRUE,
                                                                  distance = "euclidean",
                                                                  ulim.unif = ulim.0,
                                                                  vlim.unif = vlim.0,
                                                                  M = 500, L = L.mindist,
                                                                  max.it = 300,  toll = 1e-3),
                             mc.cores = n.cores, mc.preschedule = FALSE) 
    coef.results = jmdf.coef.res(sim.results)
    frail.results = jmdf.frail.res(sim.results)
    # Save
    res.path = paste0(res.dir,'/s',s,'_JMDF_unif',folder,'.Rdata')
    save(coef.results, frail.results, file = res.path)
    cat(paste0('Results saved in file: ',res.path))
    cat('\n')
    # Processing time
    end_s = Sys.time()
    cat('Processing ended after:')
    print(end_s - start_s)
    cat('\n')
  }
  end_time = Sys.time()
  cat('Total processing time:')
  print(end_time-start_time)
  
}


# Uniform simulation results
##########################################################################

# Simulations for Section 5.2 of the Manuscript
#-------------------------------------------------------------------------
# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting A, with initialization (ii) and threshold L = 1.5
run.simulations.jmdf.unif(setting = 'A', folder = 'II_L15')


# Simulations for Supplementary Material S1.2
#-------------------------------------------------------------------------
# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting A, with initialization (i) and threshold L = 1.5
run.simulations.jmdf.unif(setting = 'A', folder = 'I_L15')

# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting A, with initialization (ii) and threshold L = 1
run.simulations.jmdf.unif(setting = 'A', folder = 'II_L1')

# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting A, with initialization (ii) and threshold L = 2
run.simulations.jmdf.unif(setting = 'A', folder = 'II_L2')


# Simulations for Supplementary Material S2.3
#-------------------------------------------------------------------------
# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting B, with initialization (i) and threshold L = 1.5
run.simulations.jmdf.unif(setting = 'B', folder = 'I_L15')

# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting B, with initialization (ii) and threshold L = 1
run.simulations.jmdf.unif(setting = 'B', folder = 'II_L1')

# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting B, with initialization (ii) and threshold L = 1.5
run.simulations.jmdf.unif(setting = 'B', folder = 'II_L15')

# Estimate JMDF using Uniform initialization for data/scenarios generated 
# under setting B, with initialization (ii) and threshold L = 2
run.simulations.jmdf.unif(setting = 'B', folder = 'II_L2')

