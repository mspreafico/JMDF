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
setting = 'A' #'B'
data.dir = paste0('sim_data_',setting)
res.dir = paste0('sim_results_',setting)
if (!dir.exists(res.dir)) {
  dir.create(res.dir)
}


# Uniform simulation results
#-----------------------------------------------------------------
# Implemented initializations:
#   - Init I: ulim.unif = c(-4,4) & vlim.unif = c(-4,4)
#   - Init II: ulim.unif = c(-1.5,1.5) & vlim.unif = c(-4,4)
#   - Init III: ulim.unif = c(-2,2) & vlim.unif = c(-1.5,1.5)
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
  # Processing time
  end_s = Sys.time()
  cat('Processing ended after:')
  print(end_s - start_s)
  cat('\n')
}
end_time = Sys.time()
cat('Total processing time:')
print(end_time-start_time)




