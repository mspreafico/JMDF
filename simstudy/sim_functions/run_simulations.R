################################################################################
# Simulation workflows
#
# The functions below are used both for the full simulation study and for the
# reduced reproducibility audit. When audit = TRUE, only the first
# 'audit.reps' replications are fitted and the corresponding intermediate
# results (including the model fits) are saved in an 'audit' subdirectory.
################################################################################

nreps_audit = 5

# Function to run Ng et al.'s simulations under a given setting
#-----------------------------------------------------------------
run.simulations.ng <- function(setting, audit = FALSE, audit.reps = nreps_audit){
  
  data.dir = paste0('sim_data_',setting)
  # Create folder
  if(audit){
    res.dir = paste0('sim_results_ng_',setting,'/audit')
  }else{
    res.dir = paste0('sim_results_ng_',setting)
  }
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
    # For the audit workflow, retain only the first audit.reps replications.
    if(audit){dat.sim.df <- dat.sim.df[1:audit.reps]}
    # Fit Ng at al's JM
    ng.results <- pmclapply(dat.sim.df , function(x) joint.frailty.Ng(x, patient=x[,2], 
                                                                      theta01=0.1, theta02=0.1, rho0=0.1, 
                                                                      itmax=300), 
                            mc.cores = n.cores, mc.preschedule = FALSE) 
    sim.results = ng.jm.results(ng.results)
    # Save
    res.path = paste0(res.dir,'/s',s,'_JMNg.Rdata')
    # Store the full model results for audit checks; for the full workflow,
    # only the extracted simulation results are stored.
    if(audit){
      save(ng.results, sim.results, file = res.path)
    }else{
      save(sim.results, file = res.path)
    }
    cat(paste0('Results saved in file: ',res.path))
    cat('\n')
    # Processing time
    end_s = Sys.time()
    cat('Processing ended after: ')
    print(end_s - start_s)
    cat('\n')
  }
  end_time = Sys.time()
  cat('Total processing time: ')
  print(end_time-start_time)
}


# Function to run JMDF simulations with Gaussian initializations
# under a given setting and initialization/threshold
#-------------------------------------------------------------------------------
# parameter 'setting' = A or B:
#   - Setting A: yielding approximately 80% death and 20% administrative censoring
#   - Setting B: yielding approximately 50% death and 50% administrative censoring
#-------------------------------------------------------------------------------
# parameter 'folder' = 'I_L15', 'II_L1', 'II_L15', or 'II_L2', 
# which is a combination of
#   - implemented initialization:
#       - Init (i): mu = c(0,0) & Sigma = matrix(c(2,0.2,0.2,2), nrow = 2, ncol = 2)
#       - Init (ii): mu = c(0,0) & Sigma = matrix(c(0.2,0,0,2), nrow = 2, ncol = 2)
#   - and threshold:
#       - L = 1, 1.5, 2
#-------------------------------------------------------------------------------
run.simulations.jmdf.gauss <- function(setting, folder, audit = FALSE, audit.reps = nreps_audit){
  
  # Create directory for saving results
  data.dir = paste0('sim_data_',setting)
  res.dir = paste0('sim_results_',setting)
  if (!dir.exists(res.dir)) {
    dir.create(res.dir)
  }
  
  # Create folder
  if(audit){
    res.dir = paste0(res.dir,'/audit')
  }else{
    res.dir = paste0(res.dir,'/',folder)
  }
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
    # For the audit workflow, retain only the first audit.reps replications.
    if(audit){dat.sim.df <- dat.sim.df[1:audit.reps]}
    # Fit JMDF with Gaussian initialization
    jmdf.results <- pmclapply(dat.sim.df , function(x) JMdiscfrail(x$dataR, formulaR = '~ X1 + X2',
                                                                  x$dataD, formulaD = '~ X1 + X2',
                                                                  init.unif = FALSE,
                                                                  distance = "euclidean",
                                                                  Sigma = Sigma.0,
                                                                  mu = mu.0,
                                                                  M = 500, L = L.mindist,
                                                                  max.it = 300,  toll = 1e-3), 
                             mc.cores = n.cores, mc.preschedule = FALSE) 
    coef.results = jmdf.coef.res(jmdf.results)
    frail.results = jmdf.frail.res(jmdf.results)
    # Save
    res.path = paste0(res.dir,'/s',s,'_JMDF_gauss',folder,'.Rdata')
    # Store the full model results for audit checks; for the full workflow,
    # only the extracted simulation results are stored.
    if(audit){
      save(jmdf.results, coef.results, frail.results, file = res.path)
    }else{
      save(coef.results, frail.results, file = res.path)
    }
    cat(paste0('Results saved in file: ',res.path))
    cat('\n')
    # Processing time
    end_s = Sys.time()
    cat('Processing ended after: ')
    print(end_s - start_s)
    cat('\n')
  }
  end_time = Sys.time()
  cat('Total processing time: ')
  print(end_time-start_time)
  
}



# Function to run JMDF simulations with Uniform initializations
# under a given setting and initialization/threshold
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
run.simulations.jmdf.unif <- function(setting, folder, audit = FALSE, audit.reps = nreps_audit){
  
  # Create directory for saving results
  data.dir = paste0('sim_data_',setting)
  res.dir = paste0('sim_results_',setting)
  if (!dir.exists(res.dir)) {
    dir.create(res.dir)
  }
  
  # Create folder
  if(audit){
    res.dir = paste0(res.dir,'/audit')
  }else{
    res.dir = paste0(res.dir,'/',folder)
  }
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
    # For the audit workflow, retain only the first audit.reps replications.
    if(audit){dat.sim.df <- dat.sim.df[1:audit.reps]}
    # Fit JMDF with Uniform initalization
    jmdf.results <- pmclapply(dat.sim.df , function(x) JMdiscfrail(x$dataR, formulaR = '~ X1 + X2',
                                                                  x$dataD, formulaD = '~ X1 + X2',
                                                                  init.unif = TRUE,
                                                                  distance = "euclidean",
                                                                  ulim.unif = ulim.0,
                                                                  vlim.unif = vlim.0,
                                                                  M = 500, L = L.mindist,
                                                                  max.it = 300,  toll = 1e-3),
                             mc.cores = n.cores, mc.preschedule = FALSE) 
    coef.results = jmdf.coef.res(jmdf.results)
    frail.results = jmdf.frail.res(jmdf.results)
    # Save
    res.path = paste0(res.dir,'/s',s,'_JMDF_unif',folder,'.Rdata')
    # Store the full model results for audit checks; for the full workflow,
    # only the extracted simulation results are stored.
    if(audit){
      save(jmdf.results, coef.results, frail.results, file = res.path)
    }else{
      save(coef.results, frail.results, file = res.path)
    }
    cat(paste0('Results saved in file: ',res.path))
    cat('\n')
    # Processing time
    end_s = Sys.time()
    cat('Processing ended after: ')
    print(end_s - start_s)
    cat('\n')
  }
  end_time = Sys.time()
  cat('Total processing time: ')
  print(end_time-start_time)
  
}
