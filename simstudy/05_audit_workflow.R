################################################################################
# Audit workflow - First 5 replications, setting A, scenarios 1-9
################################################################################
# This script generates a reduced set of intermediate simulation results for
# reproducibility checks. For each method, the first five replications of each
# of the nine scenarios are fitted and saved in the corresponding audit folder.
#
# The resulting files are compared with the corresponding results obtained
# by the authors using 06_audit_reproducibility_check.R.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

# Number of cores
if (.Platform$OS.type == "windows") {
  # mclapply/pmclapply supports only mc.cores = 1 on Windows
  n_cores = 1
} else {
  available = parallel::detectCores()
  n_cores = if (available >= 5) 5 else max(1, available - 1)
}

# Load functions
source('sim_functions/sim_data_joint_frail.R')
source('../functions/JointFrailtyNg.R')
source('../functions/JMdiscfrail.R')
source('sim_functions/init_jmdf.R')
source('sim_functions/run_simulations.R')
source('sim_functions/extract_sim_results.R')

library(mcprogress)

################################################################################
# Reduced audit simulations
#
# Setting A and initialization (ii) with L = 1.5 are used as a representative
# configuration for the reproducibility audit. The first five replications
# are retained for each of the nine scenarios.
################################################################################

# JMDF Gaussian initialization (ii) with L = 1.5
#----------------------------------------------------------------------
run.simulations.jmdf.gauss(setting = 'A', folder = 'II_L15', n.cores = n_cores,
                           audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 1.25 mins
# Total processing time on 1 core: Time difference of 2.37 mins


# JMDF Uniform initialization (ii) with L = 1.5
#----------------------------------------------------------------------
run.simulations.jmdf.unif(setting = 'A', folder = 'II_L15', n.cores = n_cores,
                          audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 2.01 mins
# Total processing time on 1 core: Time difference of 3.22 mins


# Ng et al. Model
#----------------------------------------------------------------------
run.simulations.ng(setting = 'A', n.cores = n_cores,
                   audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 2.572956 hours


