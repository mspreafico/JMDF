################################################################################
# Audit workflow - First 5 replications, setting A, scenarios 1-9
################################################################################
# This script generates a reduced set of intermediate simulation results for
# reproducibility checks. For each method, the first five replications of each
# of the nine scenarios are fitted and saved in the corresponding audit folder.
#
# The resulting files are compared with the corresponding results from the
# full simulation workflow using 06_audit_reproducibility_check.R.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

# Number of cores
available = parallel::detectCores()
n.cores = if (available >= 5) 5 else max(1, available - 1)

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
run.simulations.jmdf.gauss(setting = 'A', folder = 'II_L15', 
                           audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 3.283 mins


# JMDF Uniform initialization (ii) with L = 1.5
#----------------------------------------------------------------------
run.simulations.jmdf.unif(setting = 'A', folder = 'II_L15', 
                          audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 3.202 mins


# Ng et al. Model
#----------------------------------------------------------------------
run.simulations.ng(setting = 'A', audit = TRUE, audit.reps = 5)
# Total processing time on 5 cores: Time difference of 2.484 hours


