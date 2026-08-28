################################################################################
# Simulation results: Fitting joint frailty model by Ng et al.'s
################################################################################
# This file estimates joint frailty models according to Ng et al. for data/scenarios 
# generated under settinsa A or B.
#
# When possible, this code runs in parallel using either the number of available 
# cores minus 2, or up to a maximum of 40 cores.
################################################################################


# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))

# Number of cores
available = parallel::detectCores()
n.cores = if (available >= 40) 40 else max(1, available - 2)

library(mcprogress)

# Load functions
source('sim_functions/sim_data_joint_frail.R')
source('../functions/JointFrailtyNg.R')
source('sim_functions/run_simulations.R')
source('sim_functions/extract_sim_results.R')


# Setting A: yielding approximately 80% death and 20% administrative censoring
#-----------------------------------------------------------------------------
run.simulations.ng(setting = 'A')

# Setting B: yielding approximately 50% death and 50% administrative censoring
#-----------------------------------------------------------------------------
run.simulations.ng(setting = 'B')

