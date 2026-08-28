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
source('sim_functions/run_simulations.R')
source('sim_functions/extract_sim_results.R')

library(mcprogress)



# Uniform simulation results
##########################################################################
# 'setting' = A or B:
#   - Setting A: yielding approximately 80% death and 20% administrative censoring
#   - Setting B: yielding approximately 50% death and 50% administrative censoring
#-------------------------------------------------------------------------------
# 'folder' = 'I_L15', 'II_L1', 'II_L15', or 'II_L2', 
# which is a combination of
#   - implemented initialization:
#       - Init (i): ulim.unif = c(-4,4) & vlim.unif = c(-4,4)
#       - Init (ii): ulim.unif = c(-1.5,1.5) & vlim.unif = c(-4,4)
#   - and threshold:
#       - L = 1, 1.5, 2
#-------------------------------------------------------------------------------


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

