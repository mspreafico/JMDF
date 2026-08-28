################################################################################
# Reproducibility inspection
################################################################################
# Inspect the fitted model objects for a selected scenario and replication.
# This script provides an example of how individual intermediate simulation
# results can be loaded and inspected directly.
################################################################################

# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
setwd(file.path(getwd(),'simstudy'))


# Scenario 1, Replication 1
scenario = 4
replication = 1

################################################################################
# Direct comparison with the corresponding full simulation result
################################################################################

# JMDF with Gaussian initialization (ii), L = 1.5
#------------------------------------------------------------------------------
load(paste0("sim_results_A/audit/s", scenario,"_JMDF_gaussII_L15.Rdata"))
gauss.fit = jmdf.results[[replication]]

# Fixed effects - Recurrent (betas)
gauss.fit$modelR
# Fixed effects - Terminal (gammas)
gauss.fit$modelD
# Baseline Cumulative Hazard for Recurrent event
gauss.fit$cumhazR
# Baseline Cumulative Hazard for Terminal event
gauss.fit$cumhazD
# Random effects
gauss.fit$K
gauss.fit$w
gauss.fit$P
gauss.fit$se.P
# Subjects' subgroups
head(gauss.fit$id.subgroup)
table(gauss.fit$id.subgroup$subgroup)


# JMDF with Uniform initialization (ii), L = 1.5
#------------------------------------------------------------------------------
load(paste0("sim_results_A/audit/s", scenario,"_JMDF_unifII_L15.Rdata"))
unif.fit = jmdf.results[[replication]]

# Fixed effects - Recurrent (betas)
unif.fit$modelR
# Fixed effects - Terminal (gammas)
unif.fit$modelD
# Baseline Cumulative Hazard for Recurrent event
unif.fit$cumhazR
# Baseline Cumulative Hazard for Terminal event
unif.fit$cumhazD
# Random effects
unif.fit$K
unif.fit$w
unif.fit$P
unif.fit$se.P
# Subjects' subgroups
head(unif.fit$id.subgroup)
table(unif.fit$id.subgroup$subgroup)


# Ng et al. Model
#----------------------------------------------------------------------
load(paste0("sim_results_ng_A/audit/s", scenario,"_JMNg.Rdata"))
ng.fit = ng.results[[replication]]

ng.fit

plot(ng.fit$u, ng.fit$v, xlab='u', ylab='v', main='Random Effects Estimates', pch=19)
abline(h=0, lty=2)
abline(v=0, lty=2)

