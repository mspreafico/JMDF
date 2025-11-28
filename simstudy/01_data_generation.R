################################################################################
# 5.1. Data generation
################################################################################
# This file generates data as described in Manuscript Section 5.1, under setting
# A (yielding approximately 80% death and 20% administrative censoring).
#
# To generate data under setting B (yielding approximately 50% death and 50% 
# administrative censoring), modify the value of `setting` to 'B' (lines 33-34)
# and run it.
################################################################################


# Session -> Set Working Directory -> To Source File Location
#setwd('/home/spreaficom/jmdf/simstudy')
file.path(getwd(),'simstudy')

# Simulate B data sets for each scenario
source('sim_functions/sim_data_joint_frail.R')

# Number of simulated data sets for each setting and scenario
B = 200 

# Parameters
N=500
maxev=10

beta=c(-0.6, 0.8)
gamma=c(-0.8, 0.5)

a_min_fup=60
b_max_fup=120

# SELECT setting:
#   - Setting A yielding approximately 80% death and 20% administrative censoring
#   - Setting B yielding approximately 50% death and 50% administrative censoring
setting = 'A'
# setting = 'B'

# Shape parameters
tau_R=1.3 
tau_D=1.5
# Scale parameters
if(setting=='B'){
  lambda_R=5e-3
  lambda_D=2e-3
  data.dir = paste0('sim_data_',setting)
}else{
  lambda_R=2e-2
  lambda_D=1.4e-2
  data.dir = paste0('sim_data_',setting)
}

# Create folder
if (!dir.exists(data.dir)) {
  dir.create(data.dir)
}


#------------#
# SCENARIO 1 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with bivariate normal frailty
  dat.sim[[bb]] = simDataMvNormFrail(N=N, maxrecevents=maxev,
                             beta=beta, gamma=gamma,
                             theta_u=0.8, theta_v=0.8, rho=0.8,
                             lambda_R=lambda_R, lambda_D=lambda_D,
                             tau_R=tau_R, tau_D=tau_D,
                             a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario1.Rdata'))

#------------#
# SCENARIO 2 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with bivariate normal frailty
  dat.sim[[bb]] = simDataMvNormFrail(N=N, maxrecevents=maxev,
                                     beta=beta, gamma=gamma,
                                     theta_u=0.8, theta_v=0.8, rho=0.2,
                                     lambda_R=lambda_R, lambda_D=lambda_D,
                                     tau_R=tau_R, tau_D=tau_D,
                                     a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario2.Rdata'))

#------------#
# SCENARIO 3 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with bivariate normal frailty
  dat.sim[[bb]] = simDataMvNormFrail(N=N, maxrecevents=maxev,
                                     beta=beta, gamma=gamma,
                                     theta_u=0.03, theta_v=0.02, rho=0.2,
                                     lambda_R=lambda_R, lambda_D=lambda_D,
                                     tau_R=tau_R, tau_D=tau_D,
                                     a=a_min_fup, b=b_max_fup)
}

save(dat.sim, file = paste0(data.dir,'/scenario3.Rdata'))


#------------#
# SCENARIO 4 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 4 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                     beta=beta, gamma=gamma,
                                     P_matrix = matrix(c(-1.2,1.2, -1,1), nrow=2, byrow=F),
                                     weights = c(0.5,0.5),
                                     lambda_R=lambda_R, lambda_D=lambda_D,
                                     tau_R=tau_R, tau_D=tau_D,
                                     a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario4.Rdata'))


#------------#
# SCENARIO 5 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 5 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-1.2,0.51, -1,0.43), nrow=2, byrow=F),
                                   weights = c(0.3,0.7),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario5.Rdata'))


#------------#
# SCENARIO 6 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 6 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-1.2,1.2,0.68, -1,1.2,0.29), nrow=3, byrow=F),
                                   weights = c(0.45,0.32,0.23),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario6.Rdata'))


#------------#
# SCENARIO 7 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 7 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-1.2,2.2,-0.71, -1,0.7,0.98), nrow=3, byrow=F),
                                   weights = c(0.45,0.32,0.23),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario7.Rdata'))


#------------#
# SCENARIO 8 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 8 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-2,2,2,-2, -2,2,-2,2), nrow=4, byrow=F),
                                   weights = c(0.25,0.25,0.25,0.25),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario8.Rdata'))


#------------#
# SCENARIO 9 #
#------------#
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  # Simulate data with discrete frailty (see scenario 9 in Figure 7 of the manuscript)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-1.5,0.75,2.25, 0,0,0), nrow=3, byrow=F),
                                   weights = c(0.4,0.5,0.1),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario9.Rdata'))


#-------------#
# SCENARIO 10 #
#-------------#
# Additional scenario
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  dat.sim[[bb]] = simDataMvNormFrail(N=N, maxrecevents=maxev,
                                     beta=beta, gamma=gamma,
                                     theta_u=1, theta_v=2, rho=0.8,
                                     lambda_R=lambda_R, lambda_D=lambda_D,
                                     tau_R=tau_R, tau_D=tau_D,
                                     a=a_min_fup, b=b_max_fup)
}

save(dat.sim, file = paste0(data.dir,'/scenario10.Rdata'))


#-------------#
# SCENARIO 11 #
#-------------#
# Additional scenario
set.seed(1234)
pb = txtProgressBar(min = 0, max = B, style = 3, width = 50, char = "=")
dat.sim = NULL
for(bb in 1:B){
  setTxtProgressBar(pb, bb)
  dat.sim[[bb]] = simDataDiscFrail(N=N, maxrecevents=maxev,
                                   beta=beta, gamma=gamma,
                                   P_matrix = matrix(c(-1,0,1, -2,0,2), nrow=3, byrow=F),
                                   weights = c(1/3,1/3,1/3),
                                   lambda_R=lambda_R, lambda_D=lambda_D,
                                   tau_R=tau_R, tau_D=tau_D,
                                   a=a_min_fup, b=b_max_fup)
}
save(dat.sim, file = paste0(data.dir,'/scenario11.Rdata'))

