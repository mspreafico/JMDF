library(data.table)
library(MASS)

simRecDeath.events <- function(df, maxrecevents, beta, gamma,
                               lambda_R, lambda_D, tau_R, tau_D, a, b){
  
  # Recurrence Data
  options("digits" = 4) # for a better display of data.table
  dataR = df[, .(id,X1,X2,u)]
  # Simulate number of recurrent events per patient
  dataR = dataR[, nevents := sample(x = c(0:maxrecevents), size=1, replace=TRUE,
                                    prob = rev((0:maxrecevents+1)/sum(0:maxrecevents+1))), by=c('id')]
  dataR = dataR[rep(1:.N,nevents)][,countR:=1:.N,by=c('id')]
  # Inverse cumulative hazard function
  dataR[, unifR := runif(1, 0, 1), by=c('id','countR')]
  dataR[, AR := -log(unifR)*exp(-(beta[1]*X1+beta[2]*X2+u))]
  dataR[, gaptime := (lambda_R^(-1) * AR)^(1/tau_R)]
  dataR[, deltaR := 1]
  dataR[, deltaD := 0]
  dataR[order(id,countR),]
  
  # Death Data
  options("digits" = 4) # for a better display of data.table
  dataD = df[, .(id,X1,X2,v)]
  dataD[, unifD := runif(1, 0, 1), by=c('id')]
  dataD[, AD := -log(unifD)*exp(-(gamma[1]*X1+gamma[2]*X2+v))]
  dataD[, gaptime := (lambda_D^(-1) * AD)^(1/tau_D)]
  dataD[, deltaR := 0]
  dataD[, deltaD := 1]
  dataD[, countR := NA]
  # Administrative follow-up
  dataD[, AdmFup := runif(1, a, b), by=c('id')]
  #dataD[, AdmFup := 120]
  
  # Merge Recurrent & Death Data
  data = rbind(dataR[,.(id,countR,X1,X2,gaptime,deltaR,deltaD)], 
               dataD[,.(id,countR,X1,X2,gaptime,deltaR,deltaD)])
  withR = unique(data[deltaR==1]$id)
  data[!(id %in% withR), countR := 0,]
  data[id %in% withR, maxev := max(countR, na.rm=T), by=c('id') ]
  data[id %in% withR & is.na(countR), countR:= maxev+1, by=c('id') ]
  data = merge(data, dataD[,.(id,AdmFup)], by='id')
  # Remove events after administrative follow-up
  data[, cumtime := cumsum(gaptime), by=c('id')]
  data[cumtime>=AdmFup, stop := min(countR, na.rm=T), by=c('id') ]
  data = data[is.na(stop) | countR==stop]
  # Fix censoring time & event/death indicators when administrative follow-up is reached
  data[countR==stop, c('deltaR','deltaD','cumtime') := list(0,0,AdmFup), by=c('id')]
  data[, gaptime := diff(c(0,cumtime)), by=c('id')]
  
  data.return = data[, .(id,countR,gaptime,deltaR,deltaD,X1,X2)]
  options("digits" = 4) # for a better display of data.table
  
  return(data.return)
}


simRecDeath <- function(df, maxrecevents, beta, gamma,
                        lambda_R, lambda_D, tau_R, tau_D, a, b){
  
  # Recurrence Data
  options("digits" = 4) # for a better display of data.table
  dataR = df[, .(id,X1,X2,u)]
  dataR = dataR[rep(seq_len(nrow(dataR)), maxrecevents), ]
  dataR[, countR := seq(.N), by='id']
  # Inverse cumulative hazard function
  dataR[, unifR := runif(1, 0, 1), by=c('id','countR')]
  dataR[, AR := -log(unifR)*exp(-(beta[1]*X1+beta[2]*X2+u))]
  dataR[, tR := (lambda_R^(-1) * AR)^(1/tau_R)]
  dataR[order(id,countR),]
  
  # Death Data
  options("digits" = 4) # for a better display of data.table
  dataD = df[, .(id,X1,X2,v)]
  dataD[, unifD := runif(1, 0, 1), by=c('id')]
  dataD[, AD := -log(unifD)*exp(-(gamma[1]*X1+gamma[2]*X2+v))]
  dataD[, tD := (lambda_D^(-1) * AD)^(1/tau_D)]
  # Administrative follow-up
  dataD[, AdmFup := runif(1, a, b), by=c('id')]
  #dataD[, AdmFup := 120]
  
  
  # Merge Recurrent & Death Data
  data = merge(dataR[,.(id,countR,X1,X2,tR)], dataD[,.(id,tD,AdmFup)], by='id')
  # Gaptimes & Event/Death indicators
  data[, deltaR := 0]
  data[, deltaD := 0]
  data[, gaptime := min(tR,tD), by=c('id','countR')]
  data[gaptime==tR, deltaR := 1]
  data[gaptime==tD, deltaD := 1]
  # Remove events after death
  data[deltaD==1, mincountR := min(countR), by=c('id')]
  withD = unique(data[deltaD==1]$id)
  data[id %in% withD, mincountR := min(mincountR, na.rm=T), by=c('id')]
  data = data[ !(id %in% withD) | (id %in% withD & countR<=mincountR) ]
  # Remove events after administrative follow-up
  data[, cumtime := cumsum(gaptime), by=c('id')]
  data[cumtime>=AdmFup, stop := min(countR, na.rm=T), by=c('id') ]
  data = data[is.na(stop) | countR==stop]
  # Fix censoring time & event/death indicators when administrative follow-up is reached
  data[countR==stop, c('deltaR','deltaD','cumtime') := list(0,0,AdmFup), by=c('id')]
  data[, gaptime := diff(c(0,cumtime)), by=c('id')]
  
  data.return = data[, .(id,countR,gaptime,deltaR,deltaD,X1,X2)]
  options("digits" = 4) # for a better display of data.table
  
  return(data.return)
}

dataformat.Ng<-function(dat){
  dat.ng = dat[, .(countR,id,gaptime,deltaR,deltaD,X1,X2)]
  colnames(dat.ng) = c("hosp","id","time","delta1","delta2","X1","X2")
  dat.ng = as.matrix(dat.ng)
  return(dat.ng)
}

simDataMvNormFrail <- function(N, maxrecevents,
                               beta, gamma, theta_u, theta_v, rho,
                               lambda_R, lambda_D, tau_R, tau_D, a, b){
  
  # Simulate multivariate normal random effects
  Sigma = matrix(c(theta_u^2, rho*theta_u*theta_v, rho*theta_u*theta_v, theta_v^2), 2, 2)
  q = mvrnorm(N, c(0, 0), Sigma)
  
  df = data.table(data.frame(
    'id' = seq(1,N),
    'X1' = rbinom(N, 1, 0.5),
    'X2' = rnorm(N),
    'u' = q[,1],
    'v' = q[,2]
  ))
  
  dat = simRecDeath(df, maxrecevents, beta, gamma,
                    lambda_R, lambda_D, tau_R, tau_D, a, b)
  
  return(dat)
}


simDataDiscFrail <- function(N, maxrecevents,
                             beta, gamma, P_matrix, weights,
                             lambda_R, lambda_D, tau_R, tau_D, a, b){
  
  # Simulate groups of discrete random effects
  K = nrow(P_matrix)
  groups = sample(x = c(1:K), size=N, replace=TRUE, prob=weights)
  
  df = data.table(data.frame(
    'id' = seq(1,N),
    'X1' = rbinom(N, 1, 0.5),
    'X2' = rnorm(N),
    'u' = P_matrix[groups,1],
    'v' = P_matrix[groups,2]
  ))
  
  dat = simRecDeath(df, maxrecevents, beta, gamma,
                    lambda_R, lambda_D, tau_R, tau_D, a, b)
  
  return(dat)
}


dataformat.JMDF<-function(dat){
  dat.df = dat[, .(id,countR,deltaR,deltaD,gaptime,X1,X2)]
  colnames(dat.df)[2] = "event"
  
  dat.df[, lastevent := max(event), by=c('id')]
  
  dat.dfD = dat.df[event==lastevent]
  dat.dfD[, lastevent := NULL]
  dat.df[, lastevent := NULL]
  
  return(list('dataR' = dat.df, 'dataD' = dat.dfD))
}
