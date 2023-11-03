# Function that computes the frailty-stratified BASELINE survival curves 
# for recurrent and terminal processes from a fitted JMdiscfrail.
# INPUT: jmdf.out = output of function JMdiscfrail()

stratified.base.surv = function(jmdf.out){
  
  Pu <- jmdf.out$P[,1]
  Pv <- jmdf.out$P[,2]
  K <- jmdf.out$K
  
  # Recurrent events - Uniform
  survR = cbind.data.frame('time' = jmdf.out$cumhazR$time,
                           'baseline' = exp(-jmdf.out$cumhazR$cumhaz)
  )
  survD = cbind.data.frame('time' = jmdf.out$cumhazD$time,
                           'baseline' = exp(-jmdf.out$cumhazD$cumhaz)
  )
  for(k in 1:K){
    survR = cbind.data.frame(survR,
                             'Pu' = exp(-jmdf.out$cumhazR$cumhaz*exp(Pu[k]))
    )
    survD = cbind.data.frame(survD,
                             'Pv' = exp(-jmdf.out$cumhazD$cumhaz*exp(Pv[k]))
    )
    
  } 
  row0 = c(0,1,rep(1,K))
  survR = data.table(rbind(row0,survR))
  colnames(survR)[3:(K+2)] = paste0('P',seq(1,K))
  survD = data.table(rbind(row0,survD))
  colnames(survD)[3:(K+2)] = paste0('P',seq(1,K))
  
  return(list('recurrent' = survR, 'terminal' = survD))
}