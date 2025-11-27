# Session -> Set Working Directory -> To Source File Location
library(data.table)

perf.measures.ng <- function(df.coef, df.se, true.values,
                             interval.type, N){
  n.coef = ncol(df.coef)
  n.sim = dim(df.coef)[1]
  performance = NULL
  for(i in 1:n.coef){
    avg.coef = mean(df.coef[,i])
    if(interval.type[i]=='norm'){
      # Confidence interval normal regression param
      low.ci.coef = df.coef[,i] - qnorm(0.975)*df.se[,i]
      up.ci.coef = df.coef[,i] + qnorm(0.975)*df.se[,i]
    }else if(interval.type[i]=='chisq'){
      # Confidence interval variance
      low.ci.coef = sqrt(df.coef[,i]*(N-1)/qchisq(0.975, df = N-1))
      up.ci.coef = sqrt(df.coef[,i]*(N-1)/qchisq(0.025, df = N-1))
    }else if(interval.type[i]=='rho'){
      # Confidence interval correlation coef
      low.ci.coef = df.coef[,i] - qt(0.975, df = N-2)*df.se[,i]
      up.ci.coef = df.coef[,i] + qt(0.975, df = N-2)*df.se[,i]
    }
    inside.ci = ifelse(low.ci.coef <= true.values[i] & up.ci.coef >= true.values[i], 1, 0)
    inside.ci.corrected = ifelse(low.ci.coef <= avg.coef & up.ci.coef >= avg.coef, 1, 0)
    performance = rbind.data.frame(performance,
                                   cbind.data.frame('coef' = colnames(df.coef)[i],
                                                    'true' = true.values[i],
                                                    'avg' = avg.coef,
                                                    'bias' = sum(df.coef[,i]-true.values[i])/n.sim,
                                                    'empSE' = sqrt(sum((df.coef[,i]-avg.coef)^2)/(n.sim-1)),
                                                    'MSE' = sum((df.coef[,i]-true.values[i])^2)/n.sim,
                                                    'coverage' = sum(inside.ci)/n.sim,
                                                    'bias.elim.coverage' = sum(inside.ci.corrected)/n.sim
                                   )
    )
  }
  return(performance)
}


perf.measures.ng.all.scenarios <- function(setting, scenarios,
                                           interval.type = c('norm','norm','norm','norm','chisq','chisq','rho'),
                                           beta, gamma, coef.ng.true, N){
  
  
  perf.ng.df = NULL
  for(s in scenarios){
    df.path = paste0('sim_results_ng_',setting,'/s',s,'_JMNg.Rdata')
    load(df.path)
    perf.s = perf.measures.ng(df.coef = cbind.data.frame(sim.results[,2:5],
                                                         'theta_u'= sqrt(sim.results[,6]),
                                                         'theta_v' = sqrt(sim.results[,7]),
                                                         'rho'=sim.results[,8]),
                              df.se = cbind.data.frame(sim.results[,9:12],
                                                       NA, NA,
                                                       sim.results[,15]),
                              true.values = c(beta,gamma,coef.ng.true[s,]),
                              interval.type, N)
    perf.ng.df = rbind.data.frame(perf.ng.df, cbind.data.frame('scenario' = s,
                                                            perf.s)
    )
  }
  
  return(perf.ng.df)
}
