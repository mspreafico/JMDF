library(data.table)

perf.measures.jmdf <- function(df.coef, df.se, true.values){
  n.coef = ncol(df.coef)
  n.sim = dim(df.coef)[1]
  performance = NULL
  for(i in 1:n.coef){
    avg.coef = mean(df.coef[,i])
    bias = sum(df.coef[,i]-true.values[i])/n.sim
    low.ci.coef = df.coef[,i] - qnorm(0.975)*df.se[,i]
    up.ci.coef = df.coef[,i] + qnorm(0.975)*df.se[,i]
    inside.ci = ifelse(low.ci.coef <= true.values[i] & up.ci.coef >= true.values[i], 1, 0)
    inside.ci.corrected = ifelse(low.ci.coef <= avg.coef & up.ci.coef >= avg.coef, 1, 0)
    performance = rbind.data.frame(performance,
                                   cbind.data.frame('coef' = colnames(df.coef)[i],
                                                    'true' = true.values[i],
                                                    'avg' = avg.coef,
                                                    'bias' = bias,
                                                    'empSE' = sqrt(sum((df.coef[,i]-avg.coef)^2)/(n.sim-1)),
                                                    'MSE' = sum((df.coef[,i]-true.values[i])^2)/n.sim,
                                                    'coverage' = sum(inside.ci)/n.sim,
                                                    'bias.elim.coverage' = sum(inside.ci.corrected)/n.sim
                                   )
    )
  }
  return(performance)
}

perf.measures.jmdf.all.scenarios <- function(setting, folder, init, scenarios, beta, gamma, true.K){ 
  
  res.path = paste0('sim_results_',setting,'/',folder)
  
  perf.coefs = NULL
  K.est = NULL
  ord.frail.res = NULL
  for(s in scenarios){
    df.path = paste0(res.path,'/s',s,'_JMDF_',init,folder,'.Rdata')
    load(df.path)
    perf.s = perf.measures.jmdf(df.coef = coef.results[,2:5],
                                df.se = coef.results[,6:9], true.values = c(beta,gamma))
    perf.coefs = rbind.data.frame(perf.coefs,
                                  cbind.data.frame('scenario' = s,
                                                   perf.s)
    )
    tab = c('collapse'= dim(coef.results[coef.results$collapse==TRUE,])[1],
            table(factor(coef.results[coef.results$collapse==FALSE,]$K, levels=c(2:7))))
    tab.perc = round(tab/dim(coef.results)[1]*100,1)
    tab.N.perc = paste0(tab.perc,'% (',tab,')')
    names(tab.N.perc) = names(tab)
    K.est = rbind.data.frame(K.est,
                             cbind.data.frame('scenario' = s,
                                              'true K' = true.K[s],
                                              t(tab.N.perc))
    )
    ord.frail = arrange(frail.results, rep_b, Pv, Pu)
    ord.frail = data.table(ord.frail)
    ord.frail[, ord.group := seq(.N), by = rep_b]
    ord.frail.res[[s]] = ord.frail
  }
  
  return(list('perf.coefs'=perf.coefs, 'K.est'=K.est, 'ord.frail.res' = ord.frail.res))
}


correct.k.boxplots <- function(setting, folder, init, scenarios.K, performances,
                               true.K, Pu, Pv, weights,
                               colors = c('dodgerblue3','mediumseagreen','orangered','maroon3')){
  
  K.est = performances$K.est
  ord.frail.res = performances$ord.frail.res
  
  res.fig = paste0('sim_results_',setting,'/figures')
  if (!dir.exists(res.fig)) {
    dir.create(res.fig)
  }
  init.letter = ifelse(init=='unif','U','G')
  fig.name = paste0(res.fig,'/fig',setting,'_',init.letter,folder,'.pdf')
  pdf(fig.name, width = 16, height = 9)
  layout(matrix(c(1:(3*length(scenarios.K))), 3, length(scenarios.K), byrow = FALSE), respect = TRUE)
  par(mar = c(5, 5, 1, 1), oma = c(0, 0, 4, 0))  
  
  titles = NULL
  for(s in scenarios.K){
    titles[[s-3]] = c(paste0("Scenario ",s),paste0("= ",K.est[s,2+as.numeric(true.K[s])]),true.K[s])
    if(length(ord.frail.res[[s]][K==as.numeric(true.K[s])]$Pu)>0){
      # R frailty
      boxplot(ord.frail.res[[s]][K==as.numeric(true.K[s])]$Pu ~ ord.frail.res[[s]][K==as.numeric(true.K[s])]$ord.group,
              xlab = 'Mass (k)', ylab = 'Frailty Recurrent process (u)', col=colors[1:true.K[s]],
              cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75)
      for(k in 1:true.K[s]){
        abline(h = Pu[[s]][k], lty = 2, lwd = 2, col=colors[k])
      } 
      # D frailty
      boxplot(ord.frail.res[[s]][K==as.numeric(true.K[s])]$Pv ~ ord.frail.res[[s]][K==as.numeric(true.K[s])]$ord.group,
              xlab = 'Mass (k)', ylab = 'Frailty Death process (v)', col=colors[1:true.K[s]], cex.lab = 1.5, cex.axis = 1.25)
      for(k in 1:true.K[s]){
        abline(h = Pv[[s]][k], lty = 2, lwd = 2, col=colors[k])
      }
      boxplot(ord.frail.res[[s]][K==as.numeric(true.K[s])]$w ~ ord.frail.res[[s]][K==as.numeric(true.K[s])]$ord.group,
              xlab = 'Mass (k)', ylab = 'Weights (w)', col=colors[1:true.K[s]], ylim = c(0,1), 
              cex.lab = 1.5, cex.axis = 1.25)
      for(k in 1:true.K[s]){
        abline(h = weights[[s]][k], lty = 2, lwd = 2, col=colors[k])
      }
      
    }else{
      # R frailty
      boxplot(rep(0,as.numeric(true.K[s])) ~ factor(1:as.numeric(true.K[s])), border='white',
              ylim = c(min(Pu[[s]])-0.5,max(Pu[[s]]+0.5)),
              xlab = 'Mass (k)', ylab = 'Frailty Recurrent process (u)',
              cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75)
      for(k in 1:true.K[s]){
        abline(h = Pu[[s]][k], lty = 2, lwd = 2, col=colors[k])
      }
      # D frailty
      boxplot(rep(0,as.numeric(true.K[s])) ~ factor(1:as.numeric(true.K[s])), border='white',
              ylim = c(min(Pv[[s]])-0.5,max(Pv[[s]]+0.5)),
              xlab = 'Mass (k)', ylab = 'Frailty Death process (v)', col=colors[1:true.K[s]], cex.lab = 1.5, cex.axis = 1.25)
      for(k in 1:true.K[s]){
        abline(h = Pv[[s]][k], lty = 2, lwd = 2, col=colors[k])
      }
      boxplot(rep(0,as.numeric(true.K[s])) ~ factor(1:as.numeric(true.K[s])), border='white',
              xlab = 'Mass (k)', ylab = 'Weights (w)', col=colors[1:true.K[s]], ylim = c(0,1), 
              cex.lab = 1.5, cex.axis = 1.25)
      for(k in 1:true.K[s]){
        abline(h = weights[[s]][k], lty = 2, lwd = 2, col=colors[k])
      }
      
    }
  } 
  
  for (i in 1:6) {
    mtext(bquote(.(titles[[i]][1])), outer = TRUE, side = 3,
          line = 0.75, at = (i - 0.4) / 6, cex = 1.6, font = 2)
    mtext(bquote(N[K==.(titles[[i]][3])] ~ .(titles[[i]][2])), 
          outer = TRUE, side = 3, line = -1.5, at = (i - 0.4) / 6, cex = 1.2)
  }
  
  print(paste0('Figure "fig',setting,'_',init.letter,folder,'.pdf" saved in folder ','sim_results_',setting,'/figures'))
  
  dev.off()
}
