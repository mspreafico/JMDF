library(ggplot2)
library(ggpubr)
library(data.table)

source('functions/stratified_base_surv.R')

create.data.masses <- function(gauss3, gauss2, unif3, unif2){
  
  data_masses = data.table(rbind.data.frame(cbind.data.frame('init' = rep('Gauss',3),
                                                             'K' = rep(3,3),
                                                             'u' = gauss3$P[,1],
                                                             'v' = gauss3$P[,2],
                                                             'w' = gauss3$w,
                                                             'SE_u' = gauss3$se.P[,1],
                                                             'SE_v' = gauss3$se.P[,2]),
                                            cbind.data.frame('init' = rep('Gauss',2),
                                                             'K' = rep(2,2),
                                                             'u' = gauss2$P[,1],
                                                             'v' = gauss2$P[,2],
                                                             'w' = gauss2$w,
                                                             'SE_u' = gauss2$se.P[,1],
                                                             'SE_v' = gauss2$se.P[,2]),
                                            cbind.data.frame('init' = rep('Unif',3),
                                                             'K' = rep(3,3),
                                                             'u' = unif3$P[,1],
                                                             'v' = unif3$P[,2],
                                                             'w' = unif3$w,
                                                             'SE_u' = unif3$se.P[,1],
                                                             'SE_v' = unif3$se.P[,2]),
                                            cbind.data.frame('init' = rep('Unif',2),
                                                             'K' = rep(2,2),
                                                             'u' = unif2$P[,1],
                                                             'v' = unif2$P[,2],
                                                             'w' = unif2$w,
                                                             'SE_u' = unif2$se.P[,1],
                                                             'SE_v' = unif2$se.P[,2])
  ))
  setorder(data_masses, init, -K, u) 
  
  return(data_masses)
}

plot.masses <- function(result, colors=NULL){
  
  if(is.null(colors)){colors = 1:(result$K)}
  
  masses = paste0('(',sprintf("%.3f", round(result$P[,1],3)),', ',
                  sprintf("%.3f", round(result$P[,2],3)),')')
  
  dev.new()
  plot(result$P[,1], result$P[,2], cex=result$w*10,
       xlim = range(result$P[,1])+c(-0.1,0.1), 
       ylim = range(result$P[,2])+c(-0.1,0.1),
       col=colors, pch=16,
       xlab='Recurrent (u)', ylab='Terminal (v)', 
       main='Distribution of masses',
       cex.axis=1.3, cex.lab = 1.3, cex.main=1.3)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  legend(0.1,-0.1, col=colors, 
         pch = rep(16,result$K), y.intersp=1,
         legend = paste0('P',seq(1,result$K),': ',masses,' - w',seq(1,result$K),': ',round(result$w,3)))
  
  recordPlot()
}

#----------#
# Figure 3 #
#----------#
figure.masses.combined <- function(gauss3, gauss2, unif3, unif2){
  
  data_masses = create.data.masses(gauss3, gauss2, unif3, unif2)
  data_masses$shape = ifelse(data_masses$K==2,23,21)
  data_masses$id = c(1:10)
  data_masses[, coordinates := paste0('(',sprintf("%.3f", round(u,3)),', ',
                                      sprintf("%.3f", round(v,3)),')')]
  blue1 = '#3399FF'
  blue2 = '#0066CC'
  or1 = '#CC3300'
  or2 = '#FF6600'
  or3 = '#FF9933'
  pi1 = '#CC0066'
  pi2 = '#FF33CC'
  pi3 = '#FF6699'
  colori = c(or3,or2,or1,blue2,blue1,pi3,pi2,pi1,blue2,blue1)
  
  labels_names = paste0("P",c(1:3,1:2,1:3,1:2),": ",round(data_masses$w,3))
  
  p1=ggplot(data_masses[init=='Gauss'], aes(x=u, y=v, 
                                            fill=as.factor(id), color=as.factor(id), 
                                            shape=as.factor(id), label = coordinates)) +
    geom_vline(xintercept = 0.0, linetype='dotted', color='gray20') +
    geom_hline(yintercept = 0.0, linetype='dotted', color='gray20') +
    geom_point(aes(size=w), stroke = 2) +
    geom_text(hjust=c(-0.1,-0.1,-0.1,-0.15,1.1), size=3.5, show.legend = FALSE) +
    scale_fill_manual(name = "Masses & Weights", values = colori[1:5],
                      labels = labels_names[1:5]) +
    scale_color_manual(name = "Masses & Weights", values = colori[1:5],
                       labels = labels_names[1:5]) +
    scale_shape_manual(name = "Masses & Weights", values = c(21,21,21,23,23),
                       labels = labels_names[1:5]) +
    scale_x_continuous(limits=c(-0.5,1.5)) +
    scale_y_continuous(limits=c(-3,4)) +
    scale_size_continuous(guide = "none") +
    xlab("Recurrent (u)") +
    ylab("Terminal (v)") +
    ggtitle('Mass distribution - Gaussian case') +
    theme_light() +
    theme(legend.position='right', 
          axis.text=element_text(size=rel(1.1)),
          axis.title=element_text(size=rel(1.2)),
          plot.title = element_text(face="bold", size=rel(1.3)), 
          legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))
  
  p2=ggplot(data_masses[init=='Unif'], aes(x=u, y=v, 
                                           fill=as.factor(id), color=as.factor(id), 
                                           shape=as.factor(id), label=coordinates)) + 
    geom_vline(xintercept = 0.0, linetype='dotted', color='gray20') +
    geom_hline(yintercept = 0.0, linetype='dotted', color='gray20') +
    geom_point(aes(size=w), stroke = 2) +
    geom_text(hjust=c(-0.1,1.1,1.1,-0.15,1.1), size=3.5, show.legend = FALSE) +
    scale_fill_manual(name = "Masses & Weights", values = colori[6:10],
                      labels = labels_names[6:10]) +
    scale_color_manual(name = "Masses & Weights", values = colori[6:10],
                       labels = labels_names[6:10]) +
    scale_shape_manual(name = "Masses & Weights", values = c(21,21,21,23,23),
                       labels = labels_names[6:10]) +
    scale_x_continuous(limits=c(-0.5,1.5)) +
    scale_y_continuous(limits=c(-3,4)) +
    scale_size_continuous(guide = "none") +
    xlab("Recurrent (u)") +
    ylab("Terminal (v)") +
    ggtitle('Mass distribution - Uniform case') +
    theme_light() +
    theme(legend.position='right', 
          axis.text=element_text(size=rel(1.1)),
          axis.title=element_text(size=rel(1.2)),
          plot.title = element_text(face="bold", size=rel(1.3)), 
          legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=rel(1.2)))
  
  return(ggarrange(p1,p2))
} 


#----------#
# Figure 4 #
#----------#
figure.survival.curves <- function(result, colors=NULL){
  
  if(is.null(colors)){colors = 1:(result$K)}
  
  survival = stratified.base.surv(result)
  
  par(mfrow=c(1,2))
  survR = as.data.frame(survival$recurrent)
  plot(survR[,1], survR[,2], type='l', lty=2,  ylim =c(0,1),
       xlab = 'Time', ylab='Probability of no recurrence',
       main = 'Stratified baseline survival curve for recurrent events')
  for(k in 1:(result$K)){
    points(survR[,1], survR[,(k+2)], type='l', col=colors[k])
  }
  legend('topright', col=c('black', colors), 
         lty = c(2,rep(1,result$K)),
         legend = colnames(survR[,2:(result$K+2)]))
  
  survD = as.data.frame(survival$terminal)
  plot(survD[,1], survD[,2], type='l', lty=2, ylim =c(0,1),
       xlab = 'Time', ylab='Survival Probability',
       main = 'Stratified baseline survival curve for terminal event')
  for(k in 1:(result$K)){
    points(survD[,1], survD[,(k+2)], type='l', col=colors[k])
  }
  legend('bottomright', col=c('black',colors), 
         lty = c(2,rep(1,result$K)),
         legend = colnames(survD[,2:(result$K+2)]))
  
  recordPlot()
}


#----------#
# Figure 6 #
#----------#
frailty.comparisons <- function(result.rondeau, result.ng, result, colors=NULL){
  
  if(is.null(colors)){colors = 1:(result$K)}
  
  library(latex2exp)
  
  par(mfrow=c(1,2))
  boxplot(exp(result.rondeau$frailty.pred) ~ result$id.subgroups$subgroup,
          col=colors[result$id.subgroups$subgroup],
          xlab = "Mass-point", ylab = expression("Frailty term exp"(eta[i])),
          main = 'Frailty estimates by Rondeau and others')
  legend("topright", legend = paste0("P", 1:length(colors)),
         col = colors, pch = 19, bty = "n")
  
  plot(result.ng$u, result.ng$v, col=colors[result$id.subgroups$subgroup],
       xlab = "Recurrent (u)", ylab = "Terminal (v)",
       main='Frailty estimates by Ng and others', pch=19)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  legend("topleft", legend = paste0("P", 1:length(colors)),
         col = colors, pch = 19, bty = "n")
  
  recordPlot()
}

