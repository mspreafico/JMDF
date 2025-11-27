jmdf.coef.res <- function(res){
  
  coef.res = NULL
  for(b in 1:length(res)){
    res.b = res[[b]]
    if(!inherits(res.b, "try-error") ){
      coef.res = rbind.data.frame(coef.res,
                                  cbind.data.frame('rep_b' = b,
                                                   'beta1' = as.vector(coef(res.b$modelR)[1]),
                                                   'beta2' = as.vector(coef(res.b$modelR)[2]),
                                                   'gamma1' = as.vector(coef(res.b$modelD)[1]),
                                                   'gamma2' = as.vector(coef(res.b$modelD)[2]),
                                                   'se_beta1' = sqrt(diag(res.b$modelR$var))[1],
                                                   'se_beta2' = sqrt(diag(res.b$modelR$var))[2],
                                                   'se_gamma1' = sqrt(diag(res.b$modelD$var))[1],
                                                   'se_gamma2' = sqrt(diag(res.b$modelD$var))[2],
                                                   'K' = res.b$K,
                                                   'collapse' = res.b$collapse,
                                                   'AIC' = res.b$AIC,
                                                   'n.iter' = res.b$n.iter
                                  ) )
    }
  }
  
  return(coef.res)
}

jmdf.frail.res <- function(res){
  
  frail.res = NULL
  for(b in 1:length(res)){
    res.b = res[[b]]
    if(!inherits(res.b, "try-error") ){
      frail.res = rbind.data.frame(frail.res,
                                   cbind.data.frame('rep_b' = rep(b,res.b$K),
                                                    'K' = rep(res.b$K,res.b$K),
                                                    'group' = seq(1:res.b$K),
                                                    'w' = res.b$w,
                                                    'Pu' = res.b$P[,1],
                                                    'Pv' = res.b$P[,2],
                                                    'se_Pu' = res.b$se.P[,1],
                                                    'se_Pv' = res.b$se.P[,2]
                                   )
      )
    }
  }
  return(frail.res)
}

ng.jm.results <- function(res){
  results = NULL
  for(b in 1:length(res)){
    res.b = res[[b]]
    results = rbind.data.frame(results,
                               cbind.data.frame('rep_b' = b,
                                                'beta1' = res.b$Recurrent[1,1],
                                                'beta2' = res.b$Recurrent[2,1],
                                                'gamma1' = res.b$Death[1,1],
                                                'gamma2' = res.b$Death[2,1],
                                                'theta2_u' = res.b$Frailty[1,1], 
                                                'theta2_v' = res.b$Frailty[2,1], 
                                                'rho' = res.b$Frailty[3,1],
                                                'se_beta1' = res.b$Recurrent[1,2],
                                                'se_beta2' = res.b$Recurrent[2,2],
                                                'se_gamma1' = res.b$Death[1,2],
                                                'se_gamma2' = res.b$Death[2,2],
                                                'se_theta2_u' = res.b$Frailty[1,2], 
                                                'se_theta2_v' = res.b$Frailty[2,2], 
                                                'se_rho' = res.b$Frailty[3,2],
                                                'n.iter' = res.b$n.iter
                               )
    )
  } 
  return(results)
}

