## load libraries
library(survival)
library(MASS)
library(mvtnorm)

## Support Functions for log-likelihood computation
# Log-likelihood
LogLik<-function(Z,w,E_formula1,E_formula2,E_haz1,E_haz2){
  lw <-sum(Z%*%matrix(log(w)))
  lr <- sum(Z*(E_haz1-E_formula1))
  ld <- sum(Z*(E_haz2-E_formula2))
  return (lw+lr+ld)
}
# Classification log-likelihood
classLogLik<-function(Z,E_formula1,E_formula2,E_haz1,E_haz2){
  lr <- sum(Z*(E_haz1-E_formula1))
  ld <- sum(Z*(E_haz2-E_formula2))
  return (lr+ld)
}


## Joint model with discretely-distributed non-parametric frailty
################################################################################
## INPUT PARAMETERS:
##--------------
## - dataR:     Recurrent event data 
##              It must contain variables 'id', 'deltaR', 'deltaD', and 'gaptime'
## - formulaR:  Covariate formula for recurrent model ' ~ xR_1 +... + xR_p1 '
## - dataD:     Terminal event data
##              It must contain variables 'id', 'deltaR', 'deltaD', and 'gaptime'
## - formulaD:  Covariate formula for terminal model ' ~ xD_1 + ... + xD_p2 '
## - init.unif: Initialization: Uniform if TRUE (default); Gaussian if FALSE
## - distance:  Distance measure to be used. This must be one of "euclidean" (default), 
##              "maximum", "manhattan", "canberra", "binary" or "minkowski"
## - Sigma:     Variance matrix for Gaussian Initialization
## - mu:        Mean vector for Gaussian Initialization
## - ulim.unif: Limits for initial grid using Uniform Initialization (u-axis)
## - vlim.unif: Limits for initial grid using Uniform Initialization (v-axis) 
## - M:         Initial number of support points 
## - L:         Value of threshold for the merging process
## - max.it:    Maximum number of iterations (default is 100)  
## - toll:      Stopping threshold (default is 1e-3)
################################################################################
## OUTPUT: The JMdiscfrail() function returns a list with elements:
##----------------
## - $modelR:       Estimated recurrent event model
## - $modelD:       Estimated terminal event model
## - $K:            Number of mass points
## - $w:            Estimated weights
## - $P:            Estimated discrete random effect coordinates
## - $se.P:         Estimated Standard Errors for P
## - $id.subgroups: Assigned mass-point for each subject
## - $cumhazR:      Estimated baseline cumulative hazard for recurrent events
## - $cumhazD:      Estimated baseline cumulative hazard for the terminal event
## - $LogL:         Log-likelihood
## - $classLogL:    Classification log-likelihood
## - $AIC:          Akaike Information Criterion
## - $n.iter:       Number of iterations
## - $collapse:     TRUE if the algorithm has collapsed to 1-mass.
################################################################################


JMdiscfrail = function(dataR, formulaR, dataD, formulaD, 
                       init.unif = TRUE, distance = "euclidean", P.weighted = FALSE,
                       Sigma=NULL, mu=NULL, ulim.unif=NULL, vlim.unif=NULL, 
                       M, L, max.it = 100, toll = 1e-3, seed = 210197){
  #------------------------------------------------------------------------
  ## Auxiliary 
  
  # Model Matrices
  X1 <- model.matrix(as.formula(formulaR), data = dataR)[,-1]
  X2 <- model.matrix(as.formula(formulaD), data = dataD)[,-1]
  
  # Response vectors
  time1 <- dataR$gaptime
  time2 <- dataD$gaptime
  
  # Number of patients
  ID1 = dataR$id
  ID2 = dataD$id
  N = length(unique(ID1))
  
  # Number of covariates
  ncov1 = ncol(X1)
  ncov2 = ncol(X2)
  
  # Number of events per patient   
  D1 = table(dataR$deltaR,dataR$id)[2,]
  D2 = table(dataD$deltaD,dataD$id)[2,]
  
  # Encode ID as numeric 
  groups1 <- match(ID1, unique(ID1))
  groups2 <- match(ID2, unique(ID2))
  
  # Flag for 1-mass collapse
  collapse = FALSE
  
  #------------------------------------------------------------------------
  ## INITIALIZATION
  
  ## Grid Initialization
  set.seed(seed)
  if(!init.unif){
    K <- M
    if(is.null(Sigma) | is.null(mu)){
      stop('Please define Sigma matrix and/or mu vector for Gaussian initialization')
    }
    P <- mvrnorm(K,mu,Sigma)
    w<-dmvnorm(P,mu,Sigma)
    w<-w/sum(w)
    P_show<-P
    P1mean = rep(w%*%P[,1],length(P[,1]))
    P2mean = rep(w%*%P[,2],length(P[,2]))
    P_show[,1]<-P[,1]-as.vector(P1mean)
    P_show[,2]<-P[,2]-as.vector(P2mean)
  } else {
    if(is.null(ulim.unif) | is.null(vlim.unif)){
      stop('Please define grid limits for u-axis and/or v-axis in Uniform initialization')
    }
    K <- M
    temp_u <- seq(ulim.unif[1], ulim.unif[2], length.out = round(sqrt(K),0))
    temp_v <- seq(vlim.unif[1], vlim.unif[2], length.out = round(sqrt(K),0))
    P     <- unname(data.matrix(expand.grid(temp_u,temp_v)))
    w     <- rep(1/dim(P)[1],dim(P)[1])
    P_show<-P
    P1mean = rep(w%*%P[,1],length(P[,1]))
    P2mean = rep(w%*%P[,2],length(P[,2]))
    P_show[,1]<-P[,1]-as.vector(P1mean)
    P_show[,2]<-P[,2]-as.vector(P2mean)
    K<-dim(P)[1]
  }
  
  # Initial Grid Shrinkage
  is_near<-TRUE
  while(is_near){
    D<-dist(P_show, method = distance)
    D<-as.matrix(D)
    D[upper.tri(D)]<-10
    diag(D)<-10
    out<-which(D == min(D), arr.ind = TRUE)
    if(D[out][1]<(L/3)){
      #merge
      P_show[out[1,2],]=(P_show[out[1,2],]+P_show[out[1,1],])/2
      P_show<-P_show[-out[1,1],]
      #update weights
      w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
      w<-w[-out[1,1]]
      w<-w/sum(w)
      K<-K-1
    }
    else {
      is_near=FALSE
    }
  }
  
  # Assign patient to random frailty, built frailty vectors
  P_index <- sample(1:K,size=N,replace = T, prob = w)
  P_off1  <- P_show[,1][P_index[groups1]]
  P_off2  <- P_show[,2][P_index[groups2]]
  
  
  ## Parameter Initialization
  # Estimate initial recurrent model, cumulative hazard and hazard
  formula1 = as.formula(paste0('Surv(time1,deltaR)',formulaR,'+ offset(P_off1)'))
  cox1 <- coxph(formula1, data=dataR)
  beta<-cox1$coefficients
  # cumulative hazard
  s1   <- survfit(cox1,data=dataR)
  cumhaz1 = cbind.data.frame( 'time' = s1$time, 'cumhaz' = s1$cumhaz)
  # hazard
  haz1 = cbind.data.frame( 'time' = s1$time, 'hazard' = diff(c(0,cumhaz1$cumhaz))/diff(c(0,s1$time)) )
  haz1$hazard = ifelse(haz1$hazard==0,1e-200,haz1$hazard)
  
  # Estimate initial terminal model, cumulative hazard and hazard
  formula2 = as.formula(paste0('Surv(time2,deltaD)',formulaD,'+ offset(P_off2)'))
  cox2 <- coxph(formula2, data=dataD)
  gamma <-cox2$coefficients
  # cumulative hazard
  s2   <- survfit(cox2,data=dataD)
  cumhaz2 = cbind.data.frame( 'time' = s2$time, 'cumhaz' = s2$cumhaz)
  # hazard
  haz2 = cbind.data.frame( 'time' = s2$time, 'hazard' = diff(c(0,cumhaz2$cumhaz))/diff(c(0,s2$time)) )
  haz2$hazard = ifelse(haz2$hazard==0,1e-200,haz2$hazard)
  
  #------------------------------------------------------------------------
  ## EM ALGORITHM -- ITERATIONS
  # Initialize structures for computations
  numerator <- rep( 0, K )
  Z <- E_formula1 <- E_formula2 <- E_haz1 <- E_haz2 <- matrix( 0, nrow = N, ncol = K)
  E_part1 <- E_part2<-rep( 0, N)
  
  # Start loop
  converged = FALSE
  it <- 1
  while (!converged & it <= max.it ){
    
    # Save current estimates
    w_old <- w
    Z_old <- Z
    beta_old <- beta
    gamma_old <- gamma
    
    # Support Reduction I: Grid Shrinking
    # Grid Shrinking
    D<-dist(P_show, method = distance)
    D<-as.matrix(D)
    D[upper.tri(D)]<-10
    diag(D)<-10
    out<-which(D == min(D), arr.ind = TRUE)
    if(D[out][1]<L){
      #merge
      if(P.weighted){
        P_show[out[1,2],]=(P_show[out[1,2],]*w[out[1,2]]+P_show[out[1,1],]*w[out[1,1]])/(w[out[1,2]]+w[out[1,1]])
      }else{
        P_show[out[1,2],]=(P_show[out[1,2],]+P_show[out[1,1],])/2
      }
      P_show<-P_show[-out[1,1],]
      #update weights
      w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
      w<-w[-out[1,1]]
      w<-w/sum(w)
      K<-K-1
    }
    
    if(K==1){
      collapse = TRUE
      cat("Warning: The discrete frailty component has collapsed to a single 
          value, indicating no observed heterogeneity among patients for the joint 
          model. Consider increasing the distance threshold L or fitting separate 
          models for recurrent and terminal events.")
      break
    } 
    
    # Clean Structures
    Z <- E_formula1 <- E_formula2 <- E_haz1 <- E_haz2 <- matrix( 0, nrow = N, ncol = K)
    numerator <- rep(0,K)
    
    # Expectation Step
    for(i in 1:N){
      
      current_patient1 <- groups1==i
      current_patient2 <- groups2==i
      
      ebz1 <- exp( X1[current_patient1,] %*% beta )
      ebz2 <- exp( X2[current_patient2,] %*% gamma)
      
      tRij <- match(time1[current_patient1], cumhaz1$time)
      if(sum(is.na(tRij))>0){
        na.ind = which(is.na(tRij))
        for(replace.na in na.ind){
          tRij[replace.na] = which.min(abs(cumhaz1$time - time1[current_patient1][replace.na]))
        }}
      H01t <- cumhaz1$cumhaz[tRij]
      lh01t<-  log(haz1$hazard[tRij])
      
      tDi <- match(time2[current_patient2], cumhaz2$time)
      tDi = ifelse(is.na(tDi), which.min(abs(cumhaz2$time - time2[current_patient2])), tDi)
      H02t <- cumhaz2$cumhaz[tDi]
      lh02t<-  log(haz2$hazard[tDi])
      
      E_part1[i] <- ifelse( ncov1 > 0,
                            sum( H01t*ebz1 ),
                            sum( H01t ) )
      E_part2[i] <- ifelse( ncov2 > 0,
                            H02t*ebz2 ,
                            H02t)
      
      for(l in 1:K){
        
        E_formula1[i,l] <- ifelse( ncov1 > 0,
                                   sum( H01t*ebz1*exp(P_show[l,1]+P1mean[l])),
                                   sum( H01t*exp(P_show[l,1]+P1mean[l]) ) )
        E_formula2[i,l] <- ifelse( ncov2 > 0,
                                   H02t*ebz2*exp(P_show[l,2]+P2mean[l]),
                                   H02t*exp(P_show[l,2]+P2mean[l]))
        E_haz1[i,l]     <- sum(dataR$deltaR[current_patient1]*(lh01t+log(ebz1)+P_show[l,1]+P1mean[l]))
        
        E_haz2[i,l]     <- dataD$deltaD[current_patient2]*(lh02t+log(ebz2)+P_show[l,2]+P2mean[l])
        
        pivot <- min(as.numeric(D1)[i]*(P_show[l,1]+P1mean[l]) - E_formula1[i,l] +
                       + as.numeric(D2)[i]*(P_show[l,2]+P2mean[l]) - E_formula2[i,l])
        numerator[l] <- w[l]*exp(as.numeric(D1)[i]*(P_show[l,1]+P1mean[l]) - E_formula1[i,l] +
                                   + as.numeric(D2)[i]*(P_show[l,2]+P2mean[l]) - E_formula2[i,l])
        
        if(max(numerator)==0)
          numerator <- 1e-16/((as.numeric(D1)[i]*(P_show[,1]+P1mean)-E_formula1[i,]+as.numeric(D2)[i]*(P_show[,2]+P2mean)-E_formula2[i,])/pivot)
      }
      
      Z[i,] <- numerator/sum(numerator)
      
    }
    
    # Maximization Step
    # Latent partition
    P_group <- as.numeric( apply(Z, 1, which.max) )
    
    # Support Reduction II: Grid Shrinking - Unassigned points
    t<-table(factor(P_group, levels = 1:K))
    to_elim<-which(as.numeric(t)==0)
    if(length(to_elim)>0){
      Z<-Z[,-to_elim]
      E_formula1<-E_formula1[,-to_elim]
      E_formula2<-E_formula2[,-to_elim]
      E_haz1<-E_haz1[,-to_elim]
      E_haz2<-E_haz2[,-to_elim]
      numerator <-numerator[-to_elim]
      P_show<-P_show[-to_elim,]
      K<-K-length(to_elim)
      if(K>1){ Z = Z/rowSums(Z) }
    }
    
    # Vector of proportions
    if(K>1){
      w <- (colSums(Z))/ N
    }else{
      collapse = TRUE
      cat("Warning: The discrete frailty component has collapsed to a single value, 
      indicating no observed heterogeneity among patients for the joint model.
      Consider increasing the distance threshold L or fitting separate models 
      for recurrent and terminal events.")
      w <- 1 
      P <- matrix(c(0,0),nrow = 1,ncol = 2)
      message("Exited all loops")
      break
    }
    
    ## Unconstrained Optimization
    P = P_show # to set matrix dimension
    P[,1]   <- log(( as.numeric(D1) %*% Z )/( E_part1 %*% Z))
    P[,2]   <- log(( as.numeric(D2) %*% Z )/( E_part2 %*% Z))
    
    # Standard error computation
    SE_P <- matrix(NA, nrow=K, ncol=2)
    hess1 <- - ( E_part1 %*% Z)*exp(P[,1])
    hess2 <- - ( E_part2 %*% Z)*exp(P[,2])
    SE_P[,1] <- 1/(-hess1)
    SE_P[,2] <- 1/(-hess2)
    
    # Centering mass points
    P1mean = rep(w%*%P[,1],length(P[,1]))
    P2mean = rep(w%*%P[,2],length(P[,2]))
    P_show[,1]<-P[,1]-as.vector(P1mean)
    P_show[,2]<-P[,2]-as.vector(P2mean)
    
    # Frailties update
    P_off1 <-log((Z%*%exp(matrix(P_show[,1]+P1mean)))[groups1])
    P_off2 <-log((Z%*%exp(matrix(P_show[,2]+P2mean)))[groups2])
    
    # Estimate Betas
    temp_model1 <- coxph(formula1, data=dataR, method = "breslow")
    beta       <- temp_model1$coef
    
    # Estimate Gammas
    temp_model2 <- coxph(formula2, data=dataD, method = "breslow")
    gamma       <- temp_model2$coef 
    
    # Estimate Cumulative Hazard and Hazard functions - Recurrent
    s1   <- survfit(temp_model1,data=dataR)
    cumhaz1$cumhaz = s1$cumhaz
    haz1$hazard = diff(c(0,cumhaz1$cumhaz))/diff(c(0,s1$time))
    haz1$hazard = ifelse(haz1$hazard==0,1e-200,haz1$hazard)
    
    # Estimate Cumulative Hazard and Hazard functions - Terminal
    s2   <- survfit(temp_model2,data=dataD)
    cumhaz2$cumhaz = s2$cumhaz
    haz2$hazard = diff(c(0,cumhaz2$cumhaz))/diff(c(0,s2$time))
    haz2$hazard = ifelse(haz2$hazard==0,1e-200,haz2$hazard) ####!!!!
    
    # Check for convergence
    if(length(w)==length(w_old)){
      eps = max(abs(w - w_old),na.rm=T) + max(abs(gamma - gamma_old)) + max(abs(beta - beta_old)) 
      converged = ifelse(eps < toll, TRUE, FALSE)
    } 
    
    # AIC computation
    LogL = LogLik(Z = Z, w = w, E_formula1 = E_formula1, E_formula2 = E_formula2,
                  E_haz1 = E_haz1, E_haz2 = E_haz2)
    classLogL = classLogLik(Z = Z, E_formula1 = E_formula1 ,E_formula2 = E_formula2,
                            E_haz1 = E_haz1, E_haz2 = E_haz2)
    AIC = 2*( 2*K+(K-1)+length(beta)+length(gamma) ) - 2*classLogL
    
    # Rescale baseline cumhaz by mean P
    cumhaz1_show <- cumhaz1
    cumhaz2_show <- cumhaz2
    cumhaz1_show$cumhaz <- cumhaz1$cumhaz*exp(P1mean[1])
    cumhaz2_show$cumhaz <- cumhaz2$cumhaz*exp(P2mean[1])
    
    # Save 
    temp_list <- list("modelR" = temp_model1, "modelD" = temp_model2,
                      "K" = K, "w" = w, "P" = P_show, "se.P" = SE_P,
                      "id.subgroups" = cbind.data.frame('id' = ID2, 'subgroup' = P_group),
                      "cumhazR" = cumhaz1_show, "cumhazD" = cumhaz2_show,
                      "LogL" = LogL, "classLogL" = classLogL, "AIC" = AIC,
                      "n.iter" = it)
    
    # Update iteration
    it <- it + 1
    
  }
  
  temp_list$collapse = collapse
  return(temp_list)
  
}

