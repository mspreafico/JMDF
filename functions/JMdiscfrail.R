## load libraries
library(survival)
library(MASS)
library(mvtnorm)
## load support functions for log-likelihood computations
source('functions/log_lik_fun.R')

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
## - centering: Centering cumulative baseline hazards and mass-points around 
##              (u,v)=(0,0). This helps in interpreting the results.
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
################################################################################


JMdiscfrail = function(dataR, formulaR, dataD, formulaD, 
                       init.unif = TRUE, distance = "euclidean", 
                       Sigma=NULL, mu=NULL, ulim.unif=NULL, vlim.unif=NULL, 
                       M, L, max.it = 100, toll = 1e-3, centering=TRUE){
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
  
  # Number of unique times
  nt1 = length(unique(dataR$gaptime))
  nt2 = length(unique(dataD$gaptime))
  
  # Number of observations
  nobs1 = nrow(dataR)
  nobs2 = nrow(dataD)
  
  # Number of covariates
  ncov1 = ncol(X1)
  ncov2 = ncol(X2)
  
  # Cumulative Hazards: data frame encoding
  cumhaz1 = as.data.frame( cbind( time = sort( unique( time1 )), cumhaz=rep( 0, nt1 )))
  cumhaz2 = as.data.frame( cbind( time = sort( unique( time2 )), cumhaz=rep( 0, nt2 )))
  
  # Instantaneuos hazards
  haz1 = as.data.frame( cbind( time = sort( unique( time1 )), hazard=rep( 0, nt1 )))
  haz2 = as.data.frame( cbind( time = sort( unique( time2 )), hazard=rep( 0, nt2 )))
  
  # Number of events per patient   
  D1 <- table( ID1[ dataR$deltaR == 1 ] )
  orderD1<-match(unique(ID1),unique(ID1)[order(unique(ID1))])
  D1<- D1[orderD1]
  D2 <- table( ID2[dataD$deltaD == 1])
  orderD2<-match(unique(ID1),unique(ID2)[order(unique(ID2))])
  D2<-D2[orderD2]
  
  # Risk sets
  risk_index1 <- matrix( 0, nrow = nobs1, ncol = nt1 )
  risk_index2 <- matrix(0, nrow = nobs2, ncol = nt2)
  
  # Time lists
  time_list1 <- sapply( 1:nt1, function(x) !is.na(match(time1, cumhaz1$time[x])))
  time_list2 <- sapply( 1:nt2, function(x) !is.na(match(time2, cumhaz2$time[x])))
  
  # Number of ties
  m1 <- sapply( 1:dim(time_list1)[2], function(x) sum(dataR$deltaR[time_list1[,x]]))
  m2 <- sapply( 1:dim(time_list2)[2], function(x) sum(dataD$deltaD[time_list2[,x]]))
  
  # fill risk index 
  for( l in 1:nt1 ){
    risk_index1[ which( time1 >= cumhaz1$time[ l ]), l ] <-  1 
  }
  for( l in 1:nt2 ){
    risk_index2[ which( time2 >= cumhaz2$time[ l ]), l ] <-  1 
  }
  
  # Encode ID as numeric 
  groups1 <- match(ID1, unique(ID1))
  groups2 <- match(ID2, unique(ID1))
  
  #------------------------------------------------------------------------
  ## INITIALIZATION
  
  ## Grid Initialization
  set.seed(210197)
  if(!init.unif){
    K <- M
    if(is.null(Sigma) | is.null(mu)){
      stop('Please define Sigma matrix and/or mu vector for Gaussian initialization')
    }
    P <- mvrnorm(K,mu,Sigma)
    w<-dmvnorm(P,mu,Sigma)
    w<-w/sum(w)
    P_show<-P
    P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
    P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))
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
    P_show[,1]<-P[,1]-rep(w%*%P[,1],length(P[,1]))
    P_show[,2]<-P[,2]-rep(w%*%P[,2],length(P[,2]))
    K<-dim(P)[1]
  }
  
  # Initial Grid Shrinkage
  is_near<-TRUE
  while(is_near){
    D<-dist(P, method = distance)
    D<-as.matrix(D)
    D[upper.tri(D)]<-10
    diag(D)<-10
    out<-which(D == min(D), arr.ind = TRUE)
    if(D[out][1]<(L/3)){
      #merge
      P[out[1,2],]=(P[out[1,2],]+P[out[1,1],])/2
      P<-P[-out[1,1],]
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
  P_off1  <- P[,1][P_index[groups1]]
  P_off2  <- P[,2][P_index[groups2]]
  
  
  ## Parameter Initialization
  # Estimate initial recurrent model, cumulative hazard and hazard
  formula1 = as.formula(paste0('Surv(time1,deltaR)',formulaR,'+ offset(P_off1)'))
  cox1 <- coxph(formula1, data=dataR)
  beta<-cox1$coefficients
  s1   <- survfit(cox1,data=dataR)
  cumhaz1$cumhaz = s1$cumhaz
  haz1$hazard = diff(c(0,cumhaz1$cumhaz)) 
  for(j in 1:length(haz1$hazard)){
    if(haz1$hazard[j]==0)
      haz1$hazard[j]<-haz1$hazard[j-1]
  }
  
  # Estimate initial terminal model, cumulative hazard and hazard
  formula2 = as.formula(paste0('Surv(time2,deltaD)',formulaD,'+ offset(P_off2)'))
  cox2 <- coxph(formula2, data=dataD)
  gamma <-cox2$coefficients
  s2   <- survfit(cox2,data=dataD)
  cumhaz2$cumhaz = s2$cumhaz
  haz2$hazard = diff(c(0,cumhaz2$cumhaz)) 
  for(j in 1:length(haz2$hazard)){
    if(haz2$hazard[j]==0)
      haz2$hazard[j]<-haz2$hazard[j-1]
  }
  
  
  #------------------------------------------------------------------------
  ## EM ALGORITHM -- ITERATIONS
  # Initialize structures for computations
  numerator <- rep( 0, K )
  Z <- E_formula1 <- E_formula2 <- E_haz1<-E_haz2<- matrix( 0, nrow = N, ncol = K)
  E_part1 <- E_part2<-rep( 0, N)
  
  # Start loop
  converged = FALSE
  it <- 1
  while (!converged & it <= max.it ){
    
    # Save current estimates
    w_old <- w
    Z_old <- Z
    P_old <- P
    beta_old <- beta
    gamma_old <- gamma
    
    # Support Reduction I: Grid Shrinking
    # Grid Shrinking
    D<-dist(P, method = distance)
    D<-as.matrix(D)
    D[upper.tri(D)]<-10
    diag(D)<-10
    out<-which(D == min(D), arr.ind = TRUE)
    if(D[out][1]<L){
      #merge
      P[out[1,2],]=(P[out[1,2],]+P[out[1,1],])/2
      P<-P[-out[1,1],]
      #update weights
      w[out[1,2]]<-w[out[1,2]]+w[out[1,1]]
      w<-w[-out[1,1]]
      w<-w/sum(w)
      K<-K-1
    }
    
    # Clean Structures
    Z <- E_formula1 <- E_formula2 <- E_haz1 <- E_haz2<- matrix( 0, nrow = N, ncol = K)
    numerator <- rep(0,K)
    
    # Expectation Step
    for(i in 1:N){
      
      current_patient1 <- groups1==i
      current_patient2 <- groups2==i
      
      ebz1 <- exp( X1[current_patient1,] %*% beta )
      ebz2 <- exp( X2[current_patient2,] %*% gamma)
      
      tRij <- match(time1[current_patient1], cumhaz1$time)
      H01t <- cumhaz1$cumhaz[tRij]
      lh01t<-  log(haz1$hazard[tRij])
      
      tDi <- match(time2[current_patient2], cumhaz2$time)
      H02t <- cumhaz2$cumhaz[tDi]
      lh02t<-  log(haz2$hazard[tDi])
      
      E_part1[i] <- ifelse( ncov1 > 0,
                            sum( H01t*ebz1 ),
                            sum( H01t ) )
      E_part2[i] <- ifelse( ncov2 > 0,
                            H02t*ebz2 ,
                            H02t)
      
      for(l in 1:K){
        if(K==1){P=matrix(as.vector(P), nrow=1, ncol=2)}
        E_formula1[i,l] <- ifelse( ncov1 > 0,
                                   sum( H01t*ebz1*exp(P[l,1])),
                                   sum( H01t*exp(P[l,1]) ) )
        E_formula2[i,l] <- ifelse( ncov2 > 0,
                                   H02t*ebz2*exp(P[l,2]),
                                   H02t*exp(P[l,2]))
        E_haz1[i,l]     <- sum(dataR$deltaR[current_patient1]*(lh01t+log(ebz1)+P[l,1]))
        
        E_haz2[i,l]     <- dataD$deltaD[current_patient2]*(lh02t+log(ebz2)+P[l,2])
        
        pivot <- min(as.numeric(D1)[i]*(P[l,1]) - E_formula1[i,l] +
                       + as.numeric(D2)[i]*(P[l,2]) - E_formula2[i,l])
        numerator[l] <- w[l]*exp(as.numeric(D1)[i]*(P[l,1]) - E_formula1[i,l] +
                                   + as.numeric(D2)[i]*(P[l,2]) - E_formula2[i,l])
        
        if(max(numerator)==0)
          numerator <- 1e-16/((as.numeric(D1)[i]*P[,1]-E_formula1[i,]+as.numeric(D2)[i]*P[,2]-E_formula2[i,])/pivot)
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
      P<-P[-to_elim,]
      K<-K-length(to_elim)
      Z = Z/rowSums(Z)
    }
    
    # Vector of proportions
    if(K>1){
      w <- (colSums(Z))/ N
    }else{
      w <- 1 
      P <- matrix(c(0,0),nrow = 1,ncol = 2)
    }
    
    ## Unconstrained Optimization
    P[,1]   <- log(( as.numeric(D1) %*% Z )/( E_part1 %*% Z))
    P[,2]   <- log(( as.numeric(D2) %*% Z )/( E_part2 %*% Z))
    
    # Standard error computation
    SE_P <- matrix(NA, nrow=K, ncol=2)
    hess1 <- - ( E_part1 %*% Z)*exp(P[,1])
    hess2 <- - ( E_part2 %*% Z)*exp(P[,2])
    SE_P[,1] <- 1/(-hess1)
    SE_P[,2] <- 1/(-hess2)
    
    # Frailties update
    P_off1 <-log((Z%*%exp(matrix(P[,1])))[groups1])
    P_off2 <-log((Z%*%exp(matrix(P[,2])))[groups2])
    
    # Estimate Betas
    temp_model1 <- coxph(formula1, data=dataR, method = "breslow")
    beta       <- temp_model1$coef
    
    # Estimate Gammas
    temp_model2 <- coxph(formula2, data=dataD, method = "breslow")
    gamma       <- temp_model2$coef 
    
    # Estimate Cumulative Hazard and Hazard functions - Recurrent
    s1   <- survfit(temp_model1,data=dataR)
    cumhaz1$cumhaz = s1$cumhaz
    haz1$hazard = diff(c(0,cumhaz1$cumhaz)) 
    for(j in 1:length(haz1$hazard)){
      if(haz1$hazard[j]==0)
        haz1$hazard[j]<-haz1$hazard[j-1]
    }
    
    # Estimate Cumulative Hazard and Hazard functions - Terminal
    s2   <- survfit(temp_model2,data=dataD)
    cumhaz2$cumhaz = s2$cumhaz
    haz2$hazard = diff(c(0,cumhaz2$cumhaz)) 
    for(j in 1:length(haz2$hazard)){
      if(haz2$hazard[j]==0)
        haz2$hazard[j]<-haz2$hazard[j-1]
    }
    
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
    
    # Centering baseline cumhaz and mass points for a preferable interpretation
    P_show    <- P
    cumhaz1_show <- cumhaz1
    cumhaz2_show <- cumhaz2
    if(centering){
      P_show[,1]<- P[,1] - as.vector(w %*% P[,1])
      P_show[,2]<- P[,2] - as.vector(w %*% P[,2])
      cumhaz1_show$cumhaz <- cumhaz1$cumhaz*exp(as.vector(w %*% P[,1]))
      cumhaz2_show$cumhaz <- cumhaz2$cumhaz*exp(as.vector(w %*% P[,2]))
    }
    
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
  
  return(temp_list)
  
}

