init.min.gauss <- function(input.folder){
  # Define mu, Sigma, and L according to selected "input.folder"
  if(!(input.folder %in% c('I_L15','II_L1','II_L15','II_L2','III_L15'))){
    stop('Please define Init and L')
  }else{
    L.mindist = 1.5
    mu.0 = c(0,0)
    if(input.folder=='I_L15'){
      Sigma.0 = matrix(c(2,0.2,0.2,2), nrow = 2, ncol = 2)
    }else  if(input.folder=='III_L15'){
      Sigma.0 = matrix(c(0.6,0,0,0.2), nrow = 2, ncol = 2)  ## REV 2: theta_v = 0.2 --> 0.6; theta_v = 2 --> 0.2
    }else{
      Sigma.0 = matrix(c(0.2,0,0,2), nrow = 2, ncol = 2)     
      if(input.folder=='II_L1'){
        L.mindist = 1
      }else if(input.folder=='II_L2'){
        L.mindist = 2
      }
    }
  }
  return(list(mu.0, Sigma.0, L.mindist))
}

init.min.unif <- function(input.folder){
  # Define ulim.unif, vlim.unif, and L according to selected "input.folder"
  if(!(folder %in% c('I_L15','II_L1','II_L15','II_L2','III_L15'))){
    stop('Please define Init and L')
  }else{
    L.mindist = 1.5
    if(folder=='I_L15'){
      ulim.0 = c(-4,4)
      vlim.0 = c(-4,4)
    }else if(folder=='III_L15'){
      ulim.0 = c(-2,2)
      vlim.0 = c(-1.5,1.5)
    }else{
      ulim.0 = c(-1.5,1.5) 
      vlim.0 = c(-4,4)
      if(folder=='II_L1'){
        L.mindist = 1
      }else if(folder=='II_L2'){
        L.mindist = 2
      }
    }
  }
  return(list(ulim.0, vlim.0, L.mindist))
}