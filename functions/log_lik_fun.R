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

