qte_extrapolation_bootstrap <- function(extrapolators, Y, X=NULL, D, pn, ks, CI_level=0.9, N_bootstrap, prop_scores=NULL){
  
  bootstrapped_qtes <- matrix(0, nrow=N_bootstrap, ncol=length(ks)) 
  
  if(extrapolators=="Hill"){
    
    qtes <- qte_extrapolation_hill(Y, X, D, pn, ks, CI_level, prop_scores)$qtes
    
    for (i in 1:N_bootstrap) {
      ## Sample bootstrap dataset
      selected <- sample(1:length(Y), replace=TRUE)
      Yboot <- Y[selected]
      Dboot <- D[selected]
      Xboot <- X[selected]
      
      bootstrapped_qtes[i,] <- qte_extrapolation_hill(Yboot, Xboot, Dboot, pn, ks, CI_level, prop_scores)$qtes
    }
  }
  
  if(extrapolators=="Pickands"){
    
    qtes <- qte_extrapolation_pickands(Y, X, D, pn, ks, CI_level, prop_scores)$qtes
    
    for (i in 1:N_bootstrap) {
      ## Sample bootstrap dataset
      selected <- sample(1:length(Y), replace=TRUE)
      Yboot <- Y[selected]
      Dboot <- D[selected]
      Xboot <- X[selected]
      
      bootstrapped_qtes[i,] <- qte_extrapolation_pickands(Yboot, Xboot, Dboot, pn, ks, CI_level, prop_scores)$qtes
    }
  }
  
  ## Derive bootstrap confidence intervals
  Gaussian_q <- qnorm(1-(1-CI_level)/2)
  sd_boot <- apply(bootstrapped_qtes, 2, sd)
  CIupper <- qtes + Gaussian_q * sd_boot
  CIlower <- qtes - Gaussian_q * sd_boot
  
  res <- list(qtes=qtes, CIupper=CIupper, CIlower=CIlower, sd_boot=sd_boot)
  
  return(res)
}
