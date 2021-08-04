## Estimating the extreme QTE and its CI using the Pickands type EVI estimator. 

qte_extrapolation_pickands <- function(Y, X=NULL, D, pn, ks, CI_level=0.9, prop_scores=NULL){
  
  n <- length(Y)
  
  qtes <- numeric(length(ks))           # For storing resulting qte
  q1ext <- numeric(length(ks))          # For storing extrapolated quantile q1(1-pn)
  q0ext <- numeric(length(ks))          # For storing extrapolated quantile q0(1-pn)
  p1s <- numeric(length(ks))            # For storing estimated quantile index gamma1
  p0s <- numeric(length(ks))            # For storing estimated quantile index gamma0
  variances <- numeric(length(ks))      # For storing estimated variance
  CIupper <- numeric(length(ks))        # For storing the upper bound of CI
  CIlower <- numeric(length(ks))        # For storing the lower bound of CI
  
  if(is.null(prop_scores)){ 
    ## Fit polynomial sieve
    hn <- 2*n^(1/11)
    prop.fit <- glm(D ~ poly(X, hn), family=binomial)
    prop_scores <- fitted(prop.fit)
  }
  
  ## Inverse propensity score weights
  w1 <- D / prop_scores
  w0 <- (1-D) / (1-prop_scores) 
  
  vw1 <- D/prop_scores^2                ## Weights for variance estimaton
  vw0 <- (1-D)/(1-prop_scores)^2
  
  ## Estimate intermediate counterfactual quantiles
  q1s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w1, tau=c(1-ks/n, 1-2*ks/n, 1-4*ks/n))
  q0s <- BMisc::getWeightedQuantiles(cvec=Y, weights=w0, tau=c(1-ks/n, 1-2*ks/n, 1-4*ks/n)) 
  
  ## Loop over different choices of k
  for (j in 1:length(ks)){
    k <- ks[j] 
    
    ## Extract needed intermediate quantiles
    q1 <- c(q1s[j], q1s[j+length(ks)], q1s[j+2*length(ks)])
    q0 <- c(q0s[j], q0s[j+length(ks)], q0s[j+2*length(ks)])
    
    # Pickands type estimates of the extreme value indices
    p1 <- log((q1[1]-q1[2])/(q1[2]-q1[3]))/log(2)  
    p0 <- log((q0[1]-q0[2])/(q0[2]-q0[3]))/log(2) 
    
    ## Quantile extrapolation
    dn <- k/n/pn
    qe1 <- q1[1]*dn^p1
    qe0 <- q0[1]*dn^p0
    
    p1s[j] <- p1
    p0s[j] <- p0
    qtes[j] <- qe1-qe0
    q1ext[j] <- qe1
    q0ext[j] <- qe0
    
    ########## Variance estimation #########
    ## Covariance estimation
    H1 <- matrix(nrow=3, ncol=3)
    H0 <- matrix(nrow=3, ncol=3)
    
    for (i in 1:3){
      H1[i, i] <- sum(vw1*(Y>q1[i]))/k 
      H0[i, i] <- sum(vw0*(Y>q0[i]))/k
    }
    
    H1[1,] <- H1[1,1]
    H1[,1] <- H1[1,1]
    H0[1,] <- H0[1,1]
    H0[,1] <- H0[1,1]
    H1[2,3] <- H1[2,2]
    H1[3,2] <- H1[2,2]
    H0[2,3] <- H0[2,2]
    H0[3,2] <- H0[2,2]
    
    ## Estimation of alpha
    alpha <- qe1/qe0
    
    ## Estimation of b
    b1 <- numeric(3)
    b0 <- numeric(3)
    b1[1] <- (2^(-p1)-1)/p1/log(2)
    b0[1] <- (2^(-p0)-1)/p0/log(2)
    b1[2] <- sqrt(2)*(2^p1 - 2^(-p1))/p1/log(2)  
    b0[2]<- sqrt(2)*(2^p0 - 2^(-p0))/p0/log(2)  
    b1[3]<- 2*(1 - 2^(p1))/p1/log(2)  
    b0[3]<- 2*(1 - 2^(p0))/p0/log(2)  
    
    ## Estimation of sigma^2 
    sigmasq <- 0
    for (i1 in 1:3){
      for (i2 in 1:3){
        sigmasq = sigmasq + (min(1, alpha^2)*b1[i1]*b1[i2]*H1[i1,i2] + min(1, 1/alpha^2)*b0[i1]*b0[i2]*H0[i1,i2])/sqrt(i1*i2)
      }
    }
    
    norm_factor <- sqrt(k)/(max(qe1, qe0))/log(dn)
    variances[j] <- sigmasq/norm_factor^2
    
    Gaussian_q <- qnorm(1-(1-CI_level)/2)
    CIupper[j] <- (qe1-qe0) + Gaussian_q * sqrt(sigmasq/norm_factor^2)
    CIlower[j] <- (qe1-qe0) - Gaussian_q * sqrt(sigmasq/norm_factor^2)
  }
  
  res <- list(qtes=qtes, q1ext=q1ext, q0ext=q0ext, p1s=p1s, p0s=p0s, sd=sqrt(variances), CIupper=CIupper, CIlower=CIlower)
  
  return(res)
}
