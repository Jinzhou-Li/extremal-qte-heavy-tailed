# Model H2 in the paper
library(evd)
source("Functions/LoadFunctions.R")

## Function to sample from model
h2.fun <- function(n){
  
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  C1 <- rfrechet(n, shape=2)
  C2 <- rfrechet(n, shape=3)
  Y1 <- C1*exp(X)
  Y0 <- C2*exp(X)

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

##
n <- 2000
CI_level=0.9

pn <- 1/n
ks <- n^(0.65)

## True QTE
# h2.true.qtes <- calculate_true_qtes(h2.fun, pn, Nsim=1e8)

## Sample data
sim_sample <- h2.fun(n)
Y1 <- sim_sample$Y1
Y0 <- sim_sample$Y0
D <- sim_sample$D
X <- sim_sample$X
Y <- D*Y1 + (1-D)*Y0

##
result_hill <- qte_extrapolation_hill(Y, X, D, pn, ks, CI_level)