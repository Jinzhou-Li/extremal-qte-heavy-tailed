# Model H3 in the paper
library(actuar)
source("Functions/LoadFunctions.R")

## Function to sample from model
h3.fun <- function(n){
  
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  Y1 <- rpareto(n, shape=1.75+X, scale=2)
  Y0 <- rpareto(n, shape=1.75+5*X, scale=1) 

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

##
n <- 2000
Nsim <- 1000
N_bootstrap <- 10
CI_level=0.9

pn <- 1/n
ks <- n^(0.65)

## True QTE
# h3.true.qtes <- calculate_true_qtes(h3.fun, pn, Nsim=1e8)

## Sample data
sim_sample <- h3.fun(n)
Y1 <- sim_sample$Y1
Y0 <- sim_sample$Y0
D <- sim_sample$D
X <- sim_sample$X
Y <- D*Y1 + (1-D)*Y0

##
result_hill <- qte_extrapolation_hill(Y, X, D, pn, ks, CI_level)