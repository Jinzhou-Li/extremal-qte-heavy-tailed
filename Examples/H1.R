# Model H1 in the paper
source("Functions/LoadFunctions.R")

## Function to sample from model
h1.fun <- function(n){
  X <- runif(n)
  U <- runif(n)
  D <- 1*(0.5*X^2 + 0.25 >= U)
  S <- rt(n, df=3)
  Y1 <- 5*S*(1+X)
  Y0 <- S*(1+X)

  res <- list(Y1=Y1, Y0=Y0, D=D, X=X)
  return(res)
}

##
n <- 2000
CI_level=0.9

pn <- 1/n
ks <- n^(0.65)

## True QTE
# h1.true.qtes <- calculate_true_qtes(h1.fun, pn, Nsim=1e8)

## Sample data
sim_sample <- h1.fun(n)
Y1 <- sim_sample$Y1
Y0 <- sim_sample$Y0
D <- sim_sample$D
X <- sim_sample$X
Y <- D*Y1 + (1-D)*Y0

##
result_hill <- qte_extrapolation_hill(Y, X, D, pn, ks, CI_level)
