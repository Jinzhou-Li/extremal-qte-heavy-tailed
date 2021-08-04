calculate_true_qtes <- function(sampling.fun, pn, Nsim=1e8){
  ## Approximates the true QTE by sampling from the Model and calculating the quantiles
  set.seed(1)
  samples <- sampling.fun(Nsim)
  q1 = quantile(samples$Y1, probs = 1-pn)
  q0 = quantile(samples$Y0, probs = 1-pn)
  return(q1 - q0)
}