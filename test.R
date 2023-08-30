# Why censoring ahaz can estimate the effect concisely even if the outcome model is incorrect
# This test aims to give a deep insight in the censoring, how it can bias the estimators

# Before the test, I do know that the commonly used assumption, the censoring independent.
# In my opinion, I believe the time to switch is not independent with Z conditional on L.
# However, all simulations indicates that the estimator for treatment effects are not biased.
library(ahaz)
library(MASS)
library(nleqslv)
expit = function(d) return(1/(1 + exp(-d)))


Simulation = function(beta, alpha, nrep, n, p, F_time, Censor, max_t = 3){
  # F_time: t, Z, Cov
  out = NULL
  Z_table = NULL
  event_table = vector(length = nrep)
  sigma = matrix(c(0.25, -1/6, -1/6, 0.25), 2, 2)
  correlation = vector(length = nrep)
  for (i in 1:nrep) {
    Cov = matrix(runif(p*n, 0, 1), nrow = n, ncol = p)
    U = mvrnorm(n, c(1.5, 1.5), sigma)
    gamma = rep(0.5, p)
    Z_p = expit(Cov %*% gamma)
    Z = rbinom(n, 1, Z_p)
    Z_table = cbind(Z_table, table(Z))
    time = F_time(beta, Z, Cov, n, U[, 1])
    C = Censor(alpha, Z, Cov, n, U[, 2])
    event = ifelse(time < C, 1, 0)
    time_c = ifelse(event, time, C)
    time_c = ifelse(time_c <= max_t, time_c, max_t)
    event_table[i] = mean(event)
    t = 1
    correlation[i] = mean((C >= t) * Z * (time >= t) * 0.2 * U[, 1])
    # idx = as.logical(event)
    surv = Surv(time= time_c + runif(n, 0, 0.0001), event = event, type = "right")
    m = ahaz(surv, cbind(Z, Cov))
    out = cbind(out, coef(m))
    cat("[[", "rep ", i, "\t", coef(m), "\t" , "]]\n",
        sep = " ")
  }
  return(list(out = out,
              Z_table = apply(Z_table, 1, mean),
              event_table = mean(event_table),
              correlation = mean(correlation)))
}



F_time = function(beta, Z, Cov, n, U1){
  return(rexp(n)/((0.2 + 0.1 * Z + Cov %*% beta + 0.2 * U1)))
}

 
Censor = function(alpha, Z, Cov, n, U2){
  F = function(t, alpha, Z, Cov, U2, uni){
    return(1 - exp(exp(-t)*(-0.025 - 0.1 * Z - (Cov %*% alpha) - 0.025 * U2) * t) - uni)
  }
  time = nleqslv(rep(0, n), F, alpha = alpha, Z = Z, Cov = Cov, U2 = U2, uni = runif(n))
  return(time$x)
}

Censor = function(alpha, Z, Cov, n, U2){
  return(rexp(n)/(0.025 + 0.1 * Z + (Cov %*% alpha) + 0.025 * U2))
}

nrep = 1000
n = 800
p = 2
out = Simulation(beta = rep(0.075, p), alpha = rep(0.025, p), nrep, n, p, F_time, Censor)
apply(out$out, 1, mean)
apply(out$out, 1, sd)
out$Z_table
out$event_table
out$correlation

