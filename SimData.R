simdata = function(n)
{
  r = 0.5
  L = runif(n, 0, 1)
  L = matrix(L, nrow = n)
  U = runif(n, 0, 1)
  epsilon = rnorm(n)
  Z_p = expit(r*L + epsilon)
  Z = rbinom(n, 1, Z_p)
  C = rexp(n)/(0.01 + 0.05*L)
  max_t = 6
  W = rexp(n)/(0.01 + 0.01*L + 0.1*Z)
  
  F_time = function(t, w, z, l, u, uni) {
    ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - 0.075*l*t- 0.075*u*t)-uni, 
           1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - 0.075*l*t
                   - 0.075*u*t) - uni)
  }
  time = nleqslv(rep(0, n), F_time, w = W, z = Z, l = L, u = U, uni = runif(n))$x
  time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time_c <= max_t, time_c, max_t)
  event = ifelse(time<=C & time<=max_t, 1, 0)
  stime = sort(time_c[as.logical(event)])
  k = length(stime)
  D_status = treatment_status(n, k, stime, Z, W, max_t)
  
  return(list(time = time_c, event = event, IV = Z, Covariates = L,
              D_status = D_status, stime = stime, U = U))
}