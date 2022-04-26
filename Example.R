library(ivsacim)
n = 400
event = rbinom(n, 1, 0.8)
IV = rbinom(n, 1, 0.5)
trt_init = IV
trt_shift = rep(0, n)
time = rexp(n)/(0.5 + trt_init * 0.2)
max_t = 3
max_t_bet = 3
n_sim = 0
fit <- ivsacim(time, event, IV, IV_valid = TRUE, trt_init, 
               trt_shift, covar = NULL, max_t, max_t_bet, n_sim)


max_time = max_t
weights <- rep(1, length(time))
event_new = event
event_new[time > max_time] = 0
stime = sort(time[event_new == 1])
weights <- weights/mean(weights)

iv_centered <- IV_center(instrument, covar)
expit <- function(x) {
  exp(x)/(1 + exp(x))
}
zmod <- glm(IV ~ 1, family = "binomial")
E_dot <- matrix(rep(1, length(IV)), ncol = 1) * c(fitted(zmod) * 
                                                   (1 - fitted(zmod)))
eps.theta <- as.matrix(iid(zmod))
p.dim <- dim(eps.theta)[2]
Z.c <- IV - fitted(zmod)
iv_centered <- list(Zc = Z.c, epstheta = eps.theta, Edot = E_dot, 
                  pdim = p.dim)


Zc <- Z.c
epstheta1 <- iv_centered$epstheta
Edot <- E_dot
pdim <- iv_centered$pdim
n <- length(time)
k <- length(stime)

D_status <- treatment_status(n = n, k = k, stime = stime, 
                             treatment_init = trt_init, treatment_shift_time = trt_shift, 
                             max_time = max_t)
D_status

