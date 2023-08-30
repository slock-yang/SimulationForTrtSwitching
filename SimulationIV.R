# 探究 IV 相关性问题

source("Function/Fun.R")
source("Function/Fun_parallel.R")
source("Function/integral.R")
library(ivsacim)
library(nleqslv)
library(jsonlite)
library(MASS)
library(ahaz)
# set.seed(1234)
expit = function(d) return(exp(d)/(exp(d)+1))
sigma = matrix(c(0.25, -1/6, -1/6, 0.25), 2, 2)


n = 800
gamma = 20
r = 0.5
L = runif(n, 0, 1)
L = matrix(L, nrow = n)
U = runif(n, 0, 1)
# x = mvrnorm(n = n, rep(2, 2), sigma)
# L = x[, 1, drop = F]
# U = x[, 2]
Z_p = expit((r*L))
Z = rbinom(n, 1, Z_p)
table(Z)
max_t = 6
W = rexp(n)/((0.1 + 0.01*L + 0.1 * gamma* Z +0.01*U))
cor(Z, W > 2)
F_time = function(t, w, z, l, u, uni) {
  ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - exp(0.075*l)*t - exp(0.075*u)*t)-uni, 
         1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - exp(0.05*l)*t - exp(0.075*u)*t) - uni)
}

time = nleqslv(rep(0, n), F_time, w = W, z = Z, l = L, u = U, uni = runif(n))$x
hist(time)
hist(W)
sum(time > W) / n
sum(ifelse(Z == 1, time > W, 0))/800
sum(ifelse(Z == 0, time > W, 0))/800

time_c = ifelse(time <= C, time, C)
time_c = ifelse(time_c <= max_t, time_c, max_t)
event = ifelse(time<=C & time<=max_t, 1, 0)
sum(time_c > W)/800
cor(Z, W)

# Generation --------------------------------------------------------------

# F_time = function(t, w, z, l, u, uni) {
#   ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - exp(0.075*l + 0.075*u)*t)-uni,
#          1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - exp(0.075*l
#                                                            + 0.075*u)*t) - uni)
# }

F_time = function(t, w, z, l, u, uni) {
  ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - exp(0.075*l)*t - exp(0.075*u)*t)-uni, 
         1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - exp(0.075*l)*t - exp(0.075*u)*t) - uni)
}

n = 800
nrep = 100
for (i in 1:nrep) {
  gamma = 20
  r = 0.5
  L = runif(n, 0, 1)
  L = matrix(L, nrow = n)
  U = runif(n, 0, 1)
  # x = mvrnorm(n = n, rep(2, 2), sigma)
  # L = x[, 1, drop = F]
  # U = x[, 2]
  Z_p = expit((r*L))
  Z = rbinom(n, 1, Z_p)
  C = rexp(n)/(0.01 + 0.05*L)
  # max_t = 6
  # W = rexp(n)/((0.01 + 0.01*L + gamma*0.1*Z + 0.01*U))
  W = rexp(n)/((0.1 + 0.01*L + gamma * 0.1*Z +0.01*U))
  time = nleqslv(rep(0, n), F_time, w = W, z = Z, l = L, u = U, uni = runif(n))$x
  dat = cbind(L, U, Z, C, W, time)
  write_json(dat, paste0("DataGeneration/IV1", n, "/", i, ".json"))
  if(i %% 10 == 0) cat("iter", i, "\n")
}



# Estimation --------------------------------------------------------------

n = 800
nrep = 100
parameters_dr = NULL
max_t = 6
for(i in 1:nrep){
  dat = read_json(paste0("DataGeneration/IV1", n, "/", i, ".json"), 
                  simplifyVector = T)
  Cov = dat[, 1, drop = F]
  Z = dat[, 3]
  time = dat[, 6]
  # C = dat[, 4]
  W = dat[, 5]
  # time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time <= W, time, W)
  # time_c = ifelse(time_c <= max_t, time_c, max_t)
  event = ifelse(time <= W, 1, 0)
  surv = Surv(time = time_c + runif(n, 0, 0.001), event = event, type = "right")
  out = coef(ahaz(surv, cbind(Z, Cov)))
  parameters_dr = cbind(parameters_dr, out)
  cat("[[", "rep ", i, "\t", out, "\t" , "]]\n",
      sep = " ")
}
apply(parameters_dr, 1, mean)
apply(parameters_dr - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, mean)
apply(parameters_dr - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep - 1)))
apply(parameters_dr, 1, sd)



n = 800
nrep = 100
parameters_dr = NULL
max_t = 6
for(i in 1:nrep){
  dat = read_json(paste0("DataGeneration/IV1", n, "/", i, ".json"), 
                  simplifyVector = T)
  Cov = dat[, 1, drop = F]
  Z = dat[, 3]
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  # time_c = ifelse(time <= C, time, C)
  # time_c = ifelse(time_c <= max_t, time_c, max_t)
  # time_c = ceiling(time*10)/10
  time_c = time
  # event = ifelse(time<=C & time<=max_t, 1, 0)
  event = rep(1, n)
  stime = sort(time_c[as.logical(event)])
  stime = unique(stime)
  k = length(stime)
  D_status = treatment_status(n, k, stime, Z, W, max_t)
  s1 = integral_est_cpp(c(0.1,0.075), time = time_c, event = event, IV = Z,
                        Covariates = Cov, Covariates2 = Cov, 
                        D_status = D_status, stime = stime)
  parameters_dr = cbind(parameters_dr, s1$x)
  cat("[[", "rep ", i, "\t", s1$x, "\t" ,s1$Convergence, "]]\n",
      sep = " ")
}
apply(parameters_dr, 1, mean)
apply(parameters_dr - matrix(c(0.1, 0.075), 
                             nrow = 2, ncol = nrep), 1, mean)
apply(parameters_dr - matrix(c(0.1, 0.075), 
                             nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(nrep - 1)))
apply(parameters_dr, 1, sd)

# explore -----------------------------------------------------------------

n = 800
nrep = 100
max_t = 6
out1 = vector(length = nrep)
out2 = vector(length = nrep)
out = vector(length = nrep)
for (i in 1:100) {
  dat = read_json(paste0("DataGeneration/IV1", n, "/", i, ".json"), 
                  simplifyVector = T)
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  Z = dat[, 3]
  # time_c = ifelse(time <= C, time, C)
  # time_c = ifelse(time_c <= max_t, time_c, max_t)
  # event = ifelse(time<=C & time<=max_t, 1, 0)
  out[i] = cor(Z, W > 1)
  out1[i] = sum(ifelse(Z == 1, time > W, 0))/n
  out2[i] = sum(ifelse(Z == 0, time > W, 0))/n
}
mean(out1)
mean(out2)
apply(out, 1, mean)

