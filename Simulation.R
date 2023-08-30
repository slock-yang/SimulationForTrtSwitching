source("Function/Fun.R")
source("Function/Fun_parallel.R")
source("Function/integral.R")
library(ivsacim)
library(nleqslv)
library(jsonlite)
library(ahaz)
# set.seed(1234)
expit = function(d) return(exp(d)/(exp(d)+1))
# core = detectCores()
# cl = makeCluster(getOption("cl.cores", core))
# registerDoParallel(cl)

# Setting --------------------------------------------------------------------
out = vector(length = 1000)
for(i in 1:1000){
n = 10000
L = runif(n, 0, 1)
U = runif(n, 0, 1)
k = exp(0.075*L + 0.075*U)
out[i] = coef(lm(k~L))[2]
}
mean(out)
sd(out)


n = 800
r = 0.5
L = runif(n*6, 0, 1)
L = matrix(L, nrow = n)
U = runif(n, 0, 1)
Z_p = expit((L%*%rep(0.001, 6)))
Z = rbinom(n, 1, Z_p)
# C = rexp(n)/(0.01 + 0.05*L)
C = rexp(n)/(0.1 + L[, 1:6] %*% rep(0.05, 6))
max_t = 6
# W = rexp(n)/(0.01 + 0.01*L + 0.1*Z +0.01*U)
W = rexp(n)/(0.1 + Z + L[, 1:6] %*% rep(0.01, 6) + 0.01*U)

# F_time = function(t, w, z, l, u, uni) {
#   ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - (0.075*l + 0.075*u)*t)-uni, 
#         1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - (0.075*l
#                  + 0.075*u)*t) - uni)
# }
F_time = function(t, w, z, l, u, uni) {
  ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - (l%*%rep(0.25,6) + 0.25*u)*t)-uni, 
         1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - (l%*%rep(0.25,6) + 0.25*u)*t) - uni)
}

time = nleqslv(rep(0, n), F_time, w = W, z = Z, l = L, u = U, uni = runif(n))$x
time_c = ifelse(time <= C, time, C)
time_c = ifelse(time_c <= max_t, time_c, max_t)
event = ifelse(time<=C & time<=max_t, 1, 0)
sum(time_c > W)/800
cor(Z, W)

# Discrete treatment status ---------------------------------------------------

time_c = ceiling(time_c*10)/10
stime = sort(time_c[as.logical(event)])
stime = unique(stime)
k = length(stime)
D_status = treatment_status(n, k, stime, Z, W, max_t)

# Estimation -------------------------------------------------------------------

Cov = cbind(L, U)

system.time(s1 <- integral_est_cpp(c(0.2,0.2, 0), time = time_c, event = event, IV = Z,
                               Covariates = Cov, Covariates2 = Cov[, 1,drop = F], D_status = D_status, stime = stime))
print(s1[[1]])



# Data generation ---------------------------------------------------------

F_time = function(t, w, z, l, u, uni) {
  ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - exp(0.075*l + 0.075*u)*t)-uni, 
         1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - exp(0.075*l
                                                           + 0.075*u)*t) - uni)
}

n = 800
nrep = 1000
for (i in 1:nrep) {
  r = 0.5
  L = runif(n, 0, 1)
  L = matrix(L, nrow = n)
  U = runif(n, 0, 1)
  Z_p = expit(exp(r*L))
  Z = rbinom(n, 1, Z_p)
  C = rexp(n)/(0.01 + 0.05*L)
  max_t = 6
  W = rexp(n)/(0.01 + 0.01*L + 0.1*Z + 0.01*U)
  time = nleqslv(rep(0, n), F_time, w = W, z = Z, l = L, u = U, uni = runif(n))$x
  dat = cbind(L, U, Z, C, W, time)
  write_json(dat, paste0("DataGeneration/PS", n, "/", i, ".json"))
  if(i %% 100 == 0) cat("iter", i, "\n")
}

# Numerical experiment ---------------------------------------------------------

n = 800
nrep = 1000
parameters_dr = NULL
max_t = 6
for(i in 1:nrep){
  dat = read_json(paste0("DataGeneration/PS", n, "/", i, ".json"), 
                  simplifyVector = T)
  Cov = dat[, 1, drop = F]
  Z = dat[, 3]
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time_c <= max_t, time_c, max_t)
  time_c = ceiling(time_c*10)/10
  event = ifelse(time<=C & time<=max_t, 1, 0)
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

# write_json(parameters_dr, "SimulationResults/Outcome_DR_800.json")


dat = read_json("SimulationResults/Outcome_DR_1600.json", simplifyVector = T)
apply(dat, 1, mean)
apply(dat - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, mean)
apply(dat - matrix(c(0.1, 0.075), 
                             nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(999)))
apply(dat, 1, sd)


n = 800
nrep = 1000
parameters_dr = NULL
max_t = 6
for(i in 1:nrep){
  dat = read_json(paste0("DataGeneration/Outcome", n, "/", i, ".json"), 
                  simplifyVector = T)
  Cov = dat[, 1, drop = F]
  Z = dat[, 3]
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  # time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time <= W, time, W)
  # time_c = ifelse(time_c <= max_t, time_c, max_t)
  
  # time_c = ceiling(time_c*10)/10
  event = ifelse(time <= W, 1, 0)
  surv = Surv(time = time_c + runif(n, 0, 0.001), event = event, type = "right")
  out = coef(ahaz(surv, cbind(Z, Cov)))
  parameters_dr = cbind(parameters_dr, out)
  cat("[[", "rep ", i, "\t", out, "]]\n",
      sep = " ")
}
write_json(parameters_dr, "SimulationResults/Outcome_DR_1600_ahaz_censoring.json")
dat = read_json("SimulationResults/Outcome_DR_800_ahaz_censoring.json", simplifyVector = T)
apply(dat, 1, mean)
apply(dat - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, mean)
apply(dat - matrix(c(0.1, 0.075), 
                   nrow = 2, ncol = nrep), 1, function(d)sqrt(sum(d^2)/(999)))
apply(dat, 1, sd)



# ITT & tvmethod ----------------------------------------------------------

n = 1600
nrep = 1000
out1 = NULL
out2 = NULL
max_t = 6
for(i in 1:nrep){
  dat = read_json(paste0("DataGeneration/PS", n, "/", i, ".json"), 
                  simplifyVector = T)
  Cov = dat[, 1, drop = F]
  Z = dat[, 3]
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time_c <= max_t, time_c, max_t)
  event = T_D < C & T_D <= max_t
  event_w = T_D_c > W
  tvdat = data.frame(id = 1:n, treatment = Z, X1 = Cov[, 1], event = event, 
                     start_time = 0, end_time = T_D_c)
  swdat = tvdat[event_w, ]
  tvdat$end_time[event_w] = W[event_w]
  tvdat$event[event_w] = F
  swdat$start_time = W[event_w]
  swdat$treatment = 1 - swdat$treatment
  adat = rbind(tvdat, swdat)
  adat = adat[order(adat$id, adat$start_time), ]
  if(any(adat$end_time == 0)){
    nn = length(adat$end_time[adat$end_time == 0])
    adat$end_time[adat$end_time == 0] = runif(nn, 0, 0.01)
  }
  
  mod  = aalen(Surv(start_time, end_time, event) ~ const(treatment) + const(X1), data = adat, max.time = 6, id = adat$id)
  surv = Surv(T_D_c + runif(n, 0, 0.01), event = event, type = "right")
  mod2 = ahaz(surv, cbind(Z, Cov))
  out1 = cbind(out1, coef(mod)[, 1])
  out2 = cbind(out2, coef(mod2)[1])
  cat("[[", "rep ", i, "\t", coef(mod)[1, 1], "\t" ,coef(mod2)[1], "]]\n",
      sep = " ")
}
out = list(tv = out1, itt = out2)
write_json(out, "SimulationResults/PS_DR_1600_itt_tv.json")


# Exploring ---------------------------------------------------------------



n = 800
nrep = 1000
max_t = 6
out1 = vector(length = nrep)
out2 = vector(length = nrep)
for (i in 1:nrep) {
  dat = read_json(paste0("DataGeneration/Outcome", n, "/", i, ".json"), 
                  simplifyVector = T)
  time = dat[, 6]
  C = dat[, 4]
  W = dat[, 5]
  Z = dat[, 3]
  time_c = ifelse(time <= C, time, C)
  time_c = ifelse(time_c <= max_t, time_c, max_t)
  event = ifelse(time<=C & time<=max_t, 1, 0)
  out1[i] = sum(ifelse(Z == 1, time_c > W, 0))
  out2[i] = sum(ifelse(Z == 0, time_c > W, 0))
}
mean(out1)
mean(out2)
