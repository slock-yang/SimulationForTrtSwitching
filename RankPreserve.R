# source("Function/Fun.R")
# source("Function/Fun_parallel.R")
source("Function/integral.R")
source("RankPreserve_fun.R")

library(ivsacim)
library(nleqslv)
library(jsonlite)
library(ahaz)
library(timereg)
expit = function(d) return(exp(d)/(exp(d)+1))




# unmeasured confounding --------------------------------------------------

expit = function(d) return(exp(d)/(exp(d)+1))

n =  800
p = 7
max_t = 6
Cov = matrix(runif(n*p), ncol = p, nrow = n)
beta = rep(0.25, p)
Z_p = expit((Cov %*% rep(0.001, p)))
Z = rbinom(n, 1, Z_p)
table(Z)

T0 = rexp(n)/(0.25 + Cov %*% beta)
T1 = rexp(n)/(0.35 + Cov %*% beta)
hist(T0)
hist(T1)

alpha = rep(0.25, p)
W = ifelse(as.logical(Z), rexp(n)/(exp(-T1)*(0.1 - Z + Cov %*% alpha)),
           rexp(n)/(exp(-T0)*(0.1 + (1 - Z) + Cov %*% alpha)))
# W = rexp(n)/(exp(-T0)*(0.01 + Z + Cov %*% alpha))
all(W > 0)
W[W < 0] = Inf

F_t = function(t, t0, w, Cov, Z, beta){
  return(ifelse(t < w, 
                                 0.25* t + Cov %*% beta * (t-t0) - 0.25 * t0 + 0.1 * Z * t,
                                 0.25* t + Cov %*% beta * (t-t0) - 0.25 * t0 + 0.1 * Z * w + 0.1 * (1-Z) * (t-w)))
}
T_D = nleqslv(rep(0, n), F_t, t0 = T0, w = W, Cov = Cov, Z = Z, beta = beta)$x
C = rexp(n)/(0.1 + Cov[, 1:6] %*% rep(0.05, 6))

table(Z, T_D > C)
T_D_c = ifelse(T_D < C, T_D, C)
T_D_c = ifelse(T_D_c <= max_t, T_D_c, max_t)
table(Z, T_D_c > W)

# time-varying covariates
event = T_D < C & T_D <= max_t
event_w = T_D_c > W
tvdat = data.frame(id = 1:n, treatment = Z, Cov[, 1:6], event = event, 
                   start_time = 0, end_time = T_D_c)
swdat = tvdat[event_w, ]
tvdat$end_time[event_w] = W[event_w]
tvdat$event[event_w] = F
swdat$start_time = W[event_w]
swdat$treatment = 1 - swdat$treatment
dat = rbind(tvdat, swdat)
dat = dat[order(dat$id, dat$start_time), ]

out  = aalen(Surv(start_time, end_time, event) ~ const(treatment) + const(X1) + const(X2) +
        const(X3) + const(X4) + const(X5) + const(X6), data = dat, max.time = 6, id = dat$id)

event = T_D_c < W
T_D_c = ifelse(T_D < W, T_D, W)

surv = Surv(T_D_c, event = event, type = "right")
out = ahaz(surv, cbind(Z, Cov[, 1:6]))

T_D_l = ceiling(T_D*100)/100
stime = sort(T_D_l)
stime = unique(stime)
k = length(stime)
D_status = matrix(nrow = n, ncol = k)
for (i in 1:n) {
  if(T_D[i] > W[i]){
    D_status[i, which(stime <= W[i])] = Z[i]
    D_status[i, which(stime > W[i])] = 1 - Z[i]
  } else {
    D_status[i, ] = Z[i]
  }
}
integral_est_cpp(c(0, rep(0, p)), time = T_D_l, event = rep(1, n),
                 IV = Z, Covariates = Cov, Covariates2 = Cov, 
                 D_status = D_status, stime = stime)
coef(out)



# dr_z ----------------------------------------------------------------------

expit = function(d) return(exp(d)/(exp(d)+1))

n =  800
p = 7
max_t = 6
Cov = matrix(runif(n*p), ncol = p, nrow = n)
beta = rep(0.25, p)

tran = function(d) return(ifelse(d > 0, exp(d), -exp(-d)))

Z_p = expit(tran(Cov %*% rep(c(1, -1), length.out = p)))
Z = rbinom(n, 1, Z_p)
table(Z)

T0 = rexp(n)/(0.25 + Cov %*% beta)
T1 = rexp(n)/(0.35 + Cov %*% beta)
hist(T0)
hist(T1)

alpha = rep(0.25, p)
W = ifelse(as.logical(Z), rexp(n)/(exp(-T1)*(0.1 - Z + Cov %*% alpha)),
           rexp(n)/(exp(-T0)*(0.1 + (1 - Z) + Cov %*% alpha)))
# W = rexp(n)/(exp(-T0)*(0.01 + Z + Cov %*% alpha))
all(W > 0)
W[W < 0] = Inf

F_t = function(t, t0, w, Cov, Z, beta){
  return(ifelse(t < w, 
                0.25* t + Cov %*% beta * (t-t0) - 0.25 * t0 + 0.1 * Z * t,
                0.25* t + Cov %*% beta * (t-t0) - 0.25 * t0 + 0.1 * Z * w + 0.1 * (1-Z) * (t-w)))
}
T_D = nleqslv(rep(0, n), F_t, t0 = T0, w = W, Cov = Cov, Z = Z, beta = beta)$x
C = rexp(n)/(0.1 + Cov[, 1:6] %*% rep(0.05, 6))

table(Z, T_D > C)
T_D_c = ifelse(T_D < C, T_D, C)
T_D_c = ifelse(T_D_c <= max_t, T_D_c, max_t)
table(Z, T_D_c > W)

# time-varying covariates
event = T_D < C & T_D <= max_t
event_w = T_D_c > W
tvdat = data.frame(id = 1:n, treatment = Z, Cov[, 1:6], event = event, 
                   start_time = 0, end_time = T_D_c)
swdat = tvdat[event_w, ]
tvdat$end_time[event_w] = W[event_w]
tvdat$event[event_w] = F
swdat$start_time = W[event_w]
swdat$treatment = 1 - swdat$treatment
dat = rbind(tvdat, swdat)
dat = dat[order(dat$id, dat$start_time), ]

out  = aalen(Surv(start_time, end_time, event) ~ const(treatment) + const(X1) + const(X2) +
               const(X3) + const(X4) + const(X5) + const(X6), data = dat, max.time = 6, id = dat$id)

event = T_D_c < W
T_D_c = ifelse(T_D < W, T_D, W)

surv = Surv(T_D_c, event = event, type = "right")
out = ahaz(surv, cbind(Z, Cov[, 1:6]))

T_D_l = ceiling(T_D*100)/100
stime = sort(T_D_l)
stime = unique(stime)
k = length(stime)
D_status = matrix(nrow = n, ncol = k)
for (i in 1:n) {
  if(T_D[i] > W[i]){
    D_status[i, which(stime <= W[i])] = Z[i]
    D_status[i, which(stime > W[i])] = 1 - Z[i]
  } else {
    D_status[i, ] = Z[i]
  }
}
integral_est_cpp(c(0.1, rep(0.25, p)), time = T_D_l, event = rep(1, n),
                 IV = Z, Covariates = Cov, Covariates2 = Cov, 
                 D_status = D_status, stime = stime)


coef(out)



# dr_outcome --------------------------------------------------------------

expit = function(d) return(exp(d)/(exp(d)+1))

n =  800
p = 7
max_t = 6
Cov = matrix(runif(n*p), ncol = p, nrow = n)
beta = rep(0.25, p)

tran = function(d) return(ifelse(d > 0, exp(d), -exp(-d)))

Z_p = expit((Cov %*% rep(c(1, -1), length.out = p)))
Z = rbinom(n, 1, Z_p)
table(Z)

T0 = rexp(n)/(0.25 + (exp(Cov %*% beta) - 1))
T1 = rexp(n)/(0.35 + (exp(Cov %*% beta) - 1))
hist(T0)
hist(T1)

alpha = rep(0.25, p)
W = ifelse(as.logical(Z), rexp(n)/(exp(-T1)*(0.1 - Z + exp(Cov %*% alpha) - 1)),
           rexp(n)/(exp(-T0)*(0.1 + (1 - Z) + exp(Cov %*% alpha) - 1)))
# W = rexp(n)/(exp(-T0)*(0.01 + Z + Cov %*% alpha))
all(W > 0)
W[W < 0] = Inf

F_t = function(t, t0, w, Cov, Z, beta){
  return(ifelse(t < w, 
                0.25* t + (exp(Cov %*% beta) - 1) * (t-t0) - 0.25 * t0 + 0.1 * Z * t,
                0.25* t + (exp(Cov %*% beta) - 1) * (t-t0) - 0.25 * t0 + 0.1 * Z * w + 0.1 * (1-Z) * (t-w)))
}
T_D = nleqslv(rep(0, n), F_t, t0 = T0, w = W, Cov = Cov, Z = Z, beta = beta)$x
C = rexp(n)/(0.1 + Cov[, 1:6] %*% rep(0.05, 6))

table(Z, T_D > C)
T_D_c = ifelse(T_D < C, T_D, C)
T_D_c = ifelse(T_D_c <= max_t, T_D_c, max_t)
table(Z, T_D_c > W)

# time-varying covariates
event = T_D < C & T_D <= max_t
event_w = T_D_c > W
tvdat = data.frame(id = 1:n, treatment = Z, Cov[, 1:6], event = event, 
                   start_time = 0, end_time = T_D_c)
swdat = tvdat[event_w, ]
tvdat$end_time[event_w] = W[event_w]
tvdat$event[event_w] = F
swdat$start_time = W[event_w]
swdat$treatment = 1 - swdat$treatment
dat = rbind(tvdat, swdat)
dat = dat[order(dat$id, dat$start_time), ]

out  = aalen(Surv(start_time, end_time, event) ~ const(treatment) + const(X1) + const(X2) +
               const(X3) + const(X4) + const(X5) + const(X6), data = dat, max.time = 6, id = dat$id)

event = T_D_c < W
T_D_c = ifelse(T_D < W, T_D, W)

surv = Surv(T_D_c, event = event, type = "right")
out = ahaz(surv, cbind(Z, Cov[, 1:6]))

T_D_l = ceiling(T_D*100)/100
stime = sort(T_D_l)
stime = unique(stime)
k = length(stime)
D_status = matrix(nrow = n, ncol = k)
for (i in 1:n) {
  if(T_D[i] > W[i]){
    D_status[i, which(stime <= W[i])] = Z[i]
    D_status[i, which(stime > W[i])] = 1 - Z[i]
  } else {
    D_status[i, ] = Z[i]
  }
}
integral_est_cpp(c(0.1, rep(0.25, p)), time = T_D_l, event = rep(1, n),
                 IV = Z, Covariates = Cov, Covariates2 = Cov, 
                 D_status = D_status, stime = stime)


coef(out)



# Simulation --------------------------------------------------------------

for (nn in c(800, 1600, 3200)) {
  for(ii in c("", "outcome", "IV")){
    for (zz in c("s1", "s2")) {
      
        for (nu in c(-1, 0, 1)) {
          dir.create(paste0("DataGeneration/", paste0(zz, ii, nu, nn)))
          if(ii == ""){
            trans = function(d) d
            transZ = function(d) d
          } else {
            if(ii == "outcome"){
              trans = function(d) return(exp(d) - 1)
              transZ = function(d) d
            } else {
              trans = function(d) d
              transZ = function(d) return(ifelse(d < 0, -exp(-d), exp(d)))
            }
          }
          newdir.= paste0(zz, ii, nu)
          if(zz == "s1"){
            sim_data(nrep = 1000, n = nn, p =10, nu = nu, trans = trans,
                   transZ = transZ,
                   newdir = newdir.,
                   root = root)
          } else {
            sim2_data(nrep = 1000, n = nn, p =10, nu = nu, trans = trans,
                      transZ = transZ,
                      newdir =newdir.,
                      root = root)
          }
          
          cat(paste0(zz, ii, nu, nn), "\n")
        }
        
      
    }
  }
}


for (nn in c(800, 1600)) {
  for (ii in c("", "outcome", "IV")) {
    for (zz in c("s1", "s2")){
      for(nu in c(-1, 0, 1)) {
        out = sim2(nrep = 1000, n = nn, newdir = paste0(zz, ii, nu),
                   root = root)
        write_json(out, paste0("simResults/", paste0(nn, zz, ii, nu), ",json"))
        cat(paste0(nn, zz, ii, nu), "\n")
      }
    }
  }
}

sim_data(nrep = 1000, n = 800, p =10, nu = 0, trans = trans,
         transZ = transZ,
         newdir = "s1-1",
         root = root)
out = sim2(nrep = 1000, n = 800, newdir = "s1-1", root = root)


