source("Fun.R")
library(ivsacim)
library(nleqslv)
# set.seed(1234)
expit = function(d) return(exp(d)/(exp(d)+1))

# ================================================
# ================ Basic setting =================
# ================================================
n = 500
r = 0.5
L = runif(n, 0, 1)
L = matrix(L, nrow = n)
U_1 = runif(n, 0, 1)
U_2 = runif(n, 0, 1)
Z = rbinom(n, 1, 0.5)

# ========================================================
# ==== treatment switching time and safety outcome =======
# ========================================================

# F_time = function(t, w, z, l, u, uni) {
#   return(1 - exp(-0.25*t - 0.1*z*t - 0.075*l*t- 0.075*u*t)-uni)
# }
# uni = runif(n)
# time = rep(0, n)
# for(i in 1:n){
#   temp = nleqslv(0, F_time, w = W[i], z = Z[i], l = L[i], u = U_2[i], 
#                  uni = uni[i])
#   time[i] = temp$x
# }
time = rexp(n)/(0.25 + 0.1*Z + 0.075 * L)
event = rep(1, n)


# ====================================================
# =========== discrete treatment status ==============
# ====================================================

stime = c(time[as.logical(event)])
stime = sort(stime)
k = length(stime)
max_t = max(time)
W = rep(0, n)
D_status = treatment_status(n, k, stime, Z, W, max_t)


# ====================================================
# ================   Estimation   ====================
# ====================================================


# s = nleqslv(rep(0, 3), Lin_Ying, time = time, stime = stime,
#             event_new = event, Z = Z, L = L, U = U_2,  D_status = D_status)
s1 = nleqslv(c(0,0), ConstantF, time = time, event = event, IV = Z, 
             Covariates = L, D_status = D_status, stime = stime)
print(s1$x)
s_iv = ivsacim(time, event, Z, treatment_init = Z)
print(s_iv$beta_D)


# ============================= rep = 100 ===============================

nrep = 100
n = 1000
Constant_beta = rep(0, nrep)
Constant_alpha = rep(0, nrep)
IVest = rep(0, nrep)
for(i in 1:nrep){
  L = runif(n, 0, 1)
  L = matrix(L, nrow = n)
  U_1 = runif(n, 0, 1)
  U_2 = runif(n, 0, 1)
  Z = rbinom(n, 1, 0.5)
  time = rexp(n)/(0.25 + 0.1*Z + 0.075 * L)
  event = rep(1, n)
  stime = c(time[as.logical(event)])
  stime = sort(stime)
  k = length(stime)
  max_t = max(time)
  W = rep(0, n)
  D_status = treatment_status(n, k, stime, Z, W, max_t)
  s1 = nleqslv(c(0,0), ConstantF, time = time, event = event, IV = Z, 
               Covariates = L, D_status = D_status, stime = stime)
  Constant_beta[i] = s1$x[1]
  Constant_alpha[i] = s1$x[2]
  s_iv = ivsacim(time, event, Z, treatment_init = Z)
  IVest[i] =  s_iv$beta_D
  cat("[[", "rep ", i, "  Constant: ", s1$x, "  IVest: ",s_iv$beta_D, "]]\n",
      sep = "")
}
print(mean(Constant_beta));print(mean(Constant_alpha))
print(sd(Constant_beta));print(sd(Constant_alpha))
print(mean(IVest))
print(sd(IVest))
