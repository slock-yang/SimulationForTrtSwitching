source("Function/Fun.R")
source("Function/Fun_parallel.R")
source("Function/Fun_parallel_est.R")
source("Function/integral.R")
library(ivsacim)
library(nleqslv)
library(ahaz)
# set.seed(1234)
expit = function(d) return(exp(d)/(exp(d)+1))
core = detectCores()/2
cl = makeCluster(getOption("cl.cores", core))
registerDoParallel(cl)

# ================================================
# ================ Basic setting =================
# ================================================
n = 1000
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

IV = Z
Covariates = L
init_parameters = c(0.1, 0.075)

system.time(s <- nleqslv(c(0.1,0.075), integral_cpp, time = time, 
                  event = event, IV = Z,
            Covariates = L, D_status = D_status, stime = stime, jacobian = TRUE, method = "Newton", global = "none"))
print(s$x)

system.time(ss <- integral_est_debug_cpp(c(0,0), time = time, 
                  event = event, IV = Z,
            Covariates = L, D_status = D_status, stime = stime, max_iter = 20, tol = 1e-2))
print(ss)

system.time(sss <- integral_est_cpp(c(0,0), time = time, 
                  event = event, IV = Z,
            Covariates = L, D_status = D_status, stime = stime, max_iter = 10, tol = 1e-3))
print(sss)

system.time(ssss <- ConstantF_parallel_Newtonest(c(0.1,0.075), time = time, 
                  event = event, IV = Z,
            Covariates = L, D_status = D_status, stime = stime))
print(ssss)
print(integral_cpp(sss$x, time = time, 
                  event = event, IV = Z,
            Covariates = L, D_status = D_status, stime = stime))



nrep = 100
n = 500
Constant_beta = rep(0, nrep)
Constant_alpha = rep(0, nrep)
# ahaz_beta = rep(0, nrep)
# ahaz_alpha = rep(0, nrep)
# IVest = rep(0, nrep)
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
  s1 = integral_est_cpp(c(0,0), time = time, event = event, IV = Z, 
               Covariates = L, D_status = D_status, stime = stime)
  Constant_beta[i] = s1$x[1]
  Constant_alpha[i] = s1$x[2]
#   ahaz_beta[i] = summary(s2)$coefficients[1, 1]
#   ahaz_alpha[i] = summary(s2)$coefficients[2, 1]
#   s_iv = ivsacim(time, event, Z, treatment_init = Z)
#   IVest[i] =  s_iv$beta_D
  cat("[[", "rep ", i, "  Constant: ", s1$x, "]]\n")
}
print(mean(Constant_beta));print(mean(Constant_alpha))
# print(mean(ahaz_beta)); print(mean(ahaz_alpha))
print(sd(Constant_beta));print(sd(Constant_alpha))
# print(sd(ahaz_beta)); print(sd(ahaz_alpha))
# print(mean(IVest))
# print(sd(IVest))