source("Fun.R")
library(ivsacim)
library(nleqslv)
library(ahaz)
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

s_iv = ivsacim(time, event, Z, treatment_init = Z)

s = SCSM(time, event, Z, L, D_status, stime)

f1 = stepfun(c(stime), c(0, s$BD))
f2 = stepfun(c(stime), cbind(0, s$beta))
cat(f1(1), f1(2), f1(3), "\n")
cat(f2(1), f2(2), f2(3), "\n")
print(s_iv$beta)
