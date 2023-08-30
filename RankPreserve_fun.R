library(jsonlite)
root = "/Users/liuyang/Documents/Boston/MS SecondLine/Code & Analysis/DRcomparasion/ExampleforTrtSwitching/SimulationForTrtSwitching/DataGeneration/"

trans = function(d) return(exp(d) - 1)
transZ = function(d) return(ifelse(d < 0, -exp(-d), exp(d)))


sim_data = function(nrep = 1000, n = 800, p = 10,
                    beta = rep(0.25, p),
                    alpha = rep(0.25, p),
                    gamma = rep(c(1, -1), length.out = p-1),
                    nu = 0,
                    trans = function(d) d,
                    transZ = function(d) d,
                    json_save = TRUE,
                    newdir = "",
                    root = "")
{
  F_time = function(t, w, z, Cov, uni) {
    ifelse(t < w, 1 - exp(-0.25*t - 0.1*z*t - trans(Cov %*% beta)*t)-uni, 
           1 - exp(-0.25*t - 0.1*z*w - 0.1*(1-z)*(t-w) - trans(Cov %*% beta)*t) - uni)
  }
  
  
  Z_proportion = NULL
  switching_rate_overall = NULL
  switching_rate_from_0 = NULL
  switching_rate_from_1 = NULL
  censoring_rate = NULL
  for (kk in 1:nrep) {
    Cov = matrix(runif(n*p), ncol = p, nrow = n)
    Z_p = expit(transZ(Cov[, 1:(p-1)] %*% gamma))
    Z = rbinom(n, 1, Z_p)
    C = rexp(n)/(0.1 + Cov[, 1:(p-1)] %*% rep(0.05, p-1))
    W = ifelse(as.logical(Z), rexp(n)/(0.5*(nu + Z + Cov %*% alpha)),
               rexp(n)/(0.5*(nu - (1 - Z) + Cov %*% alpha)))
    W[W <= 0] = Inf
    T_D = nleqslv(rep(0, n), F_time, w = W, z = Z, Cov = Cov, uni = runif(n))$x
    dat = cbind(Cov, Z, T_D, W, C)
    if(json_save) write_json(dat, paste0(root, newdir, n, "/", kk, ".json"))
    if(kk %% 100 == 0) cat("iter", kk, "\n")
    
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= 6, T_D_c, 6)
    event = T_D <= 6 & T_D <= C
    Z_proportion = c(Z_proportion, mean(Z))
    switching_rate_overall = c(switching_rate_overall, mean(T_D_c > W))
    switching_rate_from_0 = c(switching_rate_from_0, mean((T_D_c>W)[Z == 0]))
    switching_rate_from_1 = c(switching_rate_from_1, mean((T_D_c>W)[Z == 1]))
    censoring_rate = c(censoring_rate, mean(1-event))
  }
  list(Z_proportion = mean(Z_proportion),
       switching_rate_overall = mean(switching_rate_overall),
       switching_rate_from_0 = mean(switching_rate_from_0),
       switching_rate_from_1 = mean(switching_rate_from_1),
       censoring_rate = mean(censoring_rate))
}






sim2_data = function(nrep = 1000, n = 800, p = 10,
                     beta = rep(0.25, p),
                     alpha = rep(0.25, p),
                     gamma = rep(c(1, -1), length.out = p-1),
                     nu = 0,
                     trans = function(d) d,
                     transZ = function(d) d,
                     json_save = TRUE,
                     newdir = "",
                     root = "")
{
  Rank_preserving = function(t, t0, w, Cov, Z, beta){
    return(ifelse(t < w, 
                  0.25* t + trans(Cov %*% beta) * (t-t0) - 0.25 * t0 + 0.1 * Z * t,
                  0.25* t + trans(Cov %*% beta) * (t-t0) - 0.25 * t0 + 0.1 * Z * w + 0.1 * (1-Z) * (t-w)))
  }
  
  
  F_t0 = function(nu, Cov, beta) {
    return(rexp(nu)/(0.25 + Cov %*% beta))
  }
  
  
  Z_proportion = NULL
  switching_rate_overall = NULL
  switching_rate_from_0 = NULL
  switching_rate_from_1 = NULL
  censoring_rate = NULL
  for(kk in 1:nrep){
    Cov = matrix(runif(n*p), ncol = p, nrow = n)
    Z_p = expit(transZ(Cov[, 1:(p-1)] %*% gamma))
    Z = rbinom(n, 1, Z_p)
    T0 = F_t0(n, Cov, beta)
    T1 = nleqslv::nleqslv(rep(0, n), Rank_preserving, 
                          t0 = T0, w = Inf, Cov = Cov, Z = 1,
                          beta = beta)$x
    
    W = ifelse(as.logical(Z), rexp(n)/(exp(-T1)*(nu + Z + Cov %*% alpha)),
               rexp(n)/(exp(-T0)*(nu - (1 - Z) + Cov %*% alpha)))
    W[W <= 0] = Inf
    C = rexp(n)/(0.1 + Cov[, 1:(p-1)] %*% rep(0.05, p-1))
    T_D = nleqslv(rep(0, n), Rank_preserving, 
                  t0 = T0, w = W, Cov = Cov, Z = Z, beta = beta)$x
    dat = cbind(Cov, Z, T_D, W, C)
    if(json_save) write_json(dat, paste0(root, newdir, n, "/", kk, ".json"))
    if(kk %% 100 == 0) cat("iter", kk, "\n")
    
    T_D_c = ifelse(T_D <= C, T_D, C)
    T_D_c = ifelse(T_D_c <= 6, T_D_c, 6)
    event = T_D <= 6 & T_D <= C
    
    Z_proportion = c(Z_proportion, mean(Z))
    switching_rate_overall = c(switching_rate_overall, mean(T_D_c > W))
    switching_rate_from_0 = c(switching_rate_from_0, mean((T_D_c>W)[Z == 0]))
    switching_rate_from_1 = c(switching_rate_from_1, mean((T_D_c>W)[Z == 1]))
    censoring_rate = c(censoring_rate, mean(1-event))
  } 
  list(Z_proportion = mean(Z_proportion),
       switching_rate_overall = mean(switching_rate_overall),
       switching_rate_from_0 = mean(switching_rate_from_0),
       switching_rate_from_1 = mean(switching_rate_from_1),
       censoring_rate = mean(censoring_rate))
}


sim2 = function(nrep = 1000, n = 800, p = 10,
                max_t = 6,
                newdir = "",
                root = ""){
  out1 = NULL
  out2 = NULL
  out2_var = NULL
  out3 = NULL
  out4 = NULL
  out5 = NULL
  Con = vector(length = nrep)
  for (kk in 1:nrep) {
    dat = read_json(paste0(root, newdir, n, "/", kk, ".json"), 
                    simplifyVector = T)
    Cov = dat[, 1:p]
    Z = dat[, p+1]
    T_D = dat[, p+2]
    W = dat[, p+3]
    C = dat[, p+4]
    T_D_c = ifelse(T_D < C, T_D, C)
    T_D_c = ifelse(T_D_c <= max_t, T_D_c, max_t)
    
    event = T_D < C & T_D <= max_t
    event_w = T_D_c > W
    tvdat = data.frame(id = 1:n, treatment = Z, Cov[, 1:(p-1)], event = event, 
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
    
    mod5  = aalen(Surv(start_time, end_time, event) ~ const(treatment) + const(X1) + const(X2) +
                   const(X3) + const(X4) + const(X5) + const(X6) + const(X7) + const(X8) + const(X9), data = adat, max.time = 6, id = adat$id)
    
    # W = 0.9 * T0
    
    event = T_D < C & T_D < max_t 
    event_w = event & T_D < W
    T_D_c_w = ifelse(T_D_c < W, T_D_c, W)
    surv = Surv(T_D_c_w + runif(n, 0, 0.01), event = event_w, type = "right")
    mod1 = ahaz(surv, cbind(Z, Cov[, 1:(p-1)]))
    
    surv2 = Surv(T_D_c[event] + runif(sum(event), 0, 0.01), event = event[event], type = "right")
    mod3 = ahaz(surv2, cbind(Z[event], Cov[event, 1:(p-1)]))
    
    surv3 = Surv(T_D_c + runif(n, 0, 0.01), event = event, type = "right")
    mod4 = ahaz(surv3, cbind(Z, Cov[, 1:(p-1)]))
    
    T_D_c = ceiling(T_D_c * 100)/100
    stime = sort(T_D_c)
    stime = unique(stime)
    k = length(stime)
    D_status = matrix(nrow = n, ncol = k)
    for (i in 1:n) {
      if(T_D_c[i] > W[i]){
        D_status[i, which(stime <= W[i])] = Z[i]
        D_status[i, which(stime > W[i])] = 1 - Z[i]
      } else {
        D_status[i, ] = Z[i]
      }
    }
    mod2 = integral_est_cpp(c(0.1, rep(0.25, p-1)), time = T_D_c, event = event,
                            IV = Z, Covariates = Cov[, 1:(p-1)], Covariates2 = Cov, 
                            D_status = D_status, stime = stime)
    out1 = cbind(out1, coef(mod1))
    out3 = cbind(out3, coef(mod3))
    out4 = cbind(out4, coef(mod4))
    out5 = cbind(out5, coef(mod5)[, 1])
    out2 = cbind(out2, mod2$x)
    out2_var = c(out2_var, mod2$var)
    Con[kk] = mod2$Convergence
    cat("[[\t rep", kk, "\t", coef(mod1)[1], "\t", coef(mod3)[1], "\t", coef(mod4)[1], "\t", coef(mod5)[1, 1], "\t",mod2$x[1],"\t", mod2$var, "\t", mod2$Convergence,
        "\t", "]]\n")
  }
  return(list(ahaz = out1,
              delete = out3,
              itt = out4,
              timeVar = out5,
              dr = out2,
              dr_var = out2_var,
              Con = Con))
}
