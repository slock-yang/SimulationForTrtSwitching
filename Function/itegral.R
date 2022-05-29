# this file needs to be altered in function tmp_citeg
# the form of this function can be similar to iteg_expbD_individual
expit = function(d) return(exp(d)/(exp(d)+1))
itegral = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime, tol = 1e-3, xtol = 1e-5, iter = 20)
{
  # Covariates must be matrix
  # setting
  # ------------------------------------------------------------------------------
  
  n = length(time)
  k = length(stime)
  betaD = init_parameters[1]
  beta = init_parameters[-1]
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))

  # calculus
  # ------------------------------------------------------------------------------

  int_D = matrix(0, nrow = n, ncol = k)

  for(j in 1:k){
      if(j == 1){
      int_D[, j] = drop(IV  * stime[j])
    }else{
      int_D[, j] = int_D[, j-1] + 
        drop(D_status[, j-1]*(stime[j] - stime[j-1]))
    }
  }
  
  idx_event = sapply(time, function(d) return(which(time == stime)))
  switching = apply(D_status, 1, function(d) return(any(diff(d) != 0)))
  idx_switching = rep(length(stime), nrow(D_status))
  idx_switching[switching] = apply(D_status[switching, ], 1, 
                        function(d) return(which(diff(d) != 0) + 1))
  # idx_switching = apply(D_status, 1, function(d) return(which(diff(d) != 0)))
  switching = idx_switching <= idx_event

  int_expbd = iteg_expbD(D_status, betaD, IV, stime, switching, idx_switching)
  
  p1 = exp(beta_D * int_D[, idx_event]) * event
  p2 = int_expbd
  
}


# iteg_expbD = function(D_status, betaD, IV, time, stime, switching, idx_switching) {
#     if(betaD == 0){
#         int_expbd = time
#     } else {
#         int_expbd = ifelse(switching, 
#                         ifelse(IV == 0, stime[idx_switching] + 
#                                 (exp(betaD * (time - stime[idx_switching])) - 1)/betaD, 
#                             (exp(betaD * stime[idx_switching]) - 1)/betaD + 
#                                 (time - stime[idx_switching]) * 
#                                 exp(betaD * stime[idx_switching])),
#                         ifelse(IV == 0, time, (exp(betaD * time) - 1)/betaD))
#     }
#     return(int_expbd)
# }


iteg = function(D_status, betaD, IV, time, stime, switching, idx_switching, 
                    int_D) {
    idxD_status = cbind(1:nrow(D_status), D_status)
    int_expbD = apply(idxD_status, 1, iteg_expbD_individual, betaD = betaD, beta = beta,
                        IV = IV, time = time, stime = stime, int_D, Covariates)
    int_cexpbD_Lam_nm = apply(idxD_status, 1, iteg_cexpbD_Lam_individual, D_status = D_status,
                        betaD = betaD, beta = beta, IV = IV, time = time, stime = stime,
                        switching = switching, idx_switching = idx_switching, Covariates = Covariates)


    return(-(int_expbD + int_cexpbD_Lam_nm))
}


iteg_expbD_individual = function(idxD_status, betaD, beta, IV, time, stime, 
                                    int_D, Covariates){
    k = length(stime)
    idx_i = idxD_status[1]
    Yt = time[idx_i] >= stime
    D_status_individual = idxD_status[-1]
    time = time[idx_i]
    int_D = int_D[idx_i, ]
    IV = IV[idx_i]
    Covariates = Covariates[idx_i, , drop = FALSE]
    tmp_D = c(IV, D_status_individual[-k])
    diff_stime = diff(c(0, stime))
    if(betaD != 0){
        int_expbD = ifelse(tmp_D == 0, diff_stime*c(1, exp(betaD*int_D[-k]),
                                        diff(c(1, exp(betaD*int_D)))/betaD))
    } else {
        int_expbD = diff_stime
    }

    return(-sum(Yt*int_expbD*D_status_individual*betaD + 
        Yt*int_expbD*drop(Covariates%*%beta)))
}


iteg_cexpbD_Lam_individual = function(idxD_status, D_status, betaD, beta, IV, time, stime,
                                switching, idx_switching, event, Covariates) {
    idx_i = idxD_status[1]
    D_status_individual = idxD_status[-1]
    cD_status = t(D_status_individual + t(D_status))
    cIV = IV[idx_i] + IV
    idxcD_status = cbind(1:nrow(D_status), cD_status)

    int_cexpbD_Lam_nm_individual = apply(idxcD_status, 1, tmp_citeg, betaD = betaD,
                                time = time, stime = stime, cIV = cIV, idx_i = idx_i)
    
    return(sum(int_cexpbD_Lam_nm_individual[3, ] * event -
            int_cexpbD_Lam_nm_individual[2, ] * beta_D - 
            int_cexpbD_Lam_nm_individual[1, ] * (Covariates %*% beta)))
}


tmp_citeg = function(idxcD_status, D_status, betaD, time, stime, cIV, idx_i) {
    idx_j = idxcD_status[1]
    cD_status = idxcD_status[-1]
    cg_p = which(diff(cD_status) != 0) + 1
    time_min = min(time[idx_i], time[idx_j])
    cg_event = which(time_min == stime)
    cg_p = c(cg_p[cg_p < cg_event], cg_event)
    pri = vector(length = length(cg_p))
    D = D_status[cg_p - 1]
    int_cD = cD_status

    for(j in 1:length(stime)){
        if(j == 1) {
            int_cD[j] = drop(cIV[idx_j] * stime[j])
        } else {
            int_cD[j] = int_cD[j-1] + drop(cD_status[j-1] * 
                    (stime[j] - stime[j-1]))
        }
    }

    for(i in 1:length(cg_p)){
        if(betaD == 0){
            if(i == 1){
                pri[i] = stime[cg_p[i]]
                priD[i] = pri[i] * 
            } else {
                pri[i] = stime[cg_p[i]] - stime[cg_p[i-1]]
            }
        } else {
            if(cD_status[cg_p[i]-1] == 0){
                if(i == 1){
                    pri[i] = stime[cg_p[i]]
                } else {
                    pri[i] = exp(betaD * int_cD[cg_p[i]]) * 
                                (stime[cg_p[i]] - stime[cg_p[i-1]])
                }
            } else {
                if(i == 1){
                    pri[i] = (exp(betaD*int_cD[cg_p[i]])-1)/
                                (cD_status[cg_p[i]-1]*betaD)
                } else {
                    pri[i] = (exp(betaD*int_cD[cg_p[i]]) - 
                                    exp(betaD*int_cD[cg_p[i-1]]))/
                                (cD_status[cg_p[i]-1]*betaD)
                }
            }
        }
        
    }

    priD = pri * D

    return(c(sum(pri), sum(priD), exp(betaD * int_cD[cg_event])))
}