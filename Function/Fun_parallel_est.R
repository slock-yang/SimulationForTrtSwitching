library(parallel)
library(doParallel)


expit = function(d) return(exp(d)/(exp(d)+1))
ConstantF_parallel_est = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime, tol = 1e-3, xtol = 1e-5, iter = 20)
{
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  fn_jac_calculate = function() {
  dNt = drop(ifelse(time == stime[j], 1, 0))
  dNt = dNt*event
  Yt = drop(ifelse(time >= stime[j], 1, 0)) 
  res = rep(0, n)
  res_b = rep(0, n)
  res_a = matrix(0, nrow = n, ncol = ncol(Covariates))
  idx = Yt > 0
  if(sum(idx) == 0) break

  if(j == 1) {
      tmpDc = outer(IV[idx], IV[idx], "+")
      int_cbetaD = tmpDc * betaD * stime[j]
      if(betaD != 0){
          int_expbetaD = ifelse(IV[idx] == 1, (exp(betaD * int_D[idx, j]) - 1)/betaD, 
                                 stime[j])
          int_dexpbetaD = ifelse(IV[idx] == 1, int_D[idx, j] * exp(betaD * 
                                int_D[idx, j])/betaD, 0) - int_expbetaD/betaD
          int_cexpbetaD = ifelse(tmpDc > 0, 
          (exp(int_cbetaD)-1)/(tmpDc*betaD),
                             stime[j])
          int_cdexpbetaD = ifelse(tmpDc > 0, t(int_D[idx, j] * t(exp(int_cbetaD)))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((int_D[idx, j])^2/2, nrow = sum(idx), 
                                          ncol = sum(idx), byrow = TRUE))
      } else {
          int_expbetaD = stime[j]
          int_dexpbetaD = (int_D[idx, j])^2/2
          int_cexpbetaD = matrix(stime[j], nrow = sum(idx), ncol = sum(idx))
          int_cdexpbetaD = matrix((int_D[idx, j])^2/2, nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
      }
  } else {
      tmpDc = outer(D_status[idx, j-1], D_status[idx, j-1], "+")
      int_cbetaD = outer(int_D[idx, j], int_D[idx, j], "+") * betaD
      tmp_cbetaD = outer(int_D[idx, j-1], int_D[idx, j-1], "+") * betaD
      if(betaD != 0){
      int_expbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (exp(betaD * int_D[idx, j]) - 
                                    exp(betaD * int_D[idx, j-1]))/betaD, 
                                 (stime[j] - stime[j-1])*exp(betaD * int_D[idx, j]))
      int_dexpbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (int_D[idx, j]*exp(betaD * int_D[idx, j]) - 
                                    int_D[idx, j-1]*exp(betaD * int_D[idx, j-1]))/betaD,
                                    0) - int_expbetaD/betaD
      int_cexpbetaD = ifelse(tmpDc > 0, 
                             (exp(int_cbetaD) - 
                                exp(tmp_cbetaD))/(tmpDc*betaD),
                             (stime[j] - stime[j-1])*exp(int_cbetaD))
      int_cdexpbetaD = ifelse(tmpDc > 0, (t(int_D[idx, j] * t(exp(int_cbetaD))) - 
                                          t(int_D[idx, j-1] * t(exp(tmp_cbetaD))))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE))
    } else {
      int_expbetaD = stime[j] - stime[j-1]
      int_dexpbetaD = (stime[j] - stime[j-1]) * int_D[idx, j-1] + D_status[idx, j] * (stime[j] - stime[j-1])^2/2
      int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = sum(idx), ncol = sum(idx))
      int_cdexpbetaD = matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
    }
  }

  Covbeta = drop(Covariates %*% beta)
  SY = sum(exp(betaD*int_D[, j])*Yt)

  res_b1 = int_D[idx, j]*exp(betaD*int_D[idx, j])*dNt[idx]
  res_b2 = -int_dexpbetaD*(Covbeta[idx]+D_status[idx,j]*betaD)
  res_b3 = -int_expbetaD*D_status[idx,j]
  res_b4 = -(int_D[idx, j]*apply(exp(int_cbetaD)*dNt[idx], 2, sum) - 
              apply(int_cdexpbetaD*(D_status[idx,j]*betaD + Covbeta[idx]), 2, sum))/SY
  res_b5 = -((t(int_D[idx, j]*t(exp(betaD*int_D[idx, j])))*dNt[idx] - 
              apply(t(int_cdexpbetaD) * (Covbeta[idx]+D_status[idx,j]*betaD), 2, sum))*SY - 
              (apply(exp(int_cbetaD)*dNt[idx], 2, sum) - 
              apply(int_cexpbetaD*(Covbeta[idx] + D_status[idx, j]*betaD), 2, sum)) * 
              sum(int_D[idx, j]*exp(betaD*int_D[idx, j])*Yt))/SY^2

  res_b[idx] = res_b1 + res_b2 + res_b3 + res_b4 + res_b5

  res_a1 = -int_expbetaD*Covariates[idx, ]
  f = function(d){
    k = apply(Covariates[idx, , drop = F]*d,2,sum)
    k
  } 
  res_a2 = t(apply(int_cexpbetaD, 2, f))/SY

  res_a[idx, ] = res_a1 + res_a2

  k1 = c(sum(IV_c*res_b), apply(IV_c*res_a, 2, sum))
  k2 = cbind(t(Covariates) %*% res_b, t(Covariates)%*%res_a)
  Var = rbind(k1, k2)
        
  int_explam = drop(apply(exp(int_cbetaD)*drop(dNt[idx]), 2, sum))/SY - 
                  drop(apply(int_cexpbetaD*(D_status[idx, j]*betaD + 
                                  Covbeta[idx]), 2, sum))/SY
  res[idx] = dNt[idx]*exp(betaD*int_D[idx, j]) - 
              (D_status[idx, j]*betaD + 
                    Covbeta[idx])*int_expbetaD - 
                  int_explam

  out = cbind(c(sum(IV_c*res), drop(t(Covariates) %*% res)), Var)
  return(out)
}

  fn_calculate = function(){
  dNt = drop(ifelse(time == stime[j], 1, 0))
  dNt = dNt*event
  Yt = drop(ifelse(time >= stime[j], 1, 0)) 
  res = rep(0, n)
  idx = Yt > 0
  if(sum(idx) == 0) break
  if(j == 1) {
      tmpDc = outer(IV[idx], IV[idx], "+")
      int_cbetaD = tmpDc * betaD * stime[j]
      if(betaD != 0){
          int_expbetaD = ifelse(IV[idx] == 1, (exp(betaD * int_D[idx, j]) - 1)/betaD, 
                                 stime[j])
          int_dexpbetaD = ifelse(IV[idx] == 1, int_D[idx, j] * exp(betaD * 
                                int_D[idx, j])/betaD, 0) - int_expbetaD/betaD
          int_cexpbetaD = ifelse(tmpDc > 0, 
          (exp(int_cbetaD)-1)/(tmpDc*betaD),
                             stime[j])
          int_cdexpbetaD = ifelse(tmpDc > 0, t(int_D[idx, j] * t(exp(int_cbetaD)))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((int_D[idx, j])^2/2, nrow = sum(idx), 
                                          ncol = sum(idx), byrow = TRUE))
      } else {
          int_expbetaD = stime[j]
          int_dexpbetaD = (int_D[idx, j])^2/2
          int_cexpbetaD = matrix(stime[j], nrow = sum(idx), ncol = sum(idx))
          int_cdexpbetaD = matrix((int_D[idx, j])^2/2, nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
      }
  } else {
      tmpDc = outer(D_status[idx, j-1], D_status[idx, j-1], "+")
      int_cbetaD = outer(int_D[idx, j], int_D[idx, j], "+") * betaD
      tmp_cbetaD = outer(int_D[idx, j-1], int_D[idx, j-1], "+") * betaD
      if(betaD != 0){
      int_expbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (exp(betaD * int_D[idx, j]) - 
                                    exp(betaD * int_D[idx, j-1]))/betaD, 
                                 (stime[j] - stime[j-1])*exp(betaD * int_D[idx, j]))
      int_dexpbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (int_D[idx, j]*exp(betaD * int_D[idx, j]) - 
                                    int_D[idx, j-1]*exp(betaD * int_D[idx, j-1]))/betaD,
                                    0) - int_expbetaD/betaD
      int_cexpbetaD = ifelse(tmpDc > 0, 
                             (exp(int_cbetaD) - 
                                exp(tmp_cbetaD))/(tmpDc*betaD),
                             (stime[j] - stime[j-1])*exp(int_cbetaD))
      int_cdexpbetaD = ifelse(tmpDc > 0, (t(int_D[idx, j] * t(exp(int_cbetaD))) - 
                                          t(int_D[idx, j-1] * t(exp(tmp_cbetaD))))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE))
    } else {
      int_expbetaD = stime[j] - stime[j-1]
      int_dexpbetaD = (stime[j] - stime[j-1]) * int_D[idx, j-1] + D_status[idx, j] * (stime[j] - stime[j-1])^2/2
      int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = sum(idx), ncol = sum(idx))
      int_cdexpbetaD = matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
    }
  }

  Covbeta = drop(Covariates %*% beta)
  SY = sum(exp(betaD*int_D[, j])*Yt)
  int_explam = drop(apply(exp(int_cbetaD)*drop(dNt[idx]), 2, sum))/SY - 
                  drop(apply(int_cexpbetaD*(D_status[idx, j]*betaD + 
                                  Covbeta[idx]), 2, sum))/SY
  res[idx] = dNt[idx]*exp(betaD*int_D[idx, j]) - 
              (D_status[idx, j]*betaD + 
                    Covbeta[idx])*int_expbetaD - 
                  int_explam
  return(c(sum(IV_c*res), drop(t(Covariates) %*% res)))
}

  n = length(time)
  k = length(stime)
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))

  int_D = matrix(0, nrow = n, ncol = k)

  for(j in 1:k){
      if(j == 1){
      int_D[, j] = drop(IV * stime[j])
    }else{
      int_D[, j] = int_D[, j-1] + 
        drop(D_status[, j-1]*(stime[j] - stime[j-1]))
    }
  }

  covergence = FALSE
  for(i in 1:iter) {
    betaD = init_parameters[1]
    beta = init_parameters[-1]
    out = foreach(j = 1:k, .combine = "+") %dopar% {
        if(i == 1){
          out = fn_jac_calculate()
        } else {
          out = fn_calculate()
        }
        out
    }

    if(i == 1) {
      fn = drop(out[, 1])
      H = out[, -1]
      G = H
    } else {
      fn.new = out
      # browser()
      G = G + ((fn.new - fn - G %*% p)/sqrt(drop(t(p) %*% p))) %*% t(p)
    }
    
    p = solve(G)%*%fn
    tmp_parameters = init_parameters - p
    if(sum(abs(fn)) < tol | sum(abs(tmp_parameters - init_parameters)) < xtol){
        covergence = TRUE
        break
    }
    init_parameters = tmp_parameters
  }


  return(list(parameters = tmp_parameters, convergence = covergence, iter = i))
}


ConstantF_parallel_Newtonest = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime, tol = 1e-3, xtol = 1e-5, iter = 10)
{
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
    fn_jac_calculate = function() {
  dNt = drop(ifelse(time == stime[j], 1, 0))
  dNt = dNt*event
  Yt = drop(ifelse(time >= stime[j], 1, 0)) 
  res = rep(0, n)
  res_b = rep(0, n)
  res_a = matrix(0, nrow = n, ncol = ncol(Covariates))
  idx = Yt > 0
  if(sum(idx) == 0) break

  if(j == 1) {
      tmpDc = outer(IV[idx], IV[idx], "+")
      int_cbetaD = tmpDc * betaD * stime[j]
      if(betaD != 0){
          int_expbetaD = ifelse(IV[idx] == 1, (exp(betaD * int_D[idx, j]) - 1)/betaD, 
                                 stime[j])
          int_dexpbetaD = ifelse(IV[idx] == 1, int_D[idx, j] * exp(betaD * 
                                int_D[idx, j])/betaD, 0) - int_expbetaD/betaD
          int_cexpbetaD = ifelse(tmpDc > 0, 
          (exp(int_cbetaD)-1)/(tmpDc*betaD),
                             stime[j])
          int_cdexpbetaD = ifelse(tmpDc > 0, t(int_D[idx, j] * t(exp(int_cbetaD)))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((int_D[idx, j])^2/2, nrow = sum(idx), 
                                          ncol = sum(idx), byrow = TRUE))
      } else {
          int_expbetaD = stime[j]
          int_dexpbetaD = (int_D[idx, j])^2/2
          int_cexpbetaD = matrix(stime[j], nrow = sum(idx), ncol = sum(idx))
          int_cdexpbetaD = matrix((int_D[idx, j])^2/2, nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
      }
  } else {
      tmpDc = outer(D_status[idx, j-1], D_status[idx, j-1], "+")
      int_cbetaD = outer(int_D[idx, j], int_D[idx, j], "+") * betaD
      tmp_cbetaD = outer(int_D[idx, j-1], int_D[idx, j-1], "+") * betaD
      if(betaD != 0){
      int_expbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (exp(betaD * int_D[idx, j]) - 
                                    exp(betaD * int_D[idx, j-1]))/betaD, 
                                 (stime[j] - stime[j-1])*exp(betaD * int_D[idx, j]))
      int_dexpbetaD = ifelse(D_status[idx, j-1] == 1, 
                                 (int_D[idx, j]*exp(betaD * int_D[idx, j]) - 
                                    int_D[idx, j-1]*exp(betaD * int_D[idx, j-1]))/betaD,
                                    0) - int_expbetaD/betaD
      int_cexpbetaD = ifelse(tmpDc > 0, 
                             (exp(int_cbetaD) - 
                                exp(tmp_cbetaD))/(tmpDc*betaD),
                             (stime[j] - stime[j-1])*exp(int_cbetaD))
      int_cdexpbetaD = ifelse(tmpDc > 0, (t(int_D[idx, j] * t(exp(int_cbetaD))) - 
                                          t(int_D[idx, j-1] * t(exp(tmp_cbetaD))))/(tmpDc*betaD) - 
                                          int_cexpbetaD/(tmpDc*betaD), 
                                          matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE))
    } else {
      int_expbetaD = stime[j] - stime[j-1]
      int_dexpbetaD = (stime[j] - stime[j-1]) * int_D[idx, j-1] + D_status[idx, j] * (stime[j] - stime[j-1])^2/2
      int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = sum(idx), ncol = sum(idx))
      int_cdexpbetaD = matrix((stime[j] - stime[j-1]) * int_D[idx, j-1] + 
                                          D_status[idx, j] * (stime[j] - stime[j-1])^2/2, 
                                          nrow = sum(idx), ncol = sum(idx), byrow = TRUE)
    }
  }

  Covbeta = drop(Covariates %*% beta)
  SY = sum(exp(betaD*int_D[, j])*Yt)

  res_b1 = int_D[idx, j]*exp(betaD*int_D[idx, j])*dNt[idx]
  res_b2 = -int_dexpbetaD*(Covbeta[idx]+D_status[idx,j]*betaD)
  res_b3 = -int_expbetaD*D_status[idx,j]
  res_b4 = -(int_D[idx, j]*apply(exp(int_cbetaD)*dNt[idx], 2, sum) - 
              apply(int_cdexpbetaD*(D_status[idx,j]*betaD + Covbeta[idx]), 2, sum))/SY
  res_b5 = -((t(int_D[idx, j]*t(exp(betaD*int_D[idx, j])))*dNt[idx] - 
              apply(t(int_cdexpbetaD) * (Covbeta[idx]+D_status[idx,j]*betaD), 2, sum))*SY - 
              (apply(exp(int_cbetaD)*dNt[idx], 2, sum) - 
              apply(int_cexpbetaD*(Covbeta[idx] + D_status[idx, j]*betaD), 2, sum)) * 
              sum(int_D[idx, j]*exp(betaD*int_D[idx, j])*Yt))/SY^2

  res_b[idx] = res_b1 + res_b2 + res_b3 + res_b4 + res_b5

  res_a1 = -int_expbetaD*Covariates[idx, ]
  f = function(d){
    k = apply(Covariates[idx, , drop = F]*d,2,sum)
    k
  } 
  res_a2 = t(apply(int_cexpbetaD, 2, f))/SY

  res_a[idx, ] = res_a1 + res_a2

  k1 = c(sum(IV_c*res_b), apply(IV_c*res_a, 2, sum))
  k2 = cbind(t(Covariates) %*% res_b, t(Covariates)%*%res_a)
  Var = rbind(k1, k2)
        
  int_explam = drop(apply(exp(int_cbetaD)*drop(dNt[idx]), 2, sum))/SY - 
                  drop(apply(int_cexpbetaD*(D_status[idx, j]*betaD + 
                                  Covbeta[idx]), 2, sum))/SY
  res[idx] = dNt[idx]*exp(betaD*int_D[idx, j]) - 
              (D_status[idx, j]*betaD + 
                    Covbeta[idx])*int_expbetaD - 
                  int_explam

  out = cbind(c(sum(IV_c*res), drop(t(Covariates) %*% res)), Var)
  return(out)
}

  n = length(time)
  k = length(stime)
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))

  int_D = matrix(0, nrow = n, ncol = k)

  for(j in 1:k){
      if(j == 1){
      int_D[, j] = drop(IV * stime[j])
    }else{
      int_D[, j] = int_D[, j-1] + 
        drop(D_status[, j-1]*(stime[j] - stime[j-1]))
    }
  }

  covergence = FALSE
  for(i in 1:iter) {
    betaD = init_parameters[1]
    beta = init_parameters[-1]
    out = foreach(j = 1:k, .combine = "+") %dopar% {
        out = fn_jac_calculate()
        out
    }
    fn = drop(out[, 1])
    H = out[, -1]
    if(i == 1)print(H)
    
    p = solve(H)%*%fn
    tmp_parameters = init_parameters - p
    if(sum(abs(fn)) < tol | sum(abs(tmp_parameters - init_parameters)) < xtol){
        covergence = TRUE
        break
    }
    init_parameters = tmp_parameters
  }


  return(list(parameters = tmp_parameters, convergence = covergence, iter = i, H = H))
}