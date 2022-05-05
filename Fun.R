expit = function(d) return(exp(d)/(exp(d)+1))

ConstantF = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime)
{
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  betaD = init_parameters[1]
  beta = init_parameters[-1]
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))
  
  # calculus
  #----------------------------------------------------
  
  res = 0
  int_betaD = matrix(0, nrow = n, ncol = k)
  int_cbetaD = matrix(0, nrow = n, ncol = n)
  int_expbetaD = matrix(0, nrow = n, ncol = k)
  int_cexpbetaD = matrix(0, nrow = n, ncol = n)
  for(j in 1:k){
    dNt = drop(ifelse(time == stime[j], 1, 0))
    dNt = dNt*event
    Yt = drop(ifelse(time >= stime[j], 1, 0))
    
    if(j == 1){
      tmpDc = outer(IV, IV, "+")
      int_betaD[, j] = drop(IV * betaD * stime[j])
      int_cbetaD = tmpDc * betaD * stime[j]
      if(betaD != 0) {
        int_expbetaD[, j] = ifelse(IV == 1, (exp(int_betaD[, j]) - 1)/betaD, 
                                   stime[j])
        int_cexpbetaD = ifelse(tmpDc > 0, (exp(int_cbetaD)-1)/(tmpDc*betaD),
                               stime[j])
      } else {
        int_expbetaD[, j] = stime[j]
        int_cexpbetaD = matrix(stime[j], nrow = n, ncol = n)
      }
    }else{
      tmpDc = outer(D_status[, j-1], D_status[, j-1], "+")
      int_betaD[, j] = int_betaD[, j-1] + 
        drop(D_status[, j-1]*betaD*(stime[j] - stime[j-1]))
      new_int_cbetaD = int_cbetaD + tmpDc * betaD * (stime[j] - stime[j-1])
      if(betaD != 0){
        int_expbetaD[, j] = ifelse(D_status[, j-1] == 1, 
                                   (exp(int_betaD[, j]) - 
                                      exp(int_betaD[, j-1]))/betaD, 
                                   (stime[j] - stime[j-1])*exp(int_betaD[, j]))
        int_cexpbetaD = ifelse(tmpDc > 0, 
                               (exp(new_int_cbetaD) - 
                                  exp(int_cbetaD))/(tmpDc*betaD),
                               (stime[j] - stime[j-1])*exp(new_int_cbetaD))
        int_cbetaD = new_int_cbetaD
      } else {
        int_expbetaD[, j] = stime[j] - stime[j-1]
        int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = n, ncol = n)
        int_cbetaD = new_int_cbetaD
      }
    }
    
    SY = sum(exp(int_betaD[, j])*Yt)
    int_explam = Yt*drop(apply(exp(int_cbetaD)*drop(dNt), 2, sum))/SY - 
      Yt*drop(apply(int_cexpbetaD*Yt*(D_status[, j]*betaD + 
                                        drop(Covariates %*% beta)), 2, sum))/SY
    
    res = res + dNt*exp(int_betaD[, j]) - 
                    Yt*(D_status[, j]*betaD + 
                          drop(Covariates %*% beta))*int_expbetaD[,j] - 
      int_explam
  }
  
  return(c(sum(IV_c*res), drop(t(Covariates) %*% res)))
}



ConstantF_est = function(init_parameters, time, event, Z, 
                         Covariates, D_status, stime, maxit=100, tol=1e-3)
{
  # stime contains switching time
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  mod = glm(Z~Covariates, family = binomial(link = "logit"))
  Z_c = Z - expit(predict(mod))
  
  
  # calculus
  #----------------------------------------------------
  
  res = 0
  resd1 = 0
  resd2 = 0
  int_D = matrix(0, nrow = n, ncol = k)
  int_betaD = matrix(0, nrow = n, ncol = k)
  int_expbetaD = matrix(0, nrow = n, ncol = k)
  int_texpbetaD = matrix(0, nrow = n, ncol = k)
  for(k in 1:maxit){
    for(j in 1:k){
      betaD = init_parameters[1]
      beta = init_parameters[-1]
      if(j == 1){
        int_D[, j] = drop(Z*stime[j])
        int_betaD[, j] = drop(betaD*int_D[, j])
        if(betaD != 0) {
          int_expbetaD[, j] = ifelse(Z == 1, (exp(int_betaD[, j]) - 1)/betaD, 
                                     stime[j])
          int_texpbetaD[, j] = ifelse(Z == 1, 
                                      stime[j]*exp(int_betaD[, j])/betaD,
                                      0) - int_expbetaD[, j]
        } else {
          int_expbetaD[, j] = stime[j]
          int_texpbetaD[, j] = ifelse(Z == 1, stime[j], 0) - int_expbetaD[, j]
        }
        
      }else{
        int_D[, j] = int_D[, j-1] + D_status[, j-1]*(stime[j] - stime[j-1])
        int_betaD[, j] = int_betaD[, j-1] + 
          drop(D_status[, j-1]*betaD*(stime[j] - stime[j-1]))
        if(betaD != 0){
          int_expbetaD[, j] = ifelse(D_status[, j-1] == 1, 
                                     (exp(int_betaD[, j]) - 
                                        exp(int_betaD[, j-1]))/betaD, 
                                     (stime[j] - stime[j-1])*exp(int_betaD[,j]))
          int_texpbetaD[, j] = ifelse(D_status[, j-1] == 1,
                        (int_D[, j]*exp(int_betaD[, j])-
                           int_D[, j-1]*exp(int_betaD[, j-1]))/betaD,
                        0) - int_expbetaD[, j]
        } else {
          int_expbetaD[, j] = stime[j] - stime[j-1]
          int_texpbetaD[, j] = ifelse(D_status[, j] == 1, stime[j]-stime[j-1],
                                      0) - int_expbetaD[, j]
        }
      }
      
      dNt = ifelse(time == stime[j], 1, 0)
      dNt = dNt*event
      Yt = ifelse(time >= stime[j], 1, 0)
      res = res + dNt*exp(int_betaD[, j]) - 
        Yt*(D_status[, j]*betaD + 
              0.25+drop(Covariates %*% beta))*int_expbetaD[,j]
      resd1 = resd1 + int_D[, j]*exp(int_betaD[, j])*dNt -
        Yt*(D_status[, j]*betaD + 0.25 + 
              drop(Covariates %*% beta))*int_texpbetaD[, j] -
        Yt*D_status[, j]*int_expbetaD[, j]
      resd2 = resd2 - Yt*Covariates*int_expbetaD[, j]
    }
    score = c(sum(Z_c*res), drop(t(Covariates) %*% res))
    dscore1 = c(sum(Z_c*resd1), drop(apply(Z_c*resd2, 2, sum)))
    dscore2 = cbind(t(Covariates)%*%resd1, t(Covariates)%*%resd2)
    dscore = rbind(dscore1, dscore2)
    delta = drop(solve(dscore)%*%score)
    init_parameters = init_parameters - delta
    if(sum(abs(score)) < tol | sum(abs(delta)) < tol) break
  }
  
  return(c(init_parameters,sum(abs(score)), k))
}



SCSM = function(time, event, IV, 
    Covariates, D_status, stime)
{
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  p = ncol(Covariates)
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))


  # calculus
  #----------------------------------------------------

  res = 0
  int_betaD = matrix(0, nrow = n, ncol = k)
  int_expbetaD = matrix(0, nrow = n, ncol = k)
  dBD = rep(0, k)
  beta = matrix(0, nrow = p, ncol = k)
  init = c(0, rep(0, p))


  f = function(init){
    dBD = init[1]
    beta = init[-1]
    if(j == 1){
      span = stime[j]
    } else{
      span = stime[j] - stime[j-1]
    }

    res = exp(int_betaD[, j])*dNt - exp(int_betaD[,j])*Yt*(
      D_status[, j]*dBD + drop(Covariates%*%beta) + 0.25*span
    )
    return(c(sum(IV_c*res), drop(t(Covariates) %*% res)))
  }


  for(j in 1:k){
    dNt = drop(ifelse(time == stime[j], 1, 0))
    dNt = dNt*event
    Yt = drop(ifelse(time >= stime[j], 1, 0))

    if(j == 1){
      int_betaD[, j] = drop(IV * 0)
    }else{
      int_betaD[, j] = int_betaD[, j-1] + 
        drop(D_status[, j-1]*dBD[j-1])
    }
    
    s = nleqslv(init, f)
    dBD[j] = s$x[1]
    beta[, j] = s$x[-1]
    # cat("iter ", j, " ",s$message, "\n", sep = "")
  }

  return(list(dBD = dBD, BD = cumsum(dBD), dbeta = beta, 
  beta = t(apply(beta, 1, cumsum))))
}
