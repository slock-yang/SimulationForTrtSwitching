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
  
  # calculus
  #----------------------------------------------------
  
  res = 0
  int_betaD = matrix(0, nrow = n, ncol = k)
  int_expbetaD = matrix(0, nrow = n, ncol = k)
  for(j in 1:k){
    if(j == 1){
      int_betaD[, j] = drop(Z * betaD * stime[j])
      if(betaD != 0) {
        int_expbetaD[, j] = ifelse(Z == 1, (exp(int_betaD[, j]) - 1)/betaD, 
                                   stime[j])
      } else {
        int_expbetaD[, j] = stime[j]
      }
    }else{
      int_betaD[, j] = int_betaD[, j-1] + 
        drop(D_status[, j-1]*betaD*(stime[j] - stime[j-1]))
      if(betaD != 0){
        int_expbetaD[, j] = ifelse(D_status[, j-1] == 1, 
                                   (exp(int_betaD[, j]) - 
                                      exp(int_betaD[, j-1]))/betaD, 
                                   (stime[j] - stime[j-1])*exp(int_betaD[, j]))
      } else {
        int_expbetaD[, j] = stime[j] - stime[j-1]
      }
    }
    
    dNt = ifelse(time == stime[j], 1, 0)
    dNt = dNt*event
    Yt = ifelse(time >= stime[j], 1, 0)
    res = res + dNt*exp(int_betaD[, j]) - 
                    Yt*(D_status[, j]*betaD + 
                          0.25+drop(Covariates %*% beta))*int_expbetaD[,j]
  }
  
  return(c(sum(Z*res), drop(t(Covariates) %*% res)))
}



ConstantF_est = function(init_parameters, time, event, IV, 
                         Covariates, D_status, stime, maxit=100, tol=1e-3)
{
  # stime contains switching time
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  
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
    score = c(sum(Z*res), drop(t(Covariates) %*% res))
    dscore1 = c(sum(Z*resd1), drop(apply(Z*resd2, 2, sum)))
    dscore2 = cbind(t(Covariates)%*%resd1, t(Covariates)%*%resd2)
    dscore = rbind(dscore1, dscore2)
    delta = drop(solve(dscore)%*%score)
    init_parameters = init_parameters - delta
    if(sum(abs(score)) < tol | sum(abs(delta)) < tol) break
  }
  
  return(c(init_parameters,sum(abs(score)), k))
}



