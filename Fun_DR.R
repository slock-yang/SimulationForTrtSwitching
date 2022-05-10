expit = function(d) return(exp(d)/(exp(d)+1))

ConstantF_right = function(init_parameters, time, event, IV, 
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
  IV_c = IV
  
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

ConstantF_left = function(init_parameters, time, event, IV, 
                           Covariates, D_status, stime)
{
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  betaD = init_parameters[1]
  beta = rep(0, ncol(Covariates))
  mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
  IV_c = IV - expit(predict(mod))
  
  # calculus
  #----------------------------------------------------
  
  res = 0
  int_betaD = matrix(0, nrow = n, ncol = k)
  int_expbetaD = matrix(0, nrow = n, ncol = k)
  for(j in 1:k){
    dNt = drop(ifelse(time == stime[j], 1, 0))
    dNt = dNt*event
    Yt = drop(ifelse(time >= stime[j], 1, 0))
    
    if(j == 1){
      int_betaD[, j] = drop(IV * betaD * stime[j])
      if(betaD != 0) {
        int_expbetaD[, j] = ifelse(IV == 1, (exp(int_betaD[, j]) - 1)/betaD, 
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
    
    res = res + dNt*exp(int_betaD[, j]) - 
      Yt*(D_status[, j]*betaD)*int_expbetaD[,j]
  }
  
  return(sum(IV_c*res))
}



