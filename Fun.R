ConstantF = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime)
{
  # stime contains switching time
  # Covariates must be matrix
  # setting
  #----------------------------------------------------
  
  n = length(time)
  k = length(stime)
  betaD = init_parameters[1]
  beta = init_parameters[-1]
  time.order = order(-time)
  X = time[time.order]
  
  # tie.last = tie.first = 1:n
  # for(i in 2:n)
  # {
  #   if(dX[i-1]==0)
  #   {
  #     tie.first[i] = tie.first[i-1]
  #     tie.last[tie.first[i]:i] = i
  #   }
  # }
  # SY = cumsum(rep(1:n))[tie.last]
  
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
                                   (exp(betaD*stime[j]) - 
                                      exp(betaD*stime[j-1]))/betaD, 
                                   stime[j] - stime[j-1])
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

