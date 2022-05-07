# accelerate computing using parallel

library(parallel)
library(doParallel)



expit = function(d) return(exp(d)/(exp(d)+1))

ConstantF_parallel = function(init_parameters, time, event, IV, 
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

  int_betaD = matrix(0, nrow = n, ncol = k)

  for(j in 1:k){
      if(j == 1){
      int_betaD[, j] = drop(IV * betaD * stime[j])
    }else{
      int_betaD[, j] = int_betaD[, j-1] + 
        drop(D_status[, j-1]*betaD*(stime[j] - stime[j-1]))
    }
  }


  out = foreach(j = 1:k, .combine = cbind) %dopar% {
    dNt = drop(ifelse(time == stime[j], 1, 0))
    dNt = dNt*event
    Yt = drop(ifelse(time >= stime[j], 1, 0)) 
    
    if(j == 1) {
        tmpDc = outer(IV, IV, "+")
        int_cbetaD = tmpDc * betaD * stime[j]
        if(betaD != 0){
            int_expbetaD = ifelse(IV == 1, (exp(int_betaD[, j]) - 1)/betaD, 
                                   stime[j])
            int_cexpbetaD = ifelse(tmpDc > 0, 
            (exp(int_cbetaD)-1)/(tmpDc*betaD),
                               stime[j])
        } else {
            int_expbetaD = stime[j]
            int_cexpbetaD = matrix(stime[j], nrow = n, ncol = n)
        }
    } else {
        tmpDc = outer(D_status[, j-1], D_status[, j-1], "+")
        int_cbetaD = outer(int_betaD[, j], int_betaD[, j], "+")
        tmp_cbetaD = outer(int_betaD[, j-1], int_betaD[, j-1], "+")
        if(betaD != 0){
        int_expbetaD = ifelse(D_status[, j-1] == 1, 
                                   (exp(int_betaD[, j]) - 
                                      exp(int_betaD[, j-1]))/betaD, 
                                   (stime[j] - stime[j-1])*exp(int_betaD[, j]))
        int_cexpbetaD = ifelse(tmpDc > 0, 
                               (exp(int_cbetaD) - 
                                  exp(tmp_cbetaD))/(tmpDc*betaD),
                               (stime[j] - stime[j-1])*exp(int_cbetaD))
      } else {
        int_expbetaD = stime[j] - stime[j-1]
        int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = n, ncol = n)
      }
    }

    SY = sum(exp(int_betaD[, j])*Yt)
    int_explam = Yt*drop(apply(exp(int_cbetaD)*drop(dNt), 2, sum))/SY - 
      Yt*drop(apply(int_cexpbetaD*Yt*(D_status[, j]*betaD + 
                                        drop(Covariates %*% beta)), 2, sum))/SY
    
    res = dNt*exp(int_betaD[, j]) - 
                    Yt*(D_status[, j]*betaD + 
                          drop(Covariates %*% beta))*int_expbetaD - 
      int_explam
    res
  }

  res = drop(apply(out, 1, sum))

  return(c(sum(IV_c*res), drop(t(Covariates) %*% res)))
}