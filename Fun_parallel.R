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
    res = rep(0, n)
    idx = Yt > 0
    if(sum(idx) == 0) break
    if(j == 1) {
        tmpDc = outer(IV[idx], IV[idx], "+")
        int_cbetaD = tmpDc * betaD * stime[j]
        if(betaD != 0){
            int_expbetaD = ifelse(IV[idx] == 1, (exp(int_betaD[idx, j]) - 1)/betaD, 
                                   stime[j])
            int_cexpbetaD = ifelse(tmpDc > 0, 
            (exp(int_cbetaD)-1)/(tmpDc*betaD),
                               stime[j])
        } else {
            int_expbetaD = stime[j]
            int_cexpbetaD = matrix(stime[j], nrow = sum(idx), ncol = sum(idx))
        }
    } else {
        tmpDc = outer(D_status[idx, j-1], D_status[idx, j-1], "+")
        int_cbetaD = outer(int_betaD[idx, j], int_betaD[idx, j], "+")
        tmp_cbetaD = outer(int_betaD[idx, j-1], int_betaD[idx, j-1], "+")
        if(betaD != 0){
        int_expbetaD = ifelse(D_status[idx, j-1] == 1, 
                                   (exp(int_betaD[idx, j]) - 
                                      exp(int_betaD[idx, j-1]))/betaD, 
                                   (stime[j] - stime[j-1])*exp(int_betaD[idx, j]))
        int_cexpbetaD = ifelse(tmpDc > 0, 
                               (exp(int_cbetaD) - 
                                  exp(tmp_cbetaD))/(tmpDc*betaD),
                               (stime[j] - stime[j-1])*exp(int_cbetaD))
      } else {
        int_expbetaD = stime[j] - stime[j-1]
        int_cexpbetaD = matrix(stime[j] - stime[j-1], nrow = sum(idx), ncol = sum(idx))
      }
    }

    Covbeta = drop(Covariates %*% beta)

    SY = sum(exp(int_betaD[, j])*Yt)
    int_explam = drop(apply(exp(int_cbetaD)*drop(dNt[idx]), 2, sum))/SY - 
      drop(apply(int_cexpbetaD*(D_status[idx, j]*betaD + 
                                        Covbeta[idx]), 2, sum))/SY
    
    res[idx] = dNt[idx]*exp(int_betaD[idx, j]) - 
                    (D_status[idx, j]*betaD + 
                          Covbeta[idx])*int_expbetaD - 
      int_explam
    res
  }

  res = drop(apply(out, 1, sum))

  return(c(sum(IV_c*res), drop(t(Covariates) %*% res)))
}