library(Rcpp)
sourceCpp("Function/itegral/integral.cpp")

integral_cpp = function(init_parameters, time, event, IV, 
    Covariates, D_status, stime) 
{
    betaD = init_parameters[1]
    beta = init_parameters[-1]
    mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
    IV_c = IV - expit(predict(mod))

    fn = integral(betaD = betaD, beta = beta, time = time, event = event, IV = IV,
                IV_c = IV_c, Covariates = Covariates, D_status = D_status, stime = stime)

    return(fn)
}