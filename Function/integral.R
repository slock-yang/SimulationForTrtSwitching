library(Rcpp)
# sourceCpp("/ExampleforTrtSwitching/SimulationForTrtSwitching/Function/integral/integral.cpp")
sourceCpp("/Users/liuyang/Documents/Boston/MS\ SecondLine/Code\ &\ Analysis/DRcomparasion/ExampleforTrtSwitching/SimulationForTrtSwitching/Function/integral/integral_est.cpp")
# sourceCpp("/ExampleforTrtSwitching/SimulationForTrtSwitching/Function/integral/integral_est_debug.cpp")

# integral_cpp = function(init_parameters, time, event, IV, 
#     Covariates, D_status, stime) 
# {
#     betaD = init_parameters[1]
#     beta = init_parameters[-1]
#     mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
#     IV_c = IV - expit(predict(mod))
# 
#     fn = integral(betaD = betaD, beta = beta, time = time, event = event, IV = IV,
#                 IV_c = IV_c, Covariates = Covariates, D_status = D_status, stime = stime)
# 
#     return(fn)
# }


integral_est_cpp = function(init_parameters, time, event, IV, 
    Covariates, Covariates2, D_status, stime, max_iter = 20, tol = 1e-5,
    contraction = 0.5, eta = 1e-4)
{
    mod = glm(IV ~ Covariates2, family = binomial(link = "logit"))
    IV_c = IV - expit(predict(mod))

    out = integral_est(init_parameters = init_parameters, time = time, 
            event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
            D_status = D_status, stime = stime, max_iter = max_iter,
            tol = tol, contraction = contraction, eta = eta)

    return(out)
}


# integral_est_debug_cpp = function(init_parameters, time, event, IV, 
#     Covariates, D_status, stime, max_iter = 20, tol = 1e-5)
# {
#     mod = glm(IV ~ Covariates, family = binomial(link = "logit"))
#     IV_c = IV - expit(predict(mod))
# 
#     out = integral_est_debug(init_parameters = init_parameters, time = time, 
#             event = event, IV = IV, IV_c = IV_c, Covariates = Covariates, 
#             D_status = D_status, stime = stime, max_iter = max_iter,
#             tol = tol)
# 
#     return(out)
# }