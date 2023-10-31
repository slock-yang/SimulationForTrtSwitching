# This file is a wrap of Rcpp codes in directory "Function/integral" with function "integral_est_cpp"
# We use Newton's method to find the solution of equations, where a backtracking line search is also applied by its vulnerability.
# The function has inputs:
#
#   init_parameters:            initial value
#   time:                       Vector: censored time-to-event value
#   event:                      Vector: Indicator of censoring
#   IV:                         Vector: instrument variable
#   Covariates:                 Must be a matrix: the confounders in survival model
#   Covariates2:                A vector or a matrix: the confounders in propensity score
#   D_status:                   Matrix: the treatment therapy at each time point in stime
#   stime:                      The grid of time
#   max_iter:                   Max iterations allowed
#   tol:                        Tolerance
#   contraction:                Contractiong rate when using backtracking line search
#   eta:                        Threshold for backtracking


library(Rcpp)
sourceCpp("Function/integral/integral_est.cpp")


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
