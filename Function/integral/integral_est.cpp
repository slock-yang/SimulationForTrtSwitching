#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP integral_est(
    arma::vec init_parameters,
    arma::vec time,
    arma::vec event,
    arma::vec IV,
    arma::vec IV_c,
    arma::mat Covariates,
    arma::mat D_status,
    arma::vec stime,
    int max_iter,
    double tol
)
{
    int n = time.size();
    int k = stime.size();
    int p = init_parameters.size() - 1;
    // printf("%d %d %d", n, k, p);
    
    bool Convergence;
    double tmp, tmp_cexpbetaD, tmp_cD, old_cexpbetaD, betaD;
    arma::vec beta(p), res(n), Covbeta(n), int_cexpbetaD(n), int_cdexpbetaD_1(n), int_cdexpbetaD_2(n), fn(p+1), fn_abs(p+1), SY(k), dSY(k), new_parameters(p+1), diff(p+1);
    arma::mat int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_Lam(n, k), int_dexpbetaD(n, k), int_dexpbetaD_Lam(n, k), int_expbetaD_dLam(n, k), int_expbetaD_dLam_dalpha(n, p), dPhi(n, p+1), Hessian(p+1, p+1), inv_Hessian(p+1, p+1);

    for (int iter = 0; iter < max_iter; iter++)
    {
        betaD = init_parameters[0];
        beta = init_parameters[-0];
        for (int i = 0; i < n; i++)
        {
            for (int ii = 0; ii < p; ii++)
            {
                Covbeta[i] += Covariates(i, ii) * beta[ii];
            }
        }

        // arma::arma_print(Covbeta);
    
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                if(j == 0)
                    int_D(i, j) = IV[i] * stime[j];
                else 
                    int_D(i, j) = int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]);

                dNt(i, j) = (time[i] == stime[j]) ? 1:0;
                dNt(i, j) *= event[i];
                Yt(i, j) = (time[i] >= stime[j]) ? 1:0;
                SY[j] = SY[j] + exp(betaD * int_D(i, j)) * Yt(i, j);
                dSY[j] = dSY[j] + int_D(i, j) * exp(betaD * int_D(i, j)) * Yt(i, j);
            }
        }

        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                for (int kk = 0; kk < p; kk++)
                {
                    int_expbetaD_dLam_dalpha(i, kk) = 0;
                }
                
                if (Yt(i, j) > 0)
                {
                    if (j == 0)
                    {
                        if (betaD != 0)
                        {
                            int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
                            int_dexpbetaD(i, j) = (IV[i] > 0.5) ? (int_D(i, j) * exp(betaD * int_D(i, j))/betaD - int_expbetaD(i, j)/betaD):(- int_expbetaD(i, j)/betaD);
                        }
                        else
                        { 
                            int_expbetaD(i, j) = stime[j];
                            int_dexpbetaD(i, j) = int_D(i, j) * int_D(i, j) / 2;
                        }
                    }
                    else
                    {
                        if (betaD != 0)
                        {
                            int_expbetaD(i, j) = (D_status(i, j-1) > 0.5) ? (exp(betaD * int_D(i, j)) - exp(betaD * int_D(i, j-1)))/betaD:(stime[j] - stime[j-1]) * exp(betaD * int_D(i, j));
                            int_dexpbetaD(i, j) = (D_status(i, j-1) > 0.5) ? ((int_D(i, j) * exp(betaD * int_D(i, j)) - int_D(i, j-1) * exp(betaD * int_D(i, j-1)))/betaD - int_expbetaD(i, j)/betaD) : ( - int_expbetaD(i, j)/betaD);
                        }
                        else 
                        {
                            int_expbetaD(i, j) = stime[j] - stime[j-1];
                            int_dexpbetaD(i, j) = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                            
                        }
                    }
                    

                    for (int ii = 0; ii < n; ii++)
                    {
                        tmp_cexpbetaD = exp((int_D(i, j) + int_D(ii, j)) * betaD);
                        if (j == 0)
                        {
                            // printf("%d\n", j);
                            tmp_cD = IV[i] + IV[ii];
                            if (betaD != 0)
                            {
                                int_cexpbetaD[ii] = ((tmp_cD) > 0.5) ? (tmp_cexpbetaD - 1) / ((tmp_cD) * betaD):stime[j];
                                int_cdexpbetaD_1[ii] = (tmp_cD > 0.5) ? (int_D(i, j) * tmp_cexpbetaD/((tmp_cD)*betaD) - int_cexpbetaD[ii]/((tmp_cD)*betaD)):(int_D(i, j) * int_D(i, j)/2);
                                int_cdexpbetaD_2[ii] = (tmp_cD > 0.5) ? (int_D(ii, j) * tmp_cexpbetaD/((tmp_cD)*betaD) - int_cexpbetaD[ii]/((tmp_cD)*betaD)):(int_D(ii, j) * int_D(ii, j)/2);
                            }
                            else
                            {
                                int_cexpbetaD[ii] = stime[j];
                                int_cdexpbetaD_1[ii] = int_D(i, j) * int_D(i, j)/2;
                                int_cdexpbetaD_2[ii] = int_D(ii, j) * int_D(ii, j)/2;
                            }
                        }
                        else
                        {
                            
                            old_cexpbetaD = exp((int_D(i, j-1) + int_D(ii, j-1)) * betaD);
                            tmp_cD = D_status(i, j-1) + D_status(ii, j-1);
                            if (betaD != 0)
                            {
                                int_cexpbetaD[ii] = ((tmp_cD) > 0.5) ? 
                                                    (tmp_cexpbetaD - old_cexpbetaD) / ((tmp_cD) * betaD):
                                                    (stime[j] - stime[j-1]) * exp((int_D(i, j) + int_D(ii, j))*betaD);
                                int_cdexpbetaD_1[ii] = (tmp_cD > 0.5) ? 
                                                    ((int_D(i, j) * tmp_cexpbetaD - int_D(i, j-1) * old_cexpbetaD)/((tmp_cD)*betaD) - 
                                                    int_cexpbetaD[ii]/((tmp_cD)*betaD)):((stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2);
                                int_cdexpbetaD_2[ii] = (D_status(ii, j-1) > 0.5) ? 
                                                    ((int_D(ii, j) * tmp_cexpbetaD - int_D(ii, j-1) * old_cexpbetaD)/((tmp_cD)*betaD) - 
                                                    int_cexpbetaD[ii]/((tmp_cD)*betaD)):((stime[j] - stime[j-1]) * int_D(ii, j-1) + D_status(ii, j) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2);
                            }
                            else
                            {
                                // printf("%d\n", j);
                                int_cexpbetaD[ii] = stime[j] - stime[j-1];
                                int_cdexpbetaD_1[ii] = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                                int_cdexpbetaD_2[ii] = (stime[j] - stime[j-1]) * int_D(ii, j-1) + D_status(ii, j) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                            }
                        }
                        tmp = dNt(ii, j) * tmp_cexpbetaD - Yt(ii, j) * int_cexpbetaD[ii] * (D_status(ii, j) * betaD + Covbeta[ii]);
                        int_Lam(i, j) = int_Lam(i, j) +  tmp;
                        int_dexpbetaD_Lam(i, j)  = int_dexpbetaD_Lam(i, j) + dNt(ii, j) * int_D(i, j) * tmp_cexpbetaD - 
                                                        Yt(ii, j) * int_cdexpbetaD_1[ii] * (D_status(ii, j) * betaD + Covbeta[ii]);
                        int_expbetaD_dLam(i, j) = int_expbetaD_dLam(i, j) + (dNt(ii, j) * int_D(ii, j) * tmp_cexpbetaD - 
                                                        Yt(ii, j) * int_cdexpbetaD_2[ii] * (D_status(ii, j) * betaD + Covbeta[ii]) - 
                                                        int_cexpbetaD[ii] * Yt(ii, j) * D_status(ii, j)) * SY[j] - tmp * dSY[j];                        
                        for (int kk = 0; kk < p; kk++)
                        {
                            int_expbetaD_dLam_dalpha(i, kk) = int_expbetaD_dLam_dalpha(i, kk) + int_cexpbetaD[ii] * Yt(ii, j) * Covariates(ii, kk)/SY[j];
                        }
                    }
                    
                    int_Lam(i, j) = int_Lam(i, j) / SY[j];
                    int_dexpbetaD_Lam(i, j) = int_dexpbetaD_Lam(i, j) / SY[j];
                    int_expbetaD_dLam(i, j) = int_expbetaD_dLam(i, j) / (SY[j] * SY[j]);

                }
                res[i] = res[i] + exp(int_D(i, j) * betaD) * dNt(i, j) - Yt(i, j) * int_expbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i]) - Yt(i, j) * int_Lam(i, j);

                for (int kk = 0; kk < p+1; kk++)
                {
                    if (kk == 0)
                    {
                        
                        dPhi(i, kk) = dPhi(i, kk) + dNt(i, j) * int_D(i, j) * exp(betaD * int_D(i, j)) - 
                                        Yt(i, j) * int_dexpbetaD[i] * (D_status(i, j) * betaD + Covbeta[i]) - 
                                        Yt(i, j) * int_dexpbetaD_Lam(i, j) - 
                                        Yt(i, j) * int_expbetaD(i, j) * D_status(i, j) - 
                                        Yt(i, j) * int_expbetaD_dLam(i, j);
                    }
                    else
                    {
                        dPhi(i, kk) = dPhi(i, kk) - Yt(i, j) * Covariates(i, kk-1) * int_expbetaD(i, j) - 
                                        Yt(i, j) * int_expbetaD_dLam_dalpha(i, kk-1);
                    }
                    
                }
                
            }
            // printf("%d\n", i);
            for (int kk = 0; kk < p+1; kk++)
            {
                if (kk == 0)
                {
                    fn[kk] = fn[kk] + IV_c[i] * res[i];
                    
                    for (int kkk = 0; kkk < p+1; kkk++)
                    {
                        Hessian(kk, kkk) = Hessian(kk, kkk) + IV_c[i] * dPhi(i, kkk);                        
                    }
                }
                else
                {
                    fn[kk] = fn[kk] + Covariates(i, kk-1) * res[i];
                    for (int kkk = 0; kkk < p; kkk++)
                    {
                        Hessian(kk, kkk) = Hessian(kk, kkk) + Covariates(i, kkk) * dPhi(i, kkk);
                    }
                }
            }
        }
        
        new_parameters = init_parameters - arma::solve(Hessian, fn);
        fn_abs = arma::abs(fn);
        diff = arma::abs(new_parameters - init_parameters);
        if(arma::sum(fn) < tol || arma::sum(diff) < (tol * tol))
        {
            Convergence = true;
            break;
        }
        init_parameters = new_parameters;
        if (iter >= (max_iter - 1))
        {
            Convergence = false;
        }
    }

    return Rcpp::List::create(
        Named("x") = init_parameters,
        Named("Convergence") = Convergence,
        Named("Hessian") = Hessian,
        Named("dPhi") = dPhi,
        Named("int_dexpbetaD") = int_dexpbetaD,
        Named("int_dexpbetaD_Lam") = int_dexpbetaD_Lam,
        Named("int_expbetaD_dLam") = int_expbetaD_dLam
    );
}