#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP integral_est_debug(
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
    // arma::arma_print(IV[10]);
    
    bool Convergence;
    double betaD;
    arma::vec beta(p), res(n), Covbeta(n), fn(p+1), fn_abs(p+1), new_parameters(p+1), diff(p+1);
    arma::mat int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_dexpbetaD(n, k), dPhi(n, p+1), Hessian(p+1, p+1), inv_Hessian(p+1, p+1);

    for (int iter = 0; iter < max_iter; iter++)
    {
        Covbeta.zeros();
        res.zeros();
        dPhi.zeros();
        Hessian.zeros();
        fn.zeros();
        int_D.zeros();
        int_expbetaD.zeros();
        int_dexpbetaD.zeros();
        
        betaD = init_parameters[0];
        beta = init_parameters[1];
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
            }
        }

        
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < k; j++)
            {
                if (Yt(i, j) > 0)
                {
                    if (j == 0)
                    {
                        if (betaD != 0)
                        {
                            int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
                            int_dexpbetaD(i, j) = (IV[i] > 0.5) ? (int_D(i, j) * exp(betaD * int_D(i, j))/betaD - int_expbetaD(i, j)/betaD):0;
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
                            int_dexpbetaD(i, j) = (D_status(i, j-1) > 0.5) ? ((int_D(i, j) * exp(betaD * int_D(i, j)) - int_D(i, j-1) * exp(betaD * int_D(i, j-1)))/betaD - int_expbetaD(i, j)/betaD) : (int_D(i, j)*exp(betaD * int_D(i, j))*(stime[j] - stime[j-1]));
                        }
                        else 
                        {
                            int_expbetaD(i, j) = stime[j] - stime[j-1];
                            int_dexpbetaD(i, j) = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
                            
                        }
                    }
                }
                res[i] = res[i] + exp(int_D(i, j) * betaD) * dNt(i, j) - Yt(i, j) * int_expbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i] + 0.25);

                for (int kk = 0; kk < p+1; kk++)
                {
                    if (kk == 0)
                    {
                        
                        dPhi(i, kk) = dPhi(i, kk) + dNt(i, j) * int_D(i, j) * exp(betaD * int_D(i, j)) - 
                                        Yt(i, j) * int_dexpbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i] + 0.25) - 
                                        Yt(i, j) * int_expbetaD(i, j) * D_status(i, j);
                    }
                    else
                    {
                        dPhi(i, kk) = dPhi(i, kk) - Yt(i, j) * Covariates(i, kk-1) * int_expbetaD(i, j);
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
                    for (int kkk = 0; kkk < p+1; kkk++)
                    {
                        Hessian(kk, kkk) = Hessian(kk, kkk) + Covariates(i, kk-1) * dPhi(i, kkk);
                    }
                }
            }
        }
        
        new_parameters = init_parameters - arma::solve(Hessian, fn);
        fn_abs = arma::abs(fn);
        diff = arma::abs(new_parameters - init_parameters);
        if(arma::sum(fn_abs) < tol || arma::sum(diff) < (tol * tol))
        {
            printf("%d\n", iter);
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
        Named("Convergence") = Convergence
    );
}