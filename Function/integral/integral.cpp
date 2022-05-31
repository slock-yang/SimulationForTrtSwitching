#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP integral(
    double betaD,
    NumericVector beta,
    NumericVector time,
    NumericVector event,
    NumericVector IV,
    NumericVector IV_c,
    NumericMatrix Covariates,
    NumericMatrix D_status,
    NumericVector stime
)
{
    int n = time.size();
    int k = stime.size();
    int p = beta.size();
    
    NumericVector res(n), Covbeta(n), int_cexpbetaD(n), fn(p+1), SY(k);
    NumericMatrix int_D(n, k), dNt(n, k), Yt(n, k), int_expbetaD(n, k), int_Lam(n, k);

    for (int i = 0; i < n; i++)
    {
        for (int ii = 0; ii < p; ii++)
        {
            Covbeta[i] += Covariates(i, ii) * beta[ii];
        }
    }
    
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
                        int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
                    else 
                        int_expbetaD(i, j) = stime[j];
                }
                else
                {
                    if (betaD != 0)
                        int_expbetaD(i, j) = (D_status(i, j-1) > 0.5) ? (exp(betaD * int_D(i, j)) - exp(betaD * int_D(i, j-1)))/betaD:(stime[j] - stime[j-1]) * exp(betaD * int_D(i, j));
                    else 
                        int_expbetaD(i, j) = stime[j] - stime[j-1];
                }
                

                for (int ii = 0; ii < n; ii++)
                {
                    if (j == 0)
                    {
                        if (betaD != 0)
                            int_cexpbetaD[ii] = ((IV[i] + IV[ii]) > 0.5) ? (exp((int_D(i, j) + int_D(ii, j)) * betaD) - 1) / ((IV[i] + IV[ii]) * betaD):stime[j];
                        else
                            int_cexpbetaD[ii] = stime[j];
                    }
                    else
                    {
                        if (betaD != 0)
                            int_cexpbetaD[ii] = ((D_status(i, j-1) + D_status(ii, j-1)) > 0.5) ? (exp((int_D(i, j) + int_D(ii, j)) * betaD) - exp((int_D(i, j-1) + int_D(ii, j-1)) * betaD)) / ((D_status(i, j-1) + D_status(ii, j-1)) * betaD):(stime[j] - stime[j-1]) * exp((int_D(i, j) + int_D(ii, j))*betaD);
                        else
                            int_cexpbetaD[ii] = stime[j] - stime[j-1];
                    }
                    int_Lam(i, j) = int_Lam(i, j) +  dNt(ii, j) * exp((int_D(i, j) + int_D(ii, j)) * betaD) - Yt(ii, j) * int_cexpbetaD[ii] * (D_status(ii, j) * betaD + Covbeta[ii]);
                }

                int_Lam(i, j) = int_Lam(i, j) / SY[j];
            }
                   
            res[i] = res[i] + exp(int_D(i, j) * betaD) * dNt(i, j) - Yt(i, j) * int_expbetaD(i, j) * (D_status(i, j) * betaD + Covbeta[i]) - Yt(i, j) * int_Lam(i, j);
        }


        for (int kk = 0; kk < p+1; kk++)
        {
            if (kk == 0)
            {
                fn[kk] = fn[kk] + IV_c[i] * res[i];
            }
            else
            {
                fn[kk] = fn[kk] + Covariates(i, kk-1) * res[i];
            }
        }
    }

    return fn;
}