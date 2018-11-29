//#include <Rcpp.h>
#include "auxmix.h"
#include "wlapack.h"

//' @export
//[[Rcpp::export]]
Rcpp::NumericVector draw_hX(Rcpp::IntegerVector s, Rcpp::NumericVector y, double phi, double sigma2, double mu, Rcpp::NumericMatrix X, Rcpp::NumericVector beta){

  int T = s.length();
  Rcpp::NumericVector y_star = log(y * y );
  Rcpp::NumericVector gamma(T+1);
  Rcpp::NumericVector H((T+1)*2);
  Rcpp::NumericVector Omega((T+1)*2);

  // set up gamma
   gamma[0] = mu;
   for(int i = 1; i < T+1; i++){
     gamma[i] = mu * (1 - phi);
   }

   // X*beta + gamma
   dgemv(X, beta, gamma);

   //H'Sig_u^-1 (X*beta + gamma)
   UBMV(gamma, phi, sigma2);

   //Sig_y^-1 (y-d) + H'Sig_u^-1 (X*beta + gamma)
   for(int i = 0; i < T; i++){
     gamma[i+1] = gamma[i+1] + (y_star[i] - mix_mean[s[i]]) * mix_varinv[s[i]];
   }

   // Set up precision matrix
   Omega[0] = 1/sigma2;                                            // First element of main diagonal
   Omega[1] =  -phi / sigma2;                                      // First element of off-diagonal
   for(int i = 2; i < (T+1)*2; i = i+2 ){
     Omega[i] = mix_varinv[s[i/2 -1 ]] + (1 + phi*phi) / sigma2;   // main diagonal
     Omega[i+1] = -phi / sigma2;                                   // off-diagonal
   }
   Omega[2 * (T)] = mix_varinv[s[T-1]] + 1/sigma2;

   //Cholesky of Omega
   dpbtrf(Omega, T + 1);

   //solve La = c
   dtbtrs(Omega, gamma, 2);

   // add epsilon
   gamma = gamma + Rcpp::rnorm(T+1,0,1);

   //solve L'h = a + eps
   dtbtrsT(Omega, gamma, 2);

  return(gamma);

}
