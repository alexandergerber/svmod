// #include "auxmix.h"
// #include <lapacke.h>
// #include <iostream>
//
// //' @export
// //[[Rcpp::export]]
// Rcpp::NumericVector test_LAPACK(Rcpp::IntegerVector s, Rcpp::NumericVector y, double phi, double sigma2, double mu){
//   int T = s.length();
//   Rcpp::NumericVector y_star = log(y * y );
//   Rcpp::NumericVector Omega((T+1)*2);
//   Rcpp::NumericVector c(T+1);
//
//   Omega[0] = 1/sigma2;                                            // First element of main diagonal
//   Omega[1] =  -phi / sigma2;                                      // First element of off-diagonal
//   for(int i = 2; i < (T+1)*2; i = i+2 ){
//     c[i/2] = mix_varinv[s[i/2-1]]*( y_star[i/2-1] - mix_mean[s[i/2-1]]) + mu * (1 - phi) * (1-phi) / sigma2;
//     Omega[i] = mix_varinv[s[i/2 -1 ]] + (1 + phi*phi) / sigma2;     // main diagonal
//     Omega[i+1] = -phi / sigma2;                                   // off-diagonal
//   }
//   Omega[2 * (T)] = mix_varinv[s[T-1]] + 1/sigma2;
//
//
//   c[0]   = mu * (1 - phi)/sigma2;
//   // for(int ii = 1; ii < T; ii++){
//   //   c[ii] = mix_varinv[s[ii-1]]*( y_star[ii-1] - mix_mean[s[ii-1]]) + mu * (1 - phi) * (1-phi) / sigma2;                                      // off-diagonal
//   // }
//   c[T] = (y_star[T-1] - mix_mean[s[T-1]])*mix_varinv[s[T-1]] + mu*(1-phi)/sigma2;
//
//   // for debuggig
//   //Rcpp::NumericVector co1 = clone(c);
//
//   // Rcpp::Rcout << "Vector c after init" << std::endl;
//   // Rcpp::Rcout << c << std::endl;
//
//   // //LAPACK
//   //    lapack_int info1, info2, info3;
//   //    lapack_int N = T+1;
//   //    lapack_int KD = 1;
//   //    lapack_int LDAB = KD + 1;
//   //    info1 = LAPACKE_dpbtrf( LAPACK_COL_MAJOR, 'L', N, KD, Omega.begin(), LDAB);  // computes Cholesky Factor of Omega and replaces Omega with that
//   //    info2 = LAPACKE_dtbtrs( LAPACK_COL_MAJOR, 'L', 'N', 'N', N, KD, 1, Omega.begin(), LDAB, c.begin(), N);
//   //
//   //    // Rcpp::Rcout << "Vector a" << std::endl;
//   //    // Rcpp::Rcout << c << std::endl;
//   //
//   //   c = c + Rcpp::rnorm(T+1,0,1);
//   //
//   //   // Rcpp::Rcout << "Vector a + eps" << std::endl;
//   //   // Rcpp::Rcout << c << std::endl;
//   //
//   //   info3 = LAPACKE_dtbtrs( LAPACK_COL_MAJOR, 'L', 'T', 'N', N, KD, 1, Omega.begin(), LDAB, c.begin(), N);
//
//   // Rcpp::Rcout << "Vector h" << std::endl;
//   // Rcpp::Rcout << c << std::endl;
//
//   return(c);
//
// }
