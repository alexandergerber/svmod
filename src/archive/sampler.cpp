// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadilloExtensions/sample.h>
// #include "sampler.h"
// #include "utility.h"
// using namespace Rcpp;
//
// //' @export
// //[[Rcpp::export]]
// Rcpp::IntegerVector draw_s_cpp(Rcpp::NumericVector y, Rcpp::NumericVector h) {
//   NumericMatrix mixprob(10, y.length());
//   NumericVector y_star = log(y * y );
//   findMixCDF(&mixprob(0,0), y_star - h);
//
//   IntegerVector r(y.length());  // mixture indicators
//   invTransformSampling(&mixprob(0,0), &r(0),  y.length());
//   return  r;
// }
//
// //' @export
// //[[Rcpp::export]]
// Rcpp::NumericVector draw_h_cpp(Rcpp::NumericVector y, Rcpp::IntegerVector s, double phi, double sigma2, double mu) {
//   int T = y.length();
//   const double sigma2inv = pow(sigma2, -1);
//
//
//   NumericVector y_star = log(y * y );
//   NumericVector omega_diag(T+1);  // contains diagonal elements of precision matrix
//   double omega_offdiag;  // contains off-diag element of precision matrix (const)
//   NumericVector chol_offdiag(T), chol_diag(T+1);  // Cholesky-factor of Omega
//   NumericVector covector(T+1);  // holds covector (see McCausland et al. 2011)
//   NumericVector htmp(T+1);  // intermediate vector for sampling h
//   NumericVector hnew(T+1);  // intermediate vector for sampling h
//
//   omega_diag[0] = sigma2inv;
//   covector[0]   = mu * (1 - phi) * sigma2inv;
//
//   for (int j = 1; j < T; j++) {
//     omega_diag[j] = mix_varinv[s[j-1]] + (1+phi*phi)*sigma2inv;
//     covector[j] = (y_star[j-1] - mix_mean[s[j-1]])*mix_varinv[s[j-1]]
//     + mu*(1-phi)*(1-phi)*sigma2inv;
//   }
//   omega_diag[T] = mix_varinv[s[T-1]] + sigma2inv;
//   covector[T] = (y_star[T-1] - mix_mean[s[T-1]])*mix_varinv[s[T-1]]
//   + mu*(1-phi)*sigma2inv;
//   omega_offdiag = -phi*sigma2inv;
//
//
//   // for debuggig
//   // Rcpp::Rcout << "Vector c after init" << std::endl;
//   // Rcpp::Rcout << covector << std::endl;
//
//
//   // Cholesky decomposition
//   cholTridiag(omega_diag, omega_offdiag, &chol_diag(0), &chol_offdiag(0));
//
//   // Solution of Chol*x = covector ("forward algorithm")
//   forwardAlg(chol_diag, chol_offdiag, covector, &htmp(0));
//
//   // Rcpp::Rcout << "Vector a" << std::endl;
//   // Rcpp::Rcout << htmp << std::endl;
//
//   htmp = htmp + Rcpp::rnorm(T+1,0,1);
//   // Rcpp::Rcout << "Vector a + eps  " << std::endl;
//   // Rcpp::Rcout << htmp << std::endl;
//   // Solution of (Chol')*x = htmp ("backward algorithm")
//   backwardAlg(chol_diag, chol_offdiag, htmp, &hnew(0));
//
//   // Rcpp::Rcout << "Vector h" << std::endl;
//   // Rcpp::Rcout << hnew << std::endl;
//
//   return hnew;
// }
//
//
// // Step (b): sample mu, phi, sigma - __CENTERED__ version:
// //' @export
// //[[Rcpp::export]]
// Rcpp::NumericVector draw_theta_cpp(
//     double h0, const Rcpp::NumericVector &h,
//     double mu, double phi, double sigma,
//     double Bsigma,
//     double a0, double b0,
//     double bmu, double Bmu,
//     double B011inv, double B022inv) {
//
//   double Bh0inv = 1-phi*phi;
//
//   int T = h.length();
//   double z, CT, cT, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, tmp1,
//     BT11, BT12, BT22, bT1 = 0, bT2 = 0, chol11, chol12, chol22, phi_prop,
//     gamma_prop, tmpR, tmpR2, logR;
//   double R = -10000.;
//   double sigma2_prop = -10000.;
//   double sigma_prop = -10000.;
//   Rcpp::NumericVector innov(2);
//   Rcpp::NumericVector quant(2);
//
//   // first calculate bT and BT:
//   sum1 = h[0];
//   sum3 = h0*h[0];
//   sum4 = h[0]*h[0];
//   for (int j = 1; j < T-1; j++) {
//     sum1 += h[j];
//     sum3 += h[j-1]*h[j];
//     sum4 += h[j]*h[j];
//   }
//   sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
//   sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
//   sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
//   sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2
//
//   tmp1 = 1/(((sum4 + B011inv)*(T+B022inv)-sum1*sum1));
//   BT11 = (T + B022inv)*tmp1;
//   BT12 = -sum1*tmp1;
//   BT22 = (sum4+B011inv)*tmp1;
//
//   bT1 = BT11*sum3 + BT12*sum2;
//   bT2 = BT12*sum3 + BT22*sum2;
//
//    // draw sigma^2
//   cT = (T-1)/2.0; // we want IG(-.5,0) as proposal
//   CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
//   sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
//
//   // Draw beta
//   chol11 = sqrt(BT11);
//   chol12 = (BT12/chol11);
//   chol22 = sqrt(BT22-chol12*chol12);
//   chol11 *= sigma;
//   chol12 *= sigma;
//   chol22 *= sigma;
//
//   innov = Rcpp::rnorm(1);
//   phi_prop = bT1 + chol11*innov[0];
//   if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
//     Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
//     return ret;
//   }
//   // Draw gamma
//   else gamma_prop = bT2 + chol12*innov[0] + chol22*Rcpp::rnorm(1)[0];
//
//
//   // acceptance probability exp(R) calculated on a log scale
//   tmpR = 1-phi_prop;  // some temps used for acceptance probability
//   tmpR2 = 1-phi;
//
//   R = 0.;  // initialize R
//   // accept sigma jointly with "betas"
//   sigma_prop = sqrt(sigma2_prop);
//   R  = logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma);  // initialize R
//   R += logdnorm(h0, gamma_prop/tmpR, sigma_prop/sqrt(1-phi_prop*phi_prop));
//   R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
//   R += logdnorm(gamma_prop, bmu*tmpR, sqrt(Bmu)*tmpR);
//   R -= logdnorm(mu*tmpR2, bmu*tmpR2, sqrt(Bmu)*tmpR2);
//   R += logdbeta((phi_prop+1)/2, a0, b0);
//   R -= logdbeta((phi+1)/2, a0, b0);
//   R += logdnorm(phi, 0, sigma/sqrt(B011inv));
//   R -= logdnorm(phi_prop, 0, sigma_prop/sqrt(B011inv));
//   R += logdnorm(mu*tmpR2, 0, sigma/sqrt(B011inv));
//   R -= logdnorm(gamma_prop, 0, sigma_prop/sqrt(B011inv));
//
//     // accept/reject
//     if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
//       mu = gamma_prop/(1-phi_prop);
//       phi = phi_prop;
//       sigma = sigma_prop;
//     }
//
//   Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma * sigma);
//   return ret;
// }
//
// // Step (b): sample mu, phi, sigma - __CENTERED__ version:
// //' @export
// //[[Rcpp::export]]
// Rcpp::NumericVector draw_thetaX_cpp(
//     double h0, const Rcpp::NumericVector &h,
//     double mu, double phi, double sigma, double beta,
//     double Bsigma,
//     double a0, double b0,
//     double bmu, double Bmu,
//     double B011inv, double B022inv, double B033inv) {
//
//   double Bh0inv = 1-phi*phi;
//
//   int T = h.length();
//   double z, CT, cT, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, tmp1,
//     BT11, BT12, BT22, bT1 = 0, bT2 = 0, chol11, chol12, chol22, phi_prop,
//     gamma_prop, tmpR, tmpR2, logR;
//   double R = -10000.;
//   double sigma2_prop = -10000.;
//   double sigma_prop = -10000.;
//   Rcpp::NumericVector innov(2);
//   Rcpp::NumericVector quant(2);
//
//   // first calculate bT and BT:
//   sum1 = h[0];
//   sum3 = h0*h[0];
//   sum4 = h[0]*h[0];
//   for (int j = 1; j < T-1; j++) {
//     sum1 += h[j];
//     sum3 += h[j-1]*h[j];
//     sum4 += h[j]*h[j];
//   }
//   sum2 = sum1 + h[T-1];  // h_1 + h_2 + ... + h_T
//   sum1 += h0;            // h_0 + h_1 + ... + h_{T-1}
//   sum3 += h[T-2]*h[T-1]; // h_0*h_1 + h_1*h_2 + ... + h_{T-1}*h_T
//   sum4 += h0*h0;         // h_0^2 + h_1^2 + ... + h_{T-1}^2
//
//   tmp1 = 1/(((sum4 + B011inv)*(T+B022inv)-sum1*sum1));
//   BT11 = (T + B022inv)*tmp1;
//   BT12 = -sum1*tmp1;
//   BT22 = (sum4+B011inv)*tmp1;
//
//   bT1 = BT11*sum3 + BT12*sum2;
//   bT2 = BT12*sum3 + BT22*sum2;
//
//   // draw sigma^2
//   cT = (T-1)/2.0; // we want IG(-.5,0) as proposal
//   CT = .5*((sum4 - h0*h0 + h[T-1]*h[T-1]) - bT1*sum3 - bT2*sum2);
//   sigma2_prop = 1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT));
//
//   // Draw beta
//   chol11 = sqrt(BT11);
//   chol12 = (BT12/chol11);
//   chol22 = sqrt(BT22-chol12*chol12);
//   chol11 *= sigma;
//   chol12 *= sigma;
//   chol22 *= sigma;
//
//   innov = Rcpp::rnorm(1);
//   phi_prop = bT1 + chol11*innov[0];
//   if ((phi_prop >= 1) || (phi_prop <= -1)) { // outside the unit sphere
//     Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma);
//     return ret;
//   }
//   // Draw gamma
//   else gamma_prop = bT2 + chol12*innov[0] + chol22*Rcpp::rnorm(1)[0];
//
//
//   // acceptance probability exp(R) calculated on a log scale
//   tmpR = 1-phi_prop;  // some temps used for acceptance probability
//   tmpR2 = 1-phi;
//
//   R = 0.;  // initialize R
//   // accept sigma jointly with "betas"
//   sigma_prop = sqrt(sigma2_prop);
//   R  = logacceptrateGamma(sigma2_prop, sigma*sigma, Bsigma);  // initialize R
//   R += logdnorm(h0, gamma_prop/tmpR, sigma_prop/sqrt(1-phi_prop*phi_prop));
//   R -= logdnorm(h0, mu, sigma/sqrt(1-phi*phi));
//   R += logdnorm(gamma_prop, bmu*tmpR, sqrt(Bmu)*tmpR);
//   R -= logdnorm(mu*tmpR2, bmu*tmpR2, sqrt(Bmu)*tmpR2);
//   R += logdbeta((phi_prop+1)/2, a0, b0);
//   R -= logdbeta((phi+1)/2, a0, b0);
//   R += logdnorm(phi, 0, sigma/sqrt(B011inv));
//   R -= logdnorm(phi_prop, 0, sigma_prop/sqrt(B011inv));
//   R += logdnorm(mu*tmpR2, 0, sigma/sqrt(B011inv));
//   R -= logdnorm(gamma_prop, 0, sigma_prop/sqrt(B011inv));
//
//   // accept/reject
//   if (log(Rcpp::as<double>(Rcpp::runif(1))) < R) {
//     mu = gamma_prop/(1-phi_prop);
//     phi = phi_prop;
//     sigma = sigma_prop;
//   }
//
//   Rcpp::NumericVector ret = Rcpp::NumericVector::create(mu, phi, sigma * sigma);
//   return ret;
// }
