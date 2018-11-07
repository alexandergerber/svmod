#include "utility.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::MatrixXd;
// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(const Rcpp::NumericVector & omega_diag, double omega_offdiag, double * chol_diag, double * chol_offdiag)
{
  chol_diag[0] = sqrt(omega_diag[0]);  // maybe speed up via iterators?
  for (int j = 1; j < omega_diag.length(); j++) {
    chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
    chol_diag[j] = sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
  }
}

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag, const Rcpp::NumericVector & covector, double * htmp)
{
  htmp[0] = covector[0]/chol_diag[0];
  for (int j = 1; j < chol_diag.length(); j++) {
    htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
  }
}

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag, const Rcpp::NumericVector & htmp, double * h)
{
  int T = chol_diag.length();
  h[T-1] = htmp[T-1]/chol_diag[T-1];
  for (int j = T-2; j >= 0; j--) {
    h[j] = (htmp[j] - chol_offdiag[j]*h[j+1])/chol_diag[j];
  }
}

double logdnorm(double x, double mu, double sigma) {
  return -log(sigma)-((x-mu)*(x-mu)/(2*sigma*sigma));
}

// non-normalized log-density for Beta(a, b)
double logdbeta(double x, double a, double b) {
  return (a-1)*log(x)+(b-1)*log(1-x);
}

double logacceptrateGamma(double xnew, double xold, double Bsigma) {
  return (xold-xnew)/(2*Bsigma);
}

//' @export
//[[Rcpp::export]]
double invGamme_cpp(double cT, double CT){
  return(1/Rcpp::as<double>(Rcpp::rgamma(1, cT, 1/CT)));
}
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd inv_Eigen(Eigen::Map<Eigen::MatrixXd> M){
  return M.inverse();
}
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd chol_Eigen(Eigen::Map<Eigen::MatrixXd> M){
  Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
  Eigen::MatrixXd L = lltOfA.matrixL();
  return L;
}

//' @export
//[[Rcpp::export]]
SEXP eigen_mult1(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> X){
  return(Rcpp::wrap(A * X));
}

//' @export
//[[Rcpp::export]]
Eigen::MatrixXd eigen_mult2(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> X){
  return(A * X);
}






