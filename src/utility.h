#pragma once
#include <RcppArmadillo.h>

// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholTridiag(const Rcpp::NumericVector & omega_diag, double omega_offdiag,
                 double * chol_diag, double * chol_offdiag);

// Solves Chol*x = covector ("forward algorithm")
void forwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag,
                const Rcpp::NumericVector & covector, double * htmp);

// Solves (Chol')*x = htmp ("backward algorithm")
void backwardAlg(const Rcpp::NumericVector & chol_diag, const Rcpp::NumericVector & chol_offdiag,
                 const Rcpp::NumericVector & htmp, double * h);


double logdnorm(double x, double mu, double sigma);
double logdbeta(double x, double a, double b);
double logacceptrateGamma(double xnew, double xold, double Bsigma);
double invGamme_cpp(double cT, double CT);

