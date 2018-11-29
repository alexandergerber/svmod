#pragma once
// Ax + y
void dgemv(Rcpp::NumericMatrix A, Rcpp::NumericVector x, Rcpp::NumericVector y);
// Cholesky Decomposition for symmetrix band matrix
void dpbtrf(Rcpp::NumericVector A, int N);
// Solve Ax = b for triangular A
void dtbtrs(Rcpp::NumericVector A, Rcpp::NumericVector x, int ND);
void dsbmv(Rcpp::NumericVector A, Rcpp::NumericVector x, int N);
void UBMV(Rcpp::NumericVector& v, double phi, double sigma2);
void dtbtrsT(Rcpp::NumericVector A, Rcpp::NumericVector x, int ND);
