#include "auxmix.h"

//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector draw_s(Rcpp::NumericVector y, Rcpp::NumericVector h) {
  Rcpp::NumericMatrix mixprob(10, y.length());
  Rcpp::NumericVector y_star = log(y * y );
  findMixCDF(&mixprob(0,0), y_star - h);

    Rcpp::IntegerVector r(y.length());  // mixture indicators
 invTransformSampling(&mixprob(0,0), &r(0),  y.length());
  return  r;
}
