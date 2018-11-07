// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "sampler.h"
#include "utility.h"
using namespace Rcpp;

//' @export
//[[Rcpp::export]]
Rcpp::List  draw_hX_cpp(const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::VectorXi> s, const Eigen::Map<Eigen::MatrixXd> X,
                        const Eigen::Map<Eigen::VectorXd> beta, const double phi, const double sigma2,const double mu) {


  int T = y.size();
  Eigen::VectorXd y_star  = (y.array() * y.array()).log();

  // precision matrix
  Eigen::MatrixXd Omega    = Eigen::MatrixXd::Zero(T+1,T+1);
  Eigen::VectorXd covector = Eigen::VectorXd::Zero(T+1);


  //draw h0 from stationary distribution
  Omega(0,0) = 1 / sigma2;
  covector(0) = mu * (1-phi) / sigma2;
  // fill precision matrix / main diagonal
  for(int i = 1; i < T; ++i){
    Omega(i,i) = mix_varinv[s(i-1)] + (1 + phi * phi) / sigma2;
  }
  Omega(T,T) = mix_varinv[s(T-1)] + 1 / sigma2;
  // fill precision matrix / off-diaginal elements
  for(int i = 0; i < T - 1; ++i){
    Omega(i,i + 1) = -phi/sigma2;
    Omega(i + 1,i) = -phi/sigma2;
  }

  for(int i = 1; i < T; i++){
    covector(i) = (y_star(i-1) - mix_mean[s(i-1)]) * mix_varinv[s(i-1)] + mu * (1 - phi) * (1 - phi) * 1/sigma2;
  }
    covector(T) =  (y_star(T-1) - mix_mean[s(T-1)]) * mix_varinv[s(T-1)] + mu * (1 - phi) * 1/sigma2;

  return Rcpp::List::create(Rcpp::Named("Omega") = Omega,
                            Rcpp::Named("covector") = covector
  );
}

