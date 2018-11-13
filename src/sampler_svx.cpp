// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include "sampler.h"
#include "utility.h"
using namespace Rcpp;

//' @export
//[[Rcpp::export]]
Rcpp::List  draw_hX_cpp(const Eigen::Map<Eigen::VectorXd> y, Rcpp::IntegerVector s, const Eigen::Map<Eigen::MatrixXd> X,
                        const Eigen::Map<Eigen::VectorXd> beta, const double phi, const double sigma2, const double mu) {

  int T = y.size();

  Eigen::VectorXd XX(T+1);
  XX << 0, X;

  Eigen::VectorXd y_star  = (y.array() * y.array()).log();

  // Stuff from the  h side

  Eigen::VectorXd gamma(T+1);
  for(int i = 0; i < T+1; i++){
    if(i == 0){
      gamma(i) = mu;
    } else{
      gamma(i) = (1-phi) * mu;
    }
  }

  Eigen::SparseMatrix<double> H = lag_matrix(T+1, phi);

  Rcpp::NumericVector P_u_vec(T+1);
  for(int i = 0; i < T+1; i++){
    if(i == 0){
      P_u_vec(i) = (1-phi * phi)/sigma2;
    } else{
      P_u_vec(i) = 1/sigma2;
    }
  }

  Eigen::SparseMatrix<double> P_u = diagonal(P_u_vec);
  Eigen::SparseMatrix<double> P_h = H.transpose() * P_u * H;

  //stuff from the y side

  Eigen::VectorXd c_y(T+1);
  for(int i = 0; i < T+1; i++){
    if(i == 0){
      c_y(i) = 0;
    } else{
      c_y(i) = y_star(i-1) - mix_mean[s(i-1)];
    }
  }

  Rcpp::NumericVector P_y_vec(T+1);
  for(int i = 0; i < T+1; i++){
    if(i == 0){
      P_y_vec(i) = 0;
    } else{
      P_y_vec(i) = mix_varinv[s(i-1)];
    }
  }
  Eigen::SparseMatrix<double> P_y = diagonal(P_y_vec);

  // Overall precision Matrix
  Eigen::SparseMatrix<double> P = P_y + P_h;


  Eigen::VectorXd c = P_y * c_y + P_h * solve_tri(H, gamma +  XX * beta);

  return Rcpp::List::create(Rcpp::Named("lag_matrix") = H,
                            Rcpp::Named("P_u") = P_u,
                            Rcpp::Named("P_h") = P_h,
                            Rcpp::Named("P") = P,
                            Rcpp::Named("c") = c
  );
}


