#include "utility_Eigen.h"

//' @export
//[[Rcpp::export]]
Eigen::MatrixXd chol_Eigen(Eigen::Map<Eigen::MatrixXd> M){
  Eigen::LLT<Eigen::MatrixXd> lltOfA(M);
  Eigen::MatrixXd L = lltOfA.matrixL();
  return L;
}

//' @export
//[[Rcpp::export]]
Eigen::MatrixXd inv_Eigen(Eigen::Map<Eigen::MatrixXd> M){
  return M.inverse();
}
