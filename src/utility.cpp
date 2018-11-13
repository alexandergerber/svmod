// [[Rcpp::depends(RcppEigen)]]
#include "utility.h"

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
Eigen::MatrixXd eigen_mult(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> X){
  return(A * X);
}


Eigen::SparseMatrix<double> lag_matrix(int n, double phi){
  typedef Eigen::Triplet<double> T;
  // Create H_phi
  std::vector<T> tripletList_H;
  tripletList_H.reserve(n + n -1);
  for(int i; i < n; i++)
  {
    tripletList_H.push_back(T(i,i,1));
    if(i < n-1){
      tripletList_H.push_back(T(i+1, i, -phi));
    }
  }
  Eigen::SparseMatrix<double> H(n,n);
  H.setFromTriplets(tripletList_H.begin(), tripletList_H.end());
  return(H);
}

//' @export
//[[Rcpp::export]]
Eigen::SparseMatrix<double> diagonal(Rcpp::NumericVector x){
  typedef Eigen::Triplet<double> T;
  // Create H_phi
  std::vector<T> tripletList;
  tripletList.reserve(x.size());
  for(int i; i < x.size(); i++)
  {
    tripletList.push_back(T(i,i,x(i)));
  }
  Eigen::SparseMatrix<double> D(x.size(),x.size());
  D.setFromTriplets(tripletList.begin(), tripletList.end());
  return(D);
}

//' @export
//[[Rcpp::export]]
Eigen::VectorXd  solve_tri (Eigen::SparseMatrix<double> A, Eigen::VectorXd b){
  Eigen::VectorXd x = A.triangularView<Eigen::Lower>().solve(b);
  return(x);
}

//' @export
//[[Rcpp::export]]
Rcpp::List chol_sparse (Eigen::Map<Eigen::SparseMatrix<double> > A){
  Eigen::SimplicialLLT <Eigen::SparseMatrix<double> > cholesky;
   cholesky.analyzePattern(A);
   cholesky.factorize(A);
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("L") = cholesky.matrixL());
}







