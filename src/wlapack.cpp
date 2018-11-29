#include <Rcpp.h>
#include <R_ext/Lapack.h>
//' @export
//[[Rcpp::export]]
void dgemv(Rcpp::NumericMatrix A, Rcpp::NumericVector x, Rcpp::NumericVector y){
//computes A * x + y and modifies y
  //new on
  int M = A.nrow();
  int N = A.ncol();
  double ALPHA = 1;
  int LDA =  A.nrow();
  double BETA = 1;
  int INCX = 1;
  int INCY = 1;
  //Rcpp::NumericVector Y(M);
  F77_CALL(dgemv)( "N", &M, &N, &ALPHA, A.begin(), &LDA, x.begin(), &INCX, &BETA, y.begin(), &INCY);
  }

//' @export
//[[Rcpp::export]]
void dpbtrf(Rcpp::NumericVector A, int N){ // Cholesky of banded matrix
  int info;
  int KD = 1;
  int LDAB = KD + 1;
  int NRHS = 1;
  F77_CALL(dpbtrf)( "L", &N, &KD, A.begin(), &LDAB, &info);
}
//' @export
//[[Rcpp::export]]
void dtbtrs(Rcpp::NumericVector A, Rcpp::NumericVector x, int ND){
  int info;
  int N = A.size() / ND;
  int KD = 1;
  int LDAB = KD + 1;
  int NRHS = 1;
  F77_CALL(dtbtrs)("L", "N", "N", &N, &KD, &NRHS, A.begin(), &LDAB, x.begin(), &N, &info);
}

//' @export
//[[Rcpp::export]]
void dtbtrsT(Rcpp::NumericVector A, Rcpp::NumericVector x, int ND){
  int info;
  int N = A.size() / ND;
  int KD = 1;
  int LDAB = KD + 1;
  int NRHS = 1;
  F77_CALL(dtbtrs)("L", "T", "N", &N, &KD, &NRHS, A.begin(), &LDAB, x.begin(), &N, &info);
}

//' @export
//[[Rcpp::export]]
void dsbmv(Rcpp::NumericVector A, Rcpp::NumericVector x, int N){
  int info;
  int K = 1;
  double ALPHA = 1;
  double BETA = 1;
  int LDA = 2;
  int INCX = 1;
  int INCY = 1;
  Rcpp::NumericVector y(N);
  F77_CALL(dsbmv)("L", &N, &K, &ALPHA, A.begin(), &LDA, x.begin(),  &INCX, &BETA, y.begin(), &INCY);
}


//' @export
//[[Rcpp::export]]

void UBMV(Rcpp::NumericVector& v, double phi, double sigma2){ //compute H'Sig_u^-1 (X \beta + gamma)
  int T = v.size();
  v[0] = v[0] * (1 - phi * phi) / sigma2 + v[1] * -phi/sigma2;
  for(int i = 1; i < T-1; i++){
    v[i] = v[i] * 1/sigma2 + v[i+1] * -phi/sigma2;
  }
  v[T-1] = v[T-1] / sigma2;
}






