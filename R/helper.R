# Simulate SV Process -----------------------------------------------------
#' @export
create_sv <- function(spec_grid, X, T_ = 1000){
  phi <- spec_grid[["phi"]]
  mu_h <- spec_grid[["mu_h"]]
  mu_y <- spec_grid[["mu_y"]]
  kappa <- spec_grid[["kappa"]]
  sigma2_u <- spec_grid[["sigma2_u"]]
  delta <- spec_grid[["delta"]]
  beta <- spec_grid[["beta"]]
  u            <- rnorm(T_, sd = c(0 ,rep(sqrt(sigma2_u), T_ -1)  ) )
  H_phi        <- lag_Matrix(T_, lags = -phi)
  d            <- c(mu_h, rep(mu_h * (1-phi), T_-1))
  h_true       <- as.numeric(solve(H_phi, X %*% beta + d + u))
  epsilon      <- rnorm(T_)
  if(kappa !=0){
  q            <- rbinom(T_, 1, kappa)
  k            <- rnorm(T_, -0.5 * delta^2, delta)
  } else{
    q <- 0
    k <- 0
  }
  y            <- mu_y + exp(0.5*h_true)*epsilon + q * k
  list(y = y, h = h_true, q = q, X = X, k = k, para = c(phi = phi, beta = beta, mu_h = mu_h, sigma2_u = sigma2_u, kappa = kappa, delta = delta))
}

# Draw Normal using Precision Matrix --------------------------------------
precision_sampler <- function(mu, P){
  B_prime <- chol(P)
  Z       <- rnorm(length(mu))
  U       <- mu + solve(B_prime, Z)
  U[ ,1]
}


# Lower triangular band batrix ---------------------------------------------
#' @export
lag_Matrix <- function(T_, lags = numeric()){
  lags  <- c(1, lags)
  diag  <- lapply(lags, function(x) rep(x, T_)  )
  bandSparse(n = T_, k = 0:-(length(lags)-1), diag = diag)
}

inv_22 <- function(X) {
  inv_det <- 1/(X[1,1] * X[2,2] - X[2,1] * X[1,2])
  X2 <- rbind(c(X[2,2], - X[1,2]), c(X[1,1], -X[2,1]) )
  inv <- inv_det * X2
  inv
}

# Log kernels for some distributions --------------------------------------
lognorm <- function(x, mu, sd){
  -log(sd)-((x-mu)*(x-mu)/(2*sd*sd))
}

logdbeta <- function(x, alpha, beta){
  (alpha-1)*log(x)+(beta-1)*log(1-x)
}

loggamma <- function(x, alpha, beta){
  (alpha - 1) * log(x) - beta * x
}

logacceptrateGamma <- function(xnew, xold, Bsigma){
  (xold-xnew)/(2*Bsigma)
}


