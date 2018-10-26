# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @import Matrix
#' @import Rcpp
#' @useDynLib svmod

### MCMC
#' @export
SVJ <- function(y, draws, post = (draws/2+1):draws, verbose = FALSE, print_every = 1000,
                init = list("mu" = -8, "phi" = 0.9, "sigma2" = 0.4,
                            "h" = -9.7, "h0" = -9,   # h could be a vector
                            "q" = 0, "zeta" = 0.05,
                            "kappa" = 0.1, "delta" = 0.1
                ),
                prior = list("b_mu" = 0, "B_mu" = 100,
                             "alpha_0" = 5, "beta_0" = 1.5,
                             "B_sigma" = 1,
                             "B_0" = cbind(c(1e12,0),c(0,1e8)),
                             "mu_delta" = 0.01, "B_delta" = 1,
                             "kappa_alpha" = 2, "kappa_beta" = 100
                )
){
  T_ <- length(y)
  # Storage
  s_draws     <- matrix(nrow = draws, ncol = T_)
  h_draws  <- matrix(nrow = draws, ncol = T_ + 1) # if h and h0 are sampled jointly
  q_draws     <- matrix(nrow = draws, ncol = T_)
  zeta_draws  <- matrix(nrow = draws, ncol = T_)
  theta_draws <- matrix(nrow = draws, ncol = 3,
                        dimnames = list(NULL, c("mu", "phi", "sigma2"))
  )
  kappa_draws <- numeric(draws)
  delta_draws <- numeric(draws)
  # Define initial values
  theta_draws[1, ] <-  c("mu"     = init[["mu"]],
                         "phi"    = init[["phi"]],
                         "sigma2" = init[["sigma2"]])
  h_draws[1, ]     <- init[["h"]]
  q_draws[1, ]     <- init[["q"]]
  zeta_draws[1, ]  <- init[["zeta"]]
  kappa_draws[1]   <- init[["kappa"]]
  delta_draws[1]   <- init[["delta"]]
  # MCMC
  for(i in 2:draws){
    if(verbose == TRUE) if(i%%print_every == 0) print(i)
    s_draws[i, ] <- draw_s_cpp(y = y - q_draws[i-1, ] * (exp(zeta_draws[i-1, ]) - 1), h = h_draws[i-1, -1])
    h_draws[i, ] <- draw_h_cpp(mu = theta_draws[i-1, "mu"],
                               phi = theta_draws[i-1, "phi"], sigma2 = theta_draws[i-1, "sigma2"],
                               y = y - q_draws[i-1, ] * (exp(zeta_draws[i-1, ]) - 1), s = s_draws[i, ])
    theta_draws[i, ] <- draw_theta_cpp(h = h_draws[i, -1], h0 =  h_draws[i, 1],
                                       mu = theta_draws[i-1, "mu"], phi = theta_draws[i-1, "phi"], sigma = sqrt(theta_draws[i-1, "sigma2"]),
                                       Bsigma = prior[["B_sigma"]], a0 = prior[["alpha_0"]], b0 = prior[["beta_0"]],
                                       bmu = prior[["b_mu"]],  Bmu= prior[["B_mu"]],
                                       B011inv= 1/prior[["B_0"]][1,1],  B022inv = 1/prior[["B_0"]][2,2])
    zeta_draws[i, ] <- draw_zeta(y = y, h = h_draws[i, -1], q = q_draws[i-1, ], delta = delta_draws[i-1])
    q_draws[i, ]    <- draw_q(y = y, h = h_draws[i, -1], zeta = zeta_draws[i, ], kappa = kappa_draws[i-1])
    kappa_draws[i]  <- draw_kappa(q = q_draws[i, ], prior = prior)
    delta_draws[i]  <- draw_delta(y = y, h = h_draws[i, -1], q = q_draws[i, ], prior = prior, delta_old = delta_draws[i-1])

  }
  invisible(list("s" = s_draws[post, ], "h" =  h_draws[post, -1], "h0" =  h_draws[post, 1], "theta" = theta_draws[post, ],
                 "zeta" = zeta_draws[post, ], "q" = q_draws[post, ], "kappa" = kappa_draws[post], "delta" = delta_draws[post]))
}



### MCMC
#' @export
SV <- function(y, draws, post = (draws/2+1):draws, verbose = FALSE, print_every = 1000,
                init = list("mu" = -8, "phi" = 0.9, "sigma2" = 0.4,
                            "h" = -9.7, "h0" = -9   # h could be a vector
                ),
                prior = list("b_mu" = 0, "B_mu" = 100,
                             "alpha_0" = 5, "beta_0" = 1.5,
                             "B_sigma" = 1,
                             "B_0" = cbind(c(1e12,0),c(0,1e8))
                )
){
  T_ <- length(y)
  # Storage
  s_draws  <- matrix(nrow = draws, ncol = T_)
  h_draws  <- matrix(nrow = draws, ncol = T_+1)
  theta_draws <- matrix(nrow = draws, ncol = 3,
                        dimnames = list(NULL, c("mu", "phi", "sigma2"))
  )
  # Define initial values
  theta_draws[1, ] <-  c("mu"     = init[["mu"]],
                         "phi"    = init[["phi"]],
                         "sigma2" = init[["sigma2"]])
  h_draws[1, ]  <- init[["h"]]
  # MCMC
  for(i in 2:draws){
    s_draws[i, ] <- draw_s_cpp(y = y, h = h_draws[i-1, ])
    h_draws[i, ]        <- draw_h_cpp(mu = theta_draws[i-1, "mu"],
                               phi = theta_draws[i-1, "phi"], sigma2 = theta_draws[i-1, "sigma2"],
                               y = y, s = s_draws[i, ])
    theta_draws[i, ] <- draw_theta_cpp(h = h_draws[i, -1], h0 =  h_draws[i, 1] ,
                                       mu = theta_draws[i-1, "mu"], phi = theta_draws[i-1, "phi"], sigma = sqrt(theta_draws[i-1, "sigma2"]),
                                       Bsigma = prior[["B_sigma"]], a0 = prior[["alpha_0"]], b0 = prior[["beta_0"]],
                                       bmu = prior[["b_mu"]],  Bmu= prior[["B_mu"]],
                                       B011inv= 1/prior[["B_0"]][1,1],  B022inv = 1/prior[["B_0"]][2,2])
    if(verbose == TRUE) if(i%%print_every == 0) print(i)
  }
  list("s" = s_draws[post, ], "h" =  h_draws[post, -1], "h0" =  h_draws[post, 1], "theta" = theta_draws[post, ] )
}



#' @export
SV_my_theta <- function(y, draws, post = (draws/2+1):draws, verbose = FALSE, print_every = 1000,
                         init = list("mu" = -8, "phi" = 0.9, "sigma2" = 0.4,
                                     "h" = -9.7, "h0" = -9   # h could be a vector
                         ),
                         prior = list("b_mu" = 0, "B_mu" = 100,
                                      "alpha_0" = 5, "beta_0" = 1.5,
                                      "B_sigma" = 1,
                                      "B_0" = cbind(c(1e12,0),c(0,1e8))
                         )
){
  T_ <- length(y)
  # Storage
  s_draws  <- matrix(nrow = draws, ncol = T_)
  h_draws  <- matrix(nrow = draws, ncol = T_ + 1) # if h and h0 are sampled jointly
  theta_draws <- matrix(nrow = draws, ncol = 3,
                        dimnames = list(NULL, c("mu", "phi", "sigma2"))
  )
  # Define initial values
  theta_draws[1, ] <-  c("mu"     = init[["mu"]],
                         "phi"    = init[["phi"]],
                         "sigma2" = init[["sigma2"]])
  h_draws[1, ]  <- init[["h"]]
  # MCMC
  for(i in 2:draws){
    s_draws[i, ]  <- draw_s_cpp(y = y, h = h_draws[i-1, ])
    h_draws[i, ]  <- draw_h_cpp(mu = theta_draws[i-1, "mu"],
                         phi = theta_draws[i-1, "phi"], sigma2 = theta_draws[i-1, "sigma2"],
                         y = y, s = s_draws[i, ])
    theta_draws[i, ] <- draw_theta(h = h_draws[i, -1], h0 =  h_draws[i, 1], prior = prior, theta_old =  theta_draws[i-1, ], T_ = T_)
    if(verbose == TRUE) if(i%%print_every == 0) print(i)
  }

  list("s" = s_draws[post, ], "h" =  h_draws[post, ], "h0" =  h0_draws[post, 1], "theta" = theta_draws[post, ] )
}

#############################################################################################################
#############################################################################################################


# old Gibbs sampler functions -------------------------------------------------
draw_h0 <- function(h, mu, phi, sigma2){
  h_hat <- (mu + phi *(h[1] - mu))
  rnorm(1, mean = h_hat, sd = sigma2)
}

## old h
draw_h <- function(h0, mu, phi, sigma2, y, s, comp){
  T_     <- length(y)
  y_star <- log(y^2)
  H_phi   <- lag_Matrix(T_ = T_, lags = -phi)
  alpha   <- c(phi * h0, rep(0, T_ -1))
  c       <- rep(mu*(1-phi),T_)
  d_s     <- sapply(s, function(s)comp$mu[s]) - 1.2704
  Sigma_s <- Matrix::Diagonal(n = T_, sapply(s, function(s)comp$sig2[s]))
  P_h <- 1/sigma2 * t(H_phi)%*%H_phi

  P       <- P_h  + solve(Sigma_s)
  h_hat   <- solve(P, P_h %*% solve(H_phi,(alpha + c)) + solve(Sigma_s, y_star - d_s))
  precision_sampler(h_hat, P)
}

## old mixture sampler
draw_s <- function(y, h, comp){
  y_star <- log(y^2)
  sig_prob_raw <- matrix(ncol = 7, nrow = length(y))
  for(t in 1:length(y)){
    sig_prob_raw[t, ] <- comp$p * dnorm(y_star[t], h[t] + comp$mu - 1.2704, sqrt(comp$sig2))
  }
  sig_prob <- prop.table(sig_prob_raw, 1)
  s_draw   <- apply(sig_prob, 1, function(prob) sample(1:7, 1, prob = prob) )
  s_draw
}

