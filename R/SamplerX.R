#' @import Matrix
#' @import Rcpp
#' @useDynLib svmod


#' @export
SVJX <- function(y, X, draws, post = (draws/2+1):draws, verbose = FALSE, print_every = 1000, keep_all = FALSE,
                init = list("mu" = -8, "phi" = 0.9, "sigma2" = 0.4, "beta" = 0,
                            "h" = -9.7, "h0" = -9,   # h could be a vector
                            "q" = 0, "zeta" = 0.05,
                            "kappa" = 0.1, "delta" = 0.1
                ),
                prior = list("b_mu" = 0, "B_mu" = 100,
                             "alpha_0" = 5, "beta_0" = 1.5,
                             "B_sigma" = 1,
                             "B_0" = cbind(c(1e12,0, 0), c(0,1e8, 0), c(0,0,1e8)),
                             "mu_delta" = 0.01, "B_delta" = 1,
                             "kappa_alpha" = 2, "kappa_beta" = 100
                )
){
  T_ <- length(y)
  # Storage
  theta_draws <- matrix(nrow = draws, ncol = 4,
                        dimnames = list(NULL, c("mu", "phi", "sigma2", "beta"))
  )
  kappa_draws <- numeric(draws)
  delta_draws <- numeric(draws)
  # Define initial values
  theta_draws[1, ] <-  c("mu"     = init[["mu"]],
                         "phi"    = init[["phi"]],
                         "sigma2" = init[["sigma2"]],
                         "beta" = init[["beta"]])

  kappa_draws[1]   <- init[["kappa"]]
  delta_draws[1]   <- init[["delta"]]

# storage of smaller objects
  if(keep_all == TRUE){
    #set.seed(1)
    s_draws     <- matrix(nrow = draws, ncol = T_)
    h_draws  <- matrix(nrow = draws, ncol = T_ + 1) # if h and h0 are sampled jointly
    q_draws     <- matrix(nrow = draws, ncol = T_)
    zeta_draws  <- matrix(nrow = draws, ncol = T_)

    h_draws[1, ]     <- init[["h"]]
    q_draws[1, ]     <- init[["q"]]
    zeta_draws[1, ]  <- init[["zeta"]]


  # MCMC
  for(i in 2:draws){
   # set.seed(1)
    if(verbose == TRUE) if(i%%print_every == 0) print(i)
    s_draws[i, ] <- draw_s(y = y - q_draws[i-1, ] * (exp(zeta_draws[i-1, ]) - 1), h = h_draws[i-1, -1])
    h_draws[i, ] <- draw_hX(mu = theta_draws[i-1, "mu"], beta = as.vector(theta_draws[i-1, "beta"]),
                               phi = theta_draws[i-1, "phi"], sigma2 = theta_draws[i-1, "sigma2"],
                               y = y - q_draws[i-1, ] * (exp(zeta_draws[i-1, ]) - 1), s = s_draws[i, ], X = as.matrix(X))
    theta_draws[i, ] <- draw_thetaX(h = h_draws[i, -1], h0 =  h_draws[i, 1],
                                    prior = prior, theta_old = theta_draws[i-1, ], T_, X = X[-1, ])
    zeta_draws[i, ] <- draw_zeta(y = y, h = h_draws[i, -1], q = q_draws[i-1, ], delta = delta_draws[i-1])
    q_draws[i, ]    <- draw_q(y = y, h = h_draws[i, -1], zeta = zeta_draws[i, ], kappa = kappa_draws[i-1])
    kappa_draws[i]  <- draw_kappa(q = q_draws[i, ], prior = prior)
    delta_draws[i]  <- draw_delta(y = y, h = h_draws[i, -1], q = q_draws[i, ], prior = prior, delta_old = delta_draws[i-1])

  }
  return(invisible(list("s" = s_draws[post, ], "h" =  h_draws[post, -1], "h0" =  h_draws[post, 1], "theta" = theta_draws[post, ],
                 "zeta" = zeta_draws[post, ], "q" = q_draws[post, ], "kappa" = kappa_draws[post], "delta" = delta_draws[post]))
  )
  }

  # storage of larger objects
  if(keep_all == FALSE){

    s_draws     <- matrix(nrow = 1, ncol = T_)
    h_draws     <- matrix(nrow = 1, ncol = T_ + 1) # if h and h0 are sampled jointly
    q_draws     <- matrix(nrow = 1, ncol = T_)
    zeta_draws  <- matrix(nrow = 1, ncol = T_)

    h_draws[1, ]     <- init[["h"]]
    q_draws[1, ]     <- init[["q"]]
    zeta_draws[1, ]  <- init[["zeta"]]

    # MCMC
    for(i in 2:draws){
      if(verbose == TRUE) if(i%%print_every == 0) print(i)
      s_draws[1, ] <- draw_s(y = y - q_draws[1, ] * (exp(zeta_draws[1, ]) - 1), h = h_draws[1, -1])
      h_draws[1, ] <- draw_hX(mu = theta_draws[i-1, "mu"], beta = as.vector(theta_draws[i-1, "beta"]),
                              phi = theta_draws[i-1, "phi"], sigma2 = theta_draws[i-1, "sigma2"],
                              y = y - q_draws[1, ] * (exp(zeta_draws[1, ]) - 1), s = s_draws[1, ], X = as.matrix(X))
      theta_draws[i, ] <- draw_thetaX(h = h_draws[1, -1], h0 =  h_draws[1, 1],
                                      prior = prior, theta_old = theta_draws[i-1, ], T_, X = X[-1, ])
      zeta_draws[1, ] <- draw_zeta(y = y, h = h_draws[1, -1], q = q_draws[1, ], delta = delta_draws[i-1])
      q_draws[1, ]    <- draw_q(y = y, h = h_draws[1, -1], zeta = zeta_draws[1, ], kappa = kappa_draws[i-1])
      kappa_draws[i]  <- draw_kappa(q = q_draws[1, ], prior = prior)
      delta_draws[i]  <- draw_delta(y = y, h = h_draws[1, -1], q = q_draws[1, ], prior = prior, delta_old = delta_draws[i-1])
    }
    return(invisible(list("s" = s_draws, "h" =  h_draws[-1], "h0" =  h_draws[1], "theta" = theta_draws[post, ],
                          "zeta" = zeta_draws, "q" = q_draws, "kappa" = kappa_draws[post], "delta" = delta_draws[post]))
    )
  }


}



