# Simulate SV Process -----------------------------------------------------
#' @export
create_sv <- function(spec_grid, T_ = 1000){
  phi <- spec_grid["phi"]
  mu_h <- spec_grid["mu_h"]
  mu_y <- spec_grid["mu_y"]
  kappa <- spec_grid["kappa"]
  sigma2_u <- spec_grid["sigma2_u"]
  delta <- spec_grid["delta"]
  u            <- rnorm(T_, sd = c(0 ,rep(sqrt(sigma2_u), T_ -1)  ) )
  H_phi        <- lag_Matrix(T_, lags = -phi)
  d            <- c(mu_h, rep(mu_h * (1-phi), T_-1))
  h_true       <- as.numeric(solve(H_phi, d + u))
  epsilon      <- rnorm(T_)
  q            <- rbinom(T_, 1, kappa)
  k            <- rnorm(T_, -0.5 * delta^2, delta)
  y            <- mu_y + exp(0.5*h_true)*epsilon + q * k
  list(y = y, h = h_true, q = q, k = k, para = c("phi" = phi, "mu_h" = mu_h, "sigma2_u" = sigma2_u, "kappa" = kappa, "delta" = delta ))
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


# Function to display results compared ------------------------------------
result_plot <- function(my_results, stochvol_results){

  my_results_mu_h <- bind_rows(lapply(my_results, function(x) data.frame(mean_muh = x$results["mu_h","mean"],
                                                                         #sd_muh = x$results["mu_h","sd"],
                                                                         "qunat.05_muh" = x$results["mu_h","5%"],
                                                                         "qunat.5_muh" = x$results["mu_h","50%"],
                                                                         "qunat.95_muh" = x$results["mu_h","95%"]
  )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))

  my_results_phi <- bind_rows(lapply(my_results, function(x) data.frame(mean_phi = x$results["phi","mean"],
                                                                        #sd_phi = x$results["phi","sd"],
                                                                        "qunat.05_phi" = x$results["phi","5%"],
                                                                        "qunat.5_phi" = x$results["phi","50%"],
                                                                        "qunat.95_phi" = x$results["phi","95%"]
  )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))


  my_results_sigma2 <- bind_rows(lapply(my_results, function(x) data.frame("mean_sigma^2"  = x$results["sigma^2","mean"],
                                                                           #"sd_sigma^2" = x$results["sigma^2","sd"],
                                                                           "qunat.05_sigma^2 " = x$results["sigma^2","5%"],
                                                                           "qunat.5_sigma^2 " = x$results["sigma^2","50%"],
                                                                           "qunat.95_sigma^2 " = x$results["sigma^2","95%"]
  )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))

  stochvol_results_mu_h <- bind_rows(lapply(stochvol_results,
                                            function(x) data.frame(mean_muh = x["mu","mean"],
                                                                   #sd_muh = x["mu","sd"],
                                                                   "qunat.05_muh" = x["mu","5%"],
                                                                   "qunat.5_muh" = x["mu","50%"],
                                                                   "qunat.95_muh" = x["mu","95%"]
                                            )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))

  stochvol_results_phi <- bind_rows(lapply(stochvol_results, function(x) data.frame(mean_muh = x["phi","mean"],
                                                                                    #sd_muh = x["phi","sd"],
                                                                                    "qunat.05_muh" = x["phi","5%"],
                                                                                    "qunat.5_muh" = x["phi","50%"],
                                                                                    "qunat.95_muh" = x["phi","95%"]
  )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))

  stochvol_results_sigma2 <- bind_rows(lapply(stochvol_results, function(x) data.frame(mean_muh = x["sigma^2","mean"],
                                                                                       #sd_muh = x["sigma^2","sd"],
                                                                                       "qunat.05_sigma^2" = x["sigma^2","5%"],
                                                                                       "qunat.5_sigma^2" = x["sigma^2","50%"],
                                                                                       "qunat.95_sigma^2" = x["sigma^2","95%"]
  )
  )
  ) %>%
    gather(key = "parameter") %>%
    separate(parameter, c("summary", "parameter"), sep = "_") %>%
    group_by(summary) %>%
    mutate(spec = 1:(n()/length(unique(summary))))

  library(cowplot)
  plot_grid(nrow = 3,
            ## Phi
            ggplot() + geom_line(aes(x = spec, y = value, group = spec), col = alpha("steelblue", 1) , size = 5, data = my_results_phi %>% dplyr::filter(!summary == "sd") ) +
              geom_point(aes(x = spec, y = value, group = spec), col = alpha("steelblue", 1) , size = 5, data = my_results_phi %>% dplyr::filter(!summary == "sd") ) +
              geom_line(aes(x = spec, y = value, group = spec), col = alpha("snow4", 1),  size = 2, data = stochvol_results_phi %>% dplyr::filter(!summary == "sd") ) +
              geom_point(aes(x = spec, y = phi, group = spec), col = "red" ,data = spec_grid) + ggtitle("phi"),
            ## mu_h
            ggplot() + geom_line(aes(x = spec, y = value, group = spec), col = alpha("steelblue", 1) , size = 5,data = my_results_mu_h) +
              geom_line(aes(x = spec, y = value, group = spec), col = alpha("snow4", 0.5),  size = 2, data = stochvol_results_mu_h %>% dplyr::filter(!summary == "sd") ) +
              geom_point(aes(x = spec, y = mu_h, group = spec), col = "red" ,data = spec_grid) + ggtitle("mh"),
            ## sigma2
            ggplot() +
              geom_line(aes(x = spec, y = value, group = spec), col = alpha("steelblue", 1) , size = 5,data = my_results_sigma2 %>% dplyr::filter(!summary == "sd") ) +
              geom_point(aes(x = spec, y = value, group = spec), col = alpha("steelblue", 1) , size = 5,data = my_results_sigma2 %>% dplyr::filter(!summary == "sd") ) +
              geom_line(aes(x = spec, y = value, group = spec), col = alpha("snow4", 1),  size = 2, data = stochvol_results_sigma2 %>% dplyr::filter(!summary == "sd") ) +
              geom_point(aes(x = spec, y = sigma2_u, group = spec), col = "red" ,data = spec_grid) + ggtitle("sigma2")
  )
}
NULL
