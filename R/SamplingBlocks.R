#' @export
draw_q <- function(y, h, zeta, kappa){
  p1 <- kappa * dnorm(y, exp(zeta) - 1, sqrt(exp(h)))
  p0 <- (1 - kappa) * dnorm(y, 0, sqrt(exp(h)))
  q <- rbinom(length(y), 1, p1 / (p1 + p0))
  q
}

#' @export
draw_kappa <- function(q, prior){
  rbeta(1, prior[["kappa_alpha"]] + sum(q), prior[["kappa_beta"]] + length(q) - sum(q))
}

#' @export
draw_zeta <- function(y, h, q, delta){
  exph <- exp(h)
  zeta_var = 1/(1/delta^2 + 1/exph[q == 1])
  zeta_hat = zeta_var * (-.5 + (y[q == 1])/exph[q == 1])
  zeta_draw <- numeric(length = length(y))
  zeta_draw[q== 0] <- rnorm(length(y) - sum(q), -0.5 * delta^2, delta)
  zeta_draw[q== 1] <- rnorm(sum(q), zeta_hat, sqrt(zeta_var))

  #zeta_draw <- ifelse(q == 0, rnorm(length(y) - sum(q), -0.5 * delta^2, delta), rnorm(sum(q), zeta_hat, sqrt(zeta_var)))
  zeta_draw
}

#' @export
draw_delta <- function(y, h, q, prior, delta_old){
  exph <- exp(h)
  lp_delta <- function(x) -log(x) -.5/prior[["B_delta"]]*(log(x)-prior[["mu_delta"]])^2 -.5*sum(log(x^2*q^2+exph)) -.5*sum((y+.5*x^2*q)^2/(x^2*q^2+exph))
  prop_mean <-  optimise(lp_delta, interval = seq(0,1, length.out = 100000), maximum = TRUE)$maximum
  prop_var <- 0.01^2
  delta_prop <- rnorm(1, prop_mean ,sqrt(prop_var))
  if (delta_prop > 0){
    alpha_MH = lp_delta(delta_prop) - lp_delta(delta_old) + .5*(delta_old - prop_mean)^2/prop_var -.5*(delta_prop-prop_mean)^2/prop_var
  } else alpha_MH = -Inf
  if(alpha_MH > log(runif(1))){
    return(delta_prop)
  } else return(delta_old)
}


#' @export
draw_theta <- function(h, h0, prior, theta_old, T_){
  sigma     <- sqrt(theta_old["sigma2"])
  X         <- cbind(1, c(h0, h[-T_]))
  B_T       <- inv_Eigen(crossprod(X) + solve(prior[["B_0"]]))  ## B_0 prior for proposal density
  b_T       <- B_T %*% t(X) %*% h
  C_T       <- 1/2 * (sum(h^2) - t(b_T) %*% t(X) %*% h)
  c_T       <- (T_-1)/2
  ## draw proposal for sigma2
  sigma2_prop <- 1/rgamma(1, c_T, C_T)
  ## draw proposal for gamma and phi
  beta_prop <- b_T[2:1] + sigma * chol_Eigen_R(B_T[2:1, 2:1]) %*% rnorm(2)
  phi_prop  <- beta_prop[1]
  if(phi_prop > 1 | phi_prop < -1){  #to make sure proposal lies inside the unit circle
    return(theta_old)
  }
  gamma_prop  <- beta_prop[2]

  ## MH
  ### Posterior
  logp_h0_theta_new <- lognorm(h0, gamma_prop / (1-phi_prop), sigma2_prop/(1-phi_prop^2))
  logp_h0_theta_old <- lognorm(h0, theta_old[["mu"]], theta_old[["sigma2"]]/(1-theta_old[["phi"]]^2))

  logp_gamma_phi_new <- lognorm(gamma_prop, prior[["b_mu"]] * (1-phi_prop), sqrt(prior[["B_mu"]]) * (1-phi_prop) )
  logp_gamma_phi_old <- lognorm(gamma_prop, prior[["b_mu"]] * (1-theta_old[["phi"]]), sqrt(prior[["B_mu"]]) * (1-theta_old[["phi"]]) )

  logp_phi_new       <- logdbeta((phi_prop + 1)/2, prior[["alpha_0"]], prior[["beta_0"]])
  logp_phi_old       <- logdbeta((theta_old[["phi"]] + 1)/2, prior[["alpha_0"]], prior[["beta_0"]])

  ### Proposal
  logpprop_gamma_new <-  lognorm(gamma_prop, 0, prior[["B_0"]][1,1] )
  logpprop_gamma_old <-  lognorm(theta_old[["mu"]] * (1-theta_old[["phi"]] ), 0, prior[["B_0"]][1,1] )

  logpprop_phi_new <- lognorm(phi_prop, 0, prior[["B_0"]][2,2] )
  logpprop_phi_old <- lognorm(theta_old[["phi"]], 0, prior[["B_0"]][2,2] )

  ## All that has to do with sigma2
  logacceptratiosigma2 <- logacceptrateGamma(sigma2_prop, theta_old[["sigma2"]], prior[["B_sigma"]])



  alpha <- logacceptratiosigma2 + logp_h0_theta_new + logp_gamma_phi_new + logp_phi_new -
    logp_h0_theta_old - logp_gamma_phi_old - logp_phi_old +
    logpprop_gamma_old + logpprop_phi_old - logpprop_gamma_new - logpprop_phi_new

  if(log(runif(1)) < alpha){
    return(c("mu" = gamma_prop / (1-phi_prop),
             "phi" = phi_prop,
             "sigma2" = sigma2_prop))
  } else(
    return(theta_old))

}


#' @export
draw_thetaX <- function(h, h0, prior, theta_old, T_, X){
  sigma     <- sqrt(theta_old["sigma2"])
  X         <- cbind(1, c(h0, h[-T_]), X)
  B_T       <- inv_Eigen(crossprod(X) + solve(prior[["B_0"]]))  ## B_0 prior for proposal density
  b_T       <- B_T %*% t(X) %*% h
  C_T       <- 1/2 * (sum(h^2) - t(b_T) %*% t(X) %*% h)
  c_T       <- (T_-1)/2
  ## draw proposal for sigma2
  sigma2_prop <- 1/rgamma(1, c_T, C_T)
  ## draw proposal for gamma and phi
  beta_prop <- b_T + sigma * chol_Eigen(B_T) %*% rnorm(3)
  phi_prop  <- beta_prop[2]
  if(phi_prop > 1 | phi_prop < -1){  #to make sure proposal lies inside the unit circle
    return(theta_old)
  }
  gamma_prop  <- beta_prop[1]
  beta_prop   <- beta_prop[3]
  ## MH
  ### Posterior
  logp_h0_theta_new <- lognorm(h0, gamma_prop / (1-phi_prop), sigma2_prop/(1-phi_prop^2))
  logp_h0_theta_old <- lognorm(h0, theta_old[["mu"]], theta_old[["sigma2"]]/(1-theta_old[["phi"]]^2))

  logp_gamma_phi_new <- lognorm(gamma_prop, prior[["b_mu"]] * (1-phi_prop), sqrt(prior[["B_mu"]]) * (1-phi_prop) )
  logp_gamma_phi_old <- lognorm(gamma_prop, prior[["b_mu"]] * (1-theta_old[["phi"]]), sqrt(prior[["B_mu"]]) * (1-theta_old[["phi"]]) )

  logp_phi_new       <- logdbeta((phi_prop + 1)/2, prior[["alpha_0"]], prior[["beta_0"]])
  logp_phi_old       <- logdbeta((theta_old[["phi"]] + 1)/2, prior[["alpha_0"]], prior[["beta_0"]])

  ### Proposal
  logpprop_gamma_new <-  lognorm(gamma_prop, 0, prior[["B_0"]][1,1] )
  logpprop_gamma_old <-  lognorm(theta_old[["mu"]] * (1-theta_old[["phi"]] ), 0, prior[["B_0"]][1,1] )

  logpprop_phi_new <- lognorm(phi_prop, 0, prior[["B_0"]][2,2] )
  logpprop_phi_old <- lognorm(theta_old[["phi"]], 0, prior[["B_0"]][2,2] )

  ## All that has to do with sigma2
  lopacceptratiosigma2 <- logacceptrateGamma(sigma2_prop, theta_old[["sigma2"]], prior[["B_sigma"]])

  alpha <- lopacceptratiosigma2 + logp_h0_theta_new + logp_gamma_phi_new + logp_phi_new -
    logp_h0_theta_old - logp_gamma_phi_old - logp_phi_old +
    logpprop_gamma_old + logpprop_phi_old - logpprop_gamma_new - logpprop_phi_new

  if(log(runif(1)) < alpha){
    return(c("mu" = gamma_prop / (1-phi_prop),
             "phi" = phi_prop,
             "sigma2" = sigma2_prop,
             "beta" = beta_prop))
  } else(
    return(theta_old))

}

