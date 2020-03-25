# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic_moretrees} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for bernoulli outcomes.

update_hyperparams_logistic_moretrees <- function(X, W, y, 
                                                  outcomes_units,
                                                  ancestors,
                                                  levels,
                                                  n, K, p, m, Fg, # dsgn
                                                  prob, mu, Sigma, Sigma_det, tau_t,
                                                  delta, Omega, Omega_det, 
                                                  a_t_rho, b_t_rho, # vi_params
                                                  eta, g_eta, 
                                                  tau, omega, # hyperparams
                                                  a_rho, b_rho, # hyper_fixed
                                                  update_hyper) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Expected linear predictor
  xi <- mapply(FUN = function(prob, mu) prob * mu,
               prob = prob, mu = mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    beta_v <- Reduce(`+`, xi[ancestors[[v]]])
    theta_v <- Reduce(`+`, delta[ancestors[[v]]])
    lp[outcomes_units[[v]]] <- X[outcomes_units[[v]], ] %*% beta_v + 
      W[outcomes_units[[v]], ] %*% theta_v
  }
  # Expected linear predictor squared
  Sigma_u <- mapply(FUN = function(prob, Sigma, mu) prob * 
                      (Sigma + (1 - prob) * tcrossprod(mu)),
                    prob = prob, Sigma = Sigma, mu = mu,
                    SIMPLIFY = F)
  lp2 <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    Sigma_v <- Reduce(`+`, Sigma_u[ancestors[[v]]])
    Omega_v <- Reduce(`+`, Omega[ancestors[[v]]])
    lp2[outcomes_units[[v]]] <- quadFormByRow(Omega_v, W[outcomes_units[[v]], , drop = F]) +
      quadFormByRow(Sigma_v, X[outcomes_units[[v]], , drop = F])
  }
  lp2 <- lp2 + lp ^ 2
  
  # Expected sum of squared gammas
  expected_ss_gamma <- numeric(p)
  for (v in 1:p) {
    expected_ss_gamma[v] <- prob[v] * (sum(diag(Sigma[[v]])) + sum(mu[[v]] ^ 2)) + 
      (1 - prob[v]) * K * tau_t[v]
  }
  
  # Expected sum of squared thetas
  expected_ss_theta <- numeric(p)
  for (v in 1:p) {
    expected_ss_theta[v] <- (sum(diag(Omega[[v]])) + sum(delta[[v]] ^ 2))
  }
  
  # Expected log(rho) and log(1-rho)
  expected_l_rho <- digamma(a_t_rho) - digamma(a_t_rho + b_t_rho)
  expected_l_1m_rho <- digamma(b_t_rho) - digamma(a_t_rho + b_t_rho)
  
  # Update eta  --------------------------------------------------------------------
  eta <- as.numeric(sqrt(lp2))
  g_eta <- as.numeric(gfun(eta))
  
  # Update hyperparams  ------------------------------------------------------------
  if (update_hyper == T) {
    for (f in 1:Fg) {
      tau[f] <- sum(expected_ss_gamma[levels == f]) / (K * sum(levels == f))
    }
    if (m > 0) {
      for (f in 1:Fg) {
        omega[f] <- sum(expected_ss_theta[levels == f]) / (m * sum(levels == f))
      }
    }
  }
  
  # Compute ELBO -------------------------------------------------------------------
  # See "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- (1 / 2) * crossprod(y, lp) + sum(logexpit(eta)) - sum(eta) / 2 + g_eta %*% (eta ^ 2)
  line2 <- - g_eta %*% lp2
  line3 <- - sum(expected_ss_gamma / (2 * tau[levels])) - K * sum(log(2 * pi * tau[levels])) / 2
  line4 <- sum(expected_l_rho[levels] * prob + expected_l_1m_rho[levels] * (1 - prob))
  line5 <- - sum(expected_ss_theta / (2 * omega[levels])) - m * sum(log(2 * pi * omega[levels])) / 2
  line6 <- sum((a_rho - 1) * expected_l_rho + (b_rho - 1) * expected_l_1m_rho - mapply(lbeta, a_rho, b_rho))
  line7 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line8 <- (K / 2) * sum(1 - prob) + (K / 2) * sum(log(2 * pi * tau_t) * (1 - prob))
  line9 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                    sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line10 <- (m * p + sum(log(Omega_det)) + m * p * log(2 * pi)) / 2
  line11 <- -1 * sum((a_t_rho - 1) * expected_l_rho + (b_t_rho - 1) * expected_l_1m_rho - mapply(lbeta, a_t_rho, b_t_rho))
  
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + line9 + line10 + line11
  
  # Return -------------------------------------------------------------------------
  return(list(ELBO = as.numeric(ELBO), 
              eta = eta, g_eta = g_eta,
              tau = tau, omega = omega))
}