# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

#'   \code{update_hyperparams_logistic_moretrees} Performs hyperparameter updates and computes 
#'   current value of ELBO in VI algorithm for bernoulli outcomes.

update_hyperparams_logistic_moretrees <- function(X, W, y, 
                                        outcomes_units,
                                        ancestors,
                                        n, K, p, m, # data
                                        prob, mu, Sigma, Sigma_det, tau_t,
                                        delta, Omega, Omega_det, 
                                        eta, g_eta, # variational params
                                        omega, tau, rho, # hyperparameters
                                        model = "ss",
                                        update_hyper = T) {
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Expected linear predictor
  xi_u <- mapply(FUN = function(prob, mu) prob * mu,
                 prob = prob, mu = mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    beta_v <- Reduce(`+`, xi_u[ancestors[[v]]])
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
    lp2[outcomes_units[[v]]] <- apply(W[outcomes_units[[v]], , drop = F], MARGIN = 1, 
                  FUN = function(w, Omega) crossprod(w, Omega) %*% w,
                  Omega = Omega_v) +
        apply(X[outcomes_units[[v]], , drop = F], 1,
                 FUN = function(x, Sigma) crossprod(x, Sigma) %*% x,
                 Sigma = Sigma_v)
  }
  lp2 <- lp2 + lp ^ 2
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (v in 1:p) {
    expected_ss_gamma <- expected_ss_gamma + prob[v] *
      (sum(diag(Sigma[[v]])) + sum(mu[[v]] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + sum(K * tau_t * (1 - prob)))
  
  # Expected sum of squared thetas
  if (m == 0) {
    expected_ss_theta <- 0
  } else {
    expected_ss_theta <- sum(mapply(FUN = 
           function(Omega, delta) sum(diag(Omega)) + sum(delta ^ 2),
           Omega = Omega, delta = delta))
  }
  
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    if (m != 0) {
      omega <- expected_ss_theta / m
    }
    tau <- expected_ss_gamma / sum(K)
    rho <- mean(prob)
  }
  # Update eta  --------------------------------------------------------------------
  eta <- as.numeric(sqrt(lp2))
  g_eta <- as.numeric(gfun(eta))
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg 13 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- (1 / 2) * t(y) %*% lp + 
    sum(logexpit(eta)) - sum(eta) / 2 + g_eta %*% (eta ^ 2)
  line2 <- - g_eta %*% lp2
  line3 <- - expected_ss_gamma / (2 * tau) - 
    sum(K) * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (p - sum(prob)))
  line4 <- - expected_ss_theta / (2 * omega) - 
    (m / 2) * log(2 * pi * omega)
  line5 <- (sum(K * prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line6 <- (1 / 2) * sum(K * (1 - prob)) + 
    (1 / 2) * sum(K * log(2 * pi * tau_t) * (1 - prob))
  line7 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line8 <- sum(m + log(Omega_det) + m * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8
  # Return -------------------------------------------------------------------------
  return(list(ELBO = as.numeric(ELBO), omega = omega, tau = tau, rho = rho,
              eta = eta, g_eta = g_eta))
}