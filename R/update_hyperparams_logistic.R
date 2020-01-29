# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

update_hyperparams_logistic <- function(X, W, y, n, K, G, m, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                delta, Omega, Omega_det, 
                                eta, g_eta, # variational params
                                omega, tau, rho, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Expected linear predictor
  lp <- W %*% delta
  for (g in 1:G) {
    lp <- lp + prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K) 
  }
  # Expected linear predictor squared
  lp2 <- apply(W, MARGIN = 1, 
               FUN = function(w, Omega) t(w) %*% Omega %*% w, Omega = Omega)
  for (g in 1:G) {
    lp2 <- lp2 + prob[g] * apply(X[g, , ], MARGIN = 1, 
          FUN = function(x, Sigma) t(x) %*% Sigma %*% x, Sigma = Sigma[[g]])
  }
  lp2 <- lp2 + lp ^ 2

  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(diag(Sigma[[g]])) + sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + K * sum(tau_t * (1 - prob)))
  
  # Expected sum of squared thetas
  expected_ss_theta <- sum(diag(Omega)) + sum(delta ^ 2)
  
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    omega <- expected_ss_theta / m
    tau <- expected_ss_gamma / (K * G)
    rho <- mean(prob)
  }
  # Compute ELBO -------------------------------------------------------------------
  # See pg 13 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- (1 / 2) * t(y) %*% lp + 
    sum(log(expit(eta))) - sum(eta) / 2 + g_eta %*% eta ^ 2
  line2 <- - g_eta %*% lp2
  line3 <- - expected_ss_gamma / (2 * tau) - 
    K * G * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (G - sum(prob)))
  line4 <- expected_ss_theta / (2 * omega) - 
    (m / 2) * log(2 * pi * omega)
  line5 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line6 <- (K / 2) * (G - sum(prob)) + 
    (K / 2) * sum(log(2 * pi * tau_t) * (1 - prob))
  line7 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line8 <- (m + log(Omega_det) + m * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8
  # Update eta
  eta <- as.numeric(sqrt(lp2))
  g_eta <- as.numeric(gfun(eta))
  # Return -------------------------------------------------------------------------
  return(list(ELBO = as.numeric(ELBO), omega = omega, tau = tau, rho = rho,
              eta = eta, g_eta = g_eta))
}