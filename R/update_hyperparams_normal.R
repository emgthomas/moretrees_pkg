# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

update_hyperparams_normal <- function(X, XtX, y, n, K, G, # data
                                prob, mu, Sigma, Sigma_det, tau_t, # variational params
                                sigma2, tau, rho, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Sum of squared residuals
  lp <- numeric(n) + 0
  for (g in 1:G) {
    lp <- lp + prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
  }
  ssr <- sum( (y - lp) ^ 2 )
  # Expected sum of squared residuals
  prob_tr <- 0
  ssr_corr <- 0
  for (g in 1:G) {
    prob_tr <- prob_tr + 
      prob[g] * trace_prod(Sigma[g, , ], XtX[g, , ])
    ssr_corr <- ssr_corr + 
      prob[g] * (1 - prob[g]) * t(mu[g, ]) %*% XtX[g, , ] %*% mu[g, ]
  }
  expected_ssr <- as.numeric(prob_tr + ssr + ssr_corr)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(diag(matrix(Sigma[g, , ], nrow = K)) ) + sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- expected_ss_gamma + K * tau_t * (G - sum(prob))
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    sigma2 <- expected_ssr / n
    tau <- expected_ss_gamma / (K * G)
    rho <- mean(prob)
  }
  # Compute ELBO -------------------------------------------------------------------
  # See pg XX of manuscript; line numbers correspond to lines in equation
  line1 <- -1 / (2 * sigma2) * expected_ssr - (n / 2) * log(2 * pi * sigma2)
  line2 <- - expected_ss_gamma / (2 * tau) - 
    K * G * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (G - sum(prob)))
  line3 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line4 <- K * (G - sum(prob)) * (1 + log(2 * pi * tau_t)) / 2
  line5 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  ELBO <- line1 + line2 + line3 + line4 + line5
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO, sigma2 = sigma2, tau = tau, rho = rho))
}