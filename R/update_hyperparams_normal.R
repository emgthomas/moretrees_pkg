# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

update_hyperparams_normal <- function(X, XtX, W, WtW, y, n, K, G, m, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                delta, Omega, Omega_det, # variational params
                                omega, sigma2, tau, rho, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  # Sum of squared residuals
  lp <- W %*% delta
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
  expected_ssr <- as.numeric(trace_prod(Omega, WtW) + prob_tr + ssr + ssr_corr)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(diag(matrix(Sigma[g, , ], nrow = K))) + sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- as.numeric(expected_ss_gamma + K * sum(tau_t * (1 - prob)))
  
  # Expected sum of squared thetas
  expected_ss_theta <- sum(diag(Omega)) + sum(delta ^ 2)
  
  # Update hyperparameters ---------------------------------------------------------
  if (update_hyper) {
    omega <- expected_ss_theta / m
    sigma2 <- expected_ssr / n
    tau <- expected_ss_gamma / (K * G)
    rho <- mean(prob)
  }
  # Compute ELBO -------------------------------------------------------------------
  # See pg 5 of "variational inference for spike & slab model" document -
  # line numbers correspond to lines in equation
  line1 <- - 1 / (2 * sigma2) * expected_ssr - (n / 2) * log(2 * pi * sigma2)
  line2 <- - expected_ss_gamma / (2 * tau) - 
    K * G * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) +
    log((1 - rho) ^ (G - sum(prob)))
  line3 <- expected_ss_theta / (2 * omega) - 
    (m / 2) * log(2 * pi * omega)
  line4 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line5 <- (K / 2) * (G - sum(prob)) + 
    (K / 2) * sum(log(2 * pi * tau_t) * (1 - prob))
  line6 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line7 <- (m + log(Omega_det) + m * log(2 * pi)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO, omega = omega, sigma2 = sigma2, tau = tau, rho = rho))
}