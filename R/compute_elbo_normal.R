# --------------------------------------------------------------------------------- #
# --------------------- Computes the ELBO for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

compute_elbo_normal <- function(X, XtX, y, n, K, G, # data
                                prob, mu, Sigma, Sigma_det, tau_t, # variational parameters
                                sigma2, rho, tau, # hyperparameters
                                update.hyper = F) {
  # Computing quantities needed for ELBO --------------------------------------------
  # Sum of squared residuals
  lp <- numeric(n) + 0
  for (g in 1:G) {
    lp <- lp + prob[g] * matrix(X[g, , ], nrow=n) %*% matrix(mu[g, ], nrow=K)
  }
  ssr <- sum((y - lp) ^ 2)
  # Expected sum of squared residuals
  prob_tr <- 0
  for (g in 1:G) {
    prob_tr <- prob_tr + prob[g] * trace_prod(Sigma[g, , ], XtX[g, , ])
  }
  expected_ssr <- prob_tr + ssr
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] * (sum(diag(matrix(Sigma[g, , ], nrow=K))) + 
                                                          sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- expected_ss_gamma + K * tau_t * (G-sum(prob))
  # Compute ELBO -------------------------------------------------------------------
  # See pg XX of manuscript; line numbers correspond to lines in equation
  line1 <- -1 / (2 * sigma2) * expected_ssr - (n/2) * log(2 * pi * sigma2)
  line2 <- -expected_ss_gamma / (2 * tau) - 
    K * G * log(2 * pi * tau) / 2 + 
    log(rho ^ sum(prob)) + 
    log((1 - rho) ^ (G - sum(prob)))
  line3 <- 0.5 * (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det)))
  line4 <- (G - sum(prob)) * K * 0.5 * (1 + log(2 * pi * tau_t))
  line5 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) + 
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  ELBO <- line1 + line2 + line3 + line4 + line5
  # Update hyperparameters ---------------------------------------------------------
  if (update.hyper) {
    sigma2 <- expected_ssr / n
    tau <- expected_ss_gamma / (K * G)
    rho <- mean(prob)
    return(list(ELBO = ELBO, sigma2 = sigma2, tau = tau, rho = rho))
  } else {
    return(ELBO)
  }
}