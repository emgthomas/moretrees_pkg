# --------------------------------------------------------------------------------- #
# ---------------------- Computes the ELBO for Poisson outcome  ------------------- #
# --------------------------------------------------------------------------------- #

elbo_poisson <- function(X, xtx, y, n, K, G, sum_log_y_fac, # data
                                prob, mu, Sigma, Sigma_det, tau_t,
                                mu_alpha, tau_t_alpha, # variational params
                                tau, rho, mu_alpha, tau_alpha, # hyperparameters
                                update_hyper = T) { 
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  
  # Expected sum(Y_i * log(lambda_i))
  expected_y_log_lambda <- 0
  for (g in 1:G) {
    expected_y_log_lambda <- expected_y_log_lambda + 
      prob[g] * t(mu[g, ]) %*% t(X[g, , ]) %* % Y
  }
  
  # Expected sum of lambda_i
  expected_sum_lambda <- 0
  for(i in 1:n) {
    expected_sum_lambda_i <- 1
    for(g in 1:G) {
      A <- 1 - prob[g] + prob[g] * 
        exp(t(mu[g, ]) %*% X[g, i, ] + t(X[g, i, ]) %*% (Sigma[g, , ] - diag(tau_t, nrow = K)) %*% X[g, i, ] / 2)
      expected_sum_lambda_i <- expected_sum_lambda_i * A
    }
    expected_sum_lambda_i <-  exp(mu_alpha + tau_t_alpha / 2 + tau_t / 2 * xtx[i]) * expected_sum_lambda_i
    expected_sum_lambda <- expected_sum_lambda + expected_sum_lambda_i
  }
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(diag(matrix(Sigma[g, , ], nrow = K)) ) + sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- expected_ss_gamma + K * tau_t * (G - sum(prob))
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg XX of manuscript; line numbers correspond to lines in equation
  line1 <- mu_alpha * sum(Y) + expected_y_log_lambda - sum_log_y_fac
  line2_3 <- - expected_sum_lambda
  line4 <- - (expected_ss_gamma) / (2 * tau) - K * G * log(2 * pi * tau) / 2
  line5 <- log(rho ^ sum(prob)) + log((1 - rho) ^ (G - sum(prob)))
  line6 <- ((tau_t_alpha + mu_alpha ^ 2) / tau_alpha +log(2 * pi * tau_alpha)) / 2
  line7 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * log(Sigma_det))) / 2
  line8 <- K * (G - sum(prob)) * (1 + log(2 * pi * tau_t)) / 2
  line9 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line10 <- (1 + log(2 * pi * tau_t_alpha)) / 2
  ELBO <- line1 + line2_3 + line4 + line5 + line6 + line7 + line8 + line9 + line10
  
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO))
}