# --------------------------------------------------------------------------------- #
# ---------------------- Computes the ELBO for Poisson outcome  ------------------- #
# --------------------------------------------------------------------------------- #

elbo_poisson <- function(X, y, n, K, G, sum_log_y_fac, # data
                         u, mu, Sigma, mu_alpha, tau_t_alpha, expA, # variational params
                         rho, tau, tau_alpha) { # hyperparams
  
  prob <- 1 / (1 + exp(-u))
  
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  
  # Expected sum(Y_i * log(lambda_i))
  expected_y_log_lambda <- mu_alpha * sum(y)
  for (g in 1:G) {
    expected_y_log_lambda <- expected_y_log_lambda + 
      prob[g] * t(mu[g, ]) %*% t(X[g, , ]) %*% y
  }

  # Expected sum of lambda_i
  expected_sum_lambda <- 0
  for(i in 1:n) {
    expected_sum_lambda <- expected_sum_lambda + prod(1 - prob + prob * expA[i, ])
  }
  expected_sum_lambda <- expected_sum_lambda * exp(mu_alpha + tau_t_alpha / 2)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + sum(mu[g, ] ^ 2) + sum(diag(Sigma[g, , ]))
  }
  expected_ss_gamma <- expected_ss_gamma / (2 * tau)
  
  # Determinant of Sigma
  if(K == 1) {
    Sigma_det <- 1 / Sigma[1:G, , ]
  } else {
    Sigma_det <- plyr::aaply(.data = Sigma, .margins = 1, .fun = det)
  }
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg XX of manuscript; line numbers correspond to lines in equation
  line1 <- expected_y_log_lambda - sum_log_y_fac
  line2 <- - expected_sum_lambda
  line3 <- - expected_ss_gamma - K * G * log(2 * pi * tau) /2
  line4 <- log(rho ^ sum(prob)) + log((1 - rho) ^ (G - sum(prob)))
  line5 <- - ((mu_alpha ^ 2 + tau_t_alpha) / tau_alpha + log(2 * pi* tau_alpha)) / 2
  line6 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                   sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line7 <-  (G * K + sum(log(Sigma_det)) + G * K * log(2 * pi) + 1 + log(2 * pi* tau_t_alpha)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7
  
  # Return -------------------------------------------------------------------------
  return(as.numeric(ELBO))
}