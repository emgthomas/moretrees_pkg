# --------------------------------------------------------------------------------- #
# ---------------------- Computes the ELBO for Poisson outcome  ------------------- #
# --------------------------------------------------------------------------------- #

elbo_poisson <- function(par, # variational params
                         X, y, n, K, G, sum_log_y_fac, # data
                         tau, rho, tau_alpha # hyperparameters
                         ) {
  
  prob <- 1 / (1 + exp(-par[1:G]))
  mu <- matrix(par[(G + 1):(G + G * K)], nrow = G)
  Sigma <- exp(matrix(par[(G + G * K + 1):(G + G * K + G * K)], nrow = G))
  tau_t <- exp(par[G + G * K + G * K + 1])
  mu_alpha <- par[G + G * K + G * K + 2]
  tau_t_alpha <- exp(par[G + G * K + G * K + 3])
  
  # Computing quantities needed for ELBO and hyperparameter updates ---------------
  
  # Expected sum(Y_i * log(lambda_i))
  expected_y_log_lambda <- 0
  for (g in 1:G) {
    expected_y_log_lambda <- expected_y_log_lambda + 
      prob[g] * t(mu[g, ]) %*% t(X[g, , ]) %*% y
  }
  
  # Expected sum of lambda_i
  expected_sum_lambda <- 0
  for(i in 1:n) {
    log_expected_lambda_i <- 0
    for(g in 1:G) {
      A <- 1 - prob[g] + 
        prob[g] * exp(t(mu[g, ]) %*% X[g, i, ] + t(Sigma[g, ]) %*% X[g, i, ]^2 / 2)
      log_expected_lambda_i <- log_expected_lambda_i + log(A)
    }
    expected_lambda_i <-  exp(log_expected_lambda_i)
    expected_sum_lambda <- expected_sum_lambda + expected_lambda_i
  }
  expected_sum_lambda <- expected_sum_lambda * exp(mu_alpha + tau_t_alpha / 2)
  
  # Expected sum of squared gammas
  expected_ss_gamma <- 0
  for (g in 1:G) {
    expected_ss_gamma <- expected_ss_gamma + prob[g] *
      (sum(Sigma[g, ]) + sum(mu[g, ] ^ 2))
  }
  expected_ss_gamma <- expected_ss_gamma + K * tau_t * (G - sum(prob))
  
  # Compute ELBO -------------------------------------------------------------------
  # See pg XX of manuscript; line numbers correspond to lines in equation
  line1 <- mu_alpha * sum(y) + expected_y_log_lambda - sum_log_y_fac
  line2 <- - expected_sum_lambda
  line3 <- - (expected_ss_gamma) / (2 * tau) - K * G * log(2 * pi * tau) / 2
  line4 <- log(rho ^ sum(prob)) + log((1 - rho) ^ (G - sum(prob)))
  line5 <- ((tau_t_alpha + mu_alpha ^ 2) / tau_alpha +log(2 * pi * tau_alpha)) / 2
  line6 <- (K * sum(prob) * (1 + log(2 * pi)) + sum(prob * rowSums(log(Sigma)))) / 2
  line7 <- K * (G - sum(prob)) * (1 + log(2 * pi * tau_t)) / 2
  line8 <- -1 * (sum(prob[prob != 0] * log(prob[prob != 0])) +
                  sum((1 - prob[prob != 1]) * log(1 - prob[prob != 1])))
  line9 <- (1 + log(2 * pi * tau_t_alpha)) / 2
  ELBO <- line1 + line2 + line3 + line4 + line5 + line6 + line7 + line8 + line9 
  
  # Return -------------------------------------------------------------------------
  return(list(ELBO = ELBO))
}