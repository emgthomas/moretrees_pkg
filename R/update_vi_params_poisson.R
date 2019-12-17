# --------------------------------------------------------------------------------- #
# ---------- Performs one step in VI optimization for Poisson outcome  ------------ #
# --------------------------------------------------------------------------------- #

update_vi_params_poisson <- function(X, y, Xy, n, K, G, # data
                                    u, mu, Sigma, mu_alpha, tau_t_alpha, expA, # variational params
                                    rho, tau ) { # hyperparams
  prob <- 1 / (1 + exp(-u))
  exp_int <- exp(mu_alpha + tau_t_alpha / 2)
  # Update mu_g, Sigma_g -----------------------------------------------------------
  b <- numeric(n)
  for (g in 1:G) {
    # compute w vector
    w <- numeric(n)
    for (i in 1:n) {
      d <- prod(1 - prob[-g] + prob[-g] * expA[i, -g])
      w[i] <- expA[i, g] * d
    }
    # Sigma_g
    Sigma[g, , ] <- solve(diag(1 / tau, nrow = K) + prob[g] * exp_int *
                            t(X[g, , ]) %*% diag(w) %*% X[g, , ])
    # mu_g
    mu[g, ] <- mu[g, ] + Sigma[g, , ] %*% 
      (prob[g] * t(X[g, , ]) %*% (y - exp_int * w) - mu[g, ] / tau)
    # Now need to update expA[ , g]
    for (i in 1:n) {
      b[i] <- t(X[g, i, ]) %*% Sigma[g, , ] %*% X[g, i, ]
    }
    expA[, g] <- exp(mu[g, ] %*% t(X[g, , ])  + b / 2)
  }
  # Update mu_alpha, tau_t_alpha ---------------------------------------------------
  b <- 0
  for (i in 1:n) {
    b <- b + prod(1 - prob + prob * expA[i, ])
  }
  tau_t_alpha <- 1 / (exp_int * b + tau_t_alpha / tau_alpha)
  mu_alpha <- mu_alpha + tau_t_alpha * (sum(y) - mu_alpha / tau_alpha - exp_int * b)
  exp_int <- exp(mu_alpha + tau_t_alpha / 2)
  # Update u_g ---------------------------------------------------------------------
  for (g in 1:G) {
    a <- 0
    for (i in 1:n) {
      d <- prod(1 - prob[-g] + prob[-g] * expA[i , -g])
      a <- a + (1 - expA[i, g]) * d
    }
    u[g] <- t(mu[g , ]) %*% Xy[g, ] + log(rho / (1 - rho)) + a * exp_int
    prob[g] <- 1 / (1 + exp(-u[g]))
  }
  # Return -------------------------------------------------------------------------
  return(list(u = u, mu = mu, Sigma = Sigma,
              mu_alpha = mu_alpha, tau_t_alpha = tau_t_alpha,
              expA = expA))
}


