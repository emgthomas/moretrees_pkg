# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

update_vi_params_normal <- function(X, XtX, y, n, K, G, # data
                                 prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, # variational params
                                 sigma2, rho, tau, # hyperparams
                                 update_hyper_last) { 
  # Compute linear predictor for each unit (needed for updates below) ------------------
  pred_g <- matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
  }
  const <- log(rho / (1 - rho)) - 0.5 * K * log(tau_t)
  for (g in 1:G) {
    # Update tau_t ---------------------------------------------------------------------
    # only changes if tau has been updated; see below
    # Update Sigma ---------------------------------------------------------------------
    if (update_hyper_last) {
      # Sigma only needs to be updated if hyperparams were updated last step
      Sigma_inv[g, , ] <- XtX[g, , ] + diag(1 / tau, K)
      Sigma[g, , ] <- solve(Sigma_inv[g, , ])
      if (K == 1) {
        Sigma_det[g] <- Sigma[g, 1, 1]
      } else {
        Sigma_det[g] <- det(Sigma[g, , ])
      }
    }
    # Update mu ------------------------------------------------------------------------
    mu[g, ] <- (1 / sigma2) * matrix(Sigma[g, , ], nrow = K) %*%
      crossprod(matrix(X[g, , ], nrow = n), y - rowSums(pred_g[, -g, drop = F]))
    # Update prob (pi in manuscript) ---------------------------------------------------
    u <- crossprod(matrix(mu[g, ], nrow = K), matrix(Sigma_inv[g, , ], nrow = K)) %*%
      matrix(mu[g, ], nrow = K) + 0.5 * log(Sigma_det[g]) + const
    prob[g] <- exp(loglogit(u[1, 1]))
    # Update pred_g --------------------------------------------------------------------
    if (g != G) {
      pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
    }
  }

  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, sigma2 = sigma2, 
              rho = rho, tau = tau))
}