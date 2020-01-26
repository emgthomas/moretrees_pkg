# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

update_vi_params_logistic <- function(X, W, y, n, K, G, m, # data
                                    prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                    delta, Omega, Omega_inv, Omega_det, 
                                    eta, g_eta, # variational params
                                    omega, rho, tau) { # hyperparams
  # Update sparse coefficients ------------------------------------------------------
  Wdelta <- W %*% delta
  A_eta <- Matrix::Diagonal(n = n, g_eta)
  pred_g <- matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
  }
  # Update Sigma_g and tau_t_g
  for (g in 1:G) {
    Sigma_inv[[g]] <- 2 * t(X[g, , ]) %*% A_eta %*% X[g, , ] + diag(1 / tau, K)
    Sigma[[g]] <- solve(Sigma_inv[[g]])
    if (K == 1) {
      Sigma_det[g] <- Sigma[[g]]
    } else {
      Sigma_det[g] <- det(Sigma[[g]])
    }
    tau_t[g] <- tau
    # update mu_g
    mu[g, ] <- Sigma[[g]] %*%
      crossprod(matrix(X[g, , ], nrow = n),
                y / 2 - 2 * g_eta * (Wdelta + rowSums(pred_g[, -g, drop = F])))
    # update prob_g (pi_g in manuscript)
    u <- crossprod(matrix(mu[g, ], nrow = K), matrix(Sigma_inv[[g]],nrow = K)) %*%
      matrix(mu[g, ], nrow = K) + 0.5 * log(Sigma_det[g]) + 
      log(rho / (1 - rho)) - 0.5 * K * log(tau_t[g])
    prob[g] <- expit(as.numeric(u))
    # update pred_g
    pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K) 
  }
  # Update non-sparse coefficients ---------------------------------------------------
  # Update Omega only if hyperparameters were updated at last step
  Omega_inv <- 2 * t(W) %*% A_eta %*% W + diag(1 / omega, m)
  Omega <- solve(Omega_inv)
  if (m == 1) {
    Omega_det <- Omega[1, 1]
  } else {
    Omega_det <- det(Omega)
  }
  # Update delta
  delta <- as.numeric(Omega %*% t(W) %*% (y / 2 - 2 * g_eta * rowSums(pred_g)))
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}