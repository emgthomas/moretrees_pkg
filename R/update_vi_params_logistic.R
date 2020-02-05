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
  pred_g <- Matrix::Matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] * matrix(X[[g]], nrow = n) %*% mu[[g]]
  }
  # Update Sigma_g and tau_t_g
  for (g in 1:G) {
    Sigma_inv[[g]] <- 2 * Matrix::t(X[[g]]) %*% A_eta %*% X[[g]] + 
      Matrix::Diagonal(n = K[g], x = 1 / tau)
    Sigma[[g]] <- Matrix::solve(Sigma_inv[[g]])
    Sigma_det[g] <- Matrix::det(Sigma[[g]])
    tau_t[g] <- tau
    # update mu_g
    mu[[g]] <- Sigma[[g]] %*%
      Matrix::crossprod(X[[g]],
                y / 2 - 2 * g_eta * (Wdelta + apply(pred_g[, -g, drop = F], 1, sum)))
    # update prob_g (pi_g in manuscript)
    u <- 0.5 * Matrix::t(mu[[g]]) %*% Sigma_inv[[g]] %*% mu[[g]] +
      0.5 * log(Sigma_det[g]) + log(rho / (1 - rho)) - 0.5 * K[g] * log(tau_t[g])
    prob[g] <- expit(u[1, 1])
    # update pred_g
    pred_g[, g] <- prob[g] * X[[g]] %*% mu[[g]]
  }
  # Update non-sparse coefficients ---------------------------------------------------
  # Update Omega only if hyperparameters were updated at last step
  Omega_inv <- 2 * Matrix::t(W) %*% A_eta %*% W + Matrix::Diagonal(m, 1 / omega)
  if (m != 0) {
    Omega <- solve(Omega_inv)
  }
  if (m == 1) {
    Omega_det <- Omega[1, 1]
  } else {
    Omega_det <- Matrix::det(Omega)
  }
  # Update delta
  delta <- Omega %*% Matrix::t(W) %*% (y / 2 - 2 * g_eta * apply(pred_g, 1, sum))
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}