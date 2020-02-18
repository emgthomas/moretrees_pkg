# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

#' \code{update_vi_logistic_moretrees} Performs variational updates for bernoulli outcomes.

update_vi_params_logistic_moretrees <- function(X, W, y,
                                      outcomes_nodes,
                                      n, K, p, m, 
                                      prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, 
                                      delta, Omega, Omega_inv, Omega_det, 
                                      eta, g_eta, # variational params
                                      omega, rho, tau) { # hyperparams
  # Update sparse coefficients ------------------------------------------------------
  pred_v <- Matrix::Matrix(0, nrow = n, ncol = p)
  xi_u <- mapply(`*`, prob, mu, SIMPLIFY = F)
  Wtheta <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    beta_v <- Reduce(`+`, xi_u[ancestors[[v]]])
    theta_v <- Reduce(`+`, delta[ancestors[[v]]])
    pred_v[outcomes_units[[v]], ] <- X[outcomes_units[[v]], ] %*% beta_v 
    Wtheta[outcomes_units[[v]]] <-  W[outcomes_units[[v]], ] %*% theta_v
  }
  # Update Sigma_g and tau_t_g
  xxT_g_eta <- mapply(`*`, xxT, g_eta, SIMPLIFY = F)
  for (v in 1:p) {
    Sigma_inv[[v]] <- 2 * Reduce(`+`, xxT_g_eta[outcomes_nodes[[v]]]) + 
      diag(1 / tau, nrow = K)
    Sigma[[v]] <- Matrix::solve(Sigma_inv[[v]])
    Sigma_det[v] <- Matrix::det(Sigma[[v]])
    tau_t[v] <- tau
    # update mu_g
    mu[[v]] <- Sigma[[v]] %*%
      Matrix::crossprod( X[outcomes_nodes[[v]], , drop = F],
             y / 2 - 2 * g_eta[outcomes_nodes[[v]]] * 
            (Wtheta[outcomes_nodes[[v]]] + apply(pred_g[outcomes_nodes[[v]], -v, drop = F], 1, sum)))
    # update prob_g (pi_g in manuscript)
    u <- 0.5 * Matrix::crossprod(mu[[v]], Sigma_inv[[v]]) %*% mu[[v]] +
      0.5 * log(Sigma_det[v]) + log(rho / (1 - rho)) - 0.5 * K * log(tau_t[v])
    prob[v] <- expit(u[1, 1])
    # update pred_g
    pred_g[, v] <- prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  # Update non-sparse coefficients ---------------------------------------------------
  # Update Omega only if hyperparameters were updated at last step
  Omega_inv <- 2 * Matrix::crossprod(W, A_eta) %*% W + Matrix::Diagonal(m, 1 / omega)
  if (m != 0) {
    Omega <- solve(Omega_inv)
  }
  if (m == 1) {
    Omega_det <- Omega[1, 1]
  } else {
    Omega_det <- Matrix::det(Omega)
  }
  # Update delta
  delta <- Omega %*% Matrix::crossprod(W, y / 2 - 2 * g_eta * apply(pred_g, 1, sum))
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det))
}