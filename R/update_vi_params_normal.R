# --------------------------------------------------------------------------------- #
# --------- Performs one step in VI optimization for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

update_vi_params_normal <- function(X, XtX, y, n, K, G, # data
                                 prob, mu, Sigma, Sigma_inv, Sigma_det, tau_t, # variational params
                                 sigma2, rho, tau, # hyperparams
                                 update_hyper_last) { 
  # Update tau_t ---------------------------------------------------------------------
  # only changes if tau has been updated; see below
  # Update Sigma ---------------------------------------------------------------------
  if (update_hyper_last) {
    # Sigma only needs to be updated if hyperparams were updated last step
    Sigma_inv <- plyr::aaply(.data = XtX / sigma2, .margins = 1,
                             .fun = function(A, tau, K) A + diag(1 / tau, K),
                             tau = tau, K = K, .drop = F)
    Sigma <- plyr::aaply(.data = Sigma_inv, .margins = 1, .fun = solve, .drop = F)
    if (K == 1) {
      Sigma_det <- Sigma[, 1, 1]
    } else {
      Sigma_det <- plyr::aaply(.data = Sigma, .margins = 1, .fun = det)
    }
  }
  # Update mu -----------------------------------------------------------------------
  pred_g <- matrix(0, nrow = n, ncol = G)
  for (g in 2:G) {
    pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
  }
  for (g in 1:G) {
    mu[g, ] <- (1 / sigma2) * matrix(Sigma[g, , ], nrow = K) %*%
      crossprod(matrix(X[g, , ], nrow = n), y - rowSums(pred_g[, -g, drop = F]))
    if (g != G) {
      pred_g[, g] <- prob[g] * matrix(X[g, , ], nrow = n) %*% matrix(mu[g, ], nrow = K)
    }
  }
  # Update prob (pi in manuscript) --------------------------------------------------
  const <- log(rho / (1 - rho)) - 0.5 * K * log(tau_t)
  for (g in 1:G) {
    u <- crossprod(matrix(mu[g, ], nrow = K), matrix(Sigma_inv[g, , ],nrow = K)) %*%
      matrix(mu[g, ], nrow = K) + 0.5 * log(Sigma_det[g]) + const
    prob[g] <- exp(loglogit(u[1, 1]))
  }

  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, sigma2 = sigma2, 
              rho = rho, tau = tau))
}