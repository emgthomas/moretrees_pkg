# ------------------------------------------------------------ #
# --------- Performs one step in VI optimization  ------------ #
# ------------------------------------------------------------ #

#' \code{update_vi_logistic_moretrees} Performs variational updates in VI algorithm.
#' 
#' @param X,W,y,outcomes_units,outcomes_nodes,ancestors,levels outputs from 
#' \code{moretrees_design_tree}
#' @param xxT,wwT computed from \code{X} and \code{W} in \code{spike_and_slab_logisitic_moretrees()}
#' @param n,K,p,pL,m,Fg computed from the data by \code{spike_and_slab_logistic_moretrees()}
#' @param prob,mu,Sigma,Sigma_inv,Sigma_det,tau_t,delta,Omega,Omega_inv,Omega_det,a_t,b_t 
#' variational parameters updated by \code{update_vi_params_logistic_moretrees()}
#' @param eta,g_eta,tau,omega parameters updated by \code{update_hyperparams_logistic_moretrees()}
#' @param a,b fixed hyperparameters
#' @family Internal VI functions

update_vi_params_logistic_moretrees <- function(X, W, y,
                                                outcomes_nodes, outcomes_units,
                                                ancestors, levels,
                                                xxT, wwT,
                                                n, p, pL, K, m, Fg, # dsgn
                                                prob, mu, Sigma, 
                                                Sigma_inv, Sigma_det, tau_t, 
                                                delta, Omega, Omega_inv, Omega_det, 
                                                a_t, b_t, # vi_params
                                                eta, g_eta, 
                                                tau, omega, # hyperparams
                                                a, b) { # hyper_fixed

  # Update sparse coefficients ------------------------------------------------------
  xi <- mapply(`*`, prob, mu, SIMPLIFY = F)
  Wtheta <- numeric(n) + 0
  for (u in 1:pL) {
    theta_u <- Reduce(`+`, delta[ancestors[[u]]])
    Wtheta[outcomes_units[[u]]] <- W[outcomes_units[[u]], ] %*% theta_u
  }
  xxT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                      xxT = xxT, g_eta = g_eta, K = K)
  for (v in 1:p) {
    leaf_descendants <- outcomes_nodes[[v]]
    # Update Sigma_v and tau_t_v
    tau_t[v] <- tau[levels[v]]
    Sigma_inv[[v]] <- 2 * Reduce(`+`, xxT_g_eta[leaf_descendants]) + 
      diag(1 / tau_t[v], nrow = K)
    Sigma[[v]] <- solve(Sigma_inv[[v]])
    Sigma_det[v] <- det(Sigma[[v]])
    # Update mu_v
    mu[[v]] <- matrix(0, nrow = K, ncol = 1)
    for (u in leaf_descendants) {
      anc_u_mv <- setdiff(ancestors[[u]], v)
      beta_u_mv <- Reduce(`+`, xi[anc_u_mv])
      units_u <- outcomes_units[[u]]
      mu[[v]] <- mu[[v]] + crossprod(X[units_u, , drop = FALSE],
                       (y[units_u] / 2 - 2 * g_eta[units_u] * 
                       (X[units_u, , drop = FALSE] %*% beta_u_mv + Wtheta[units_u]))
      )
    }
    mu[[v]] <- Sigma[[v]] %*% mu[[v]]
    # Update u_v
    u_v <-  digamma(a_t[levels[v]]) - digamma(b_t[levels[v]]) +
      crossprod(mu[[v]], Sigma_inv[[v]]) %*% mu[[v]] / 2 +
      (log(Sigma_det[v]) - K * log(tau_t[v])) / 2
    prob[v] <- expit(u_v)
    # Update xi
    xi[[v]] <- prob[v] * mu[[v]]
  }
  
  # Update rho ------------------------------------------------------------------------
  for (f in 1:Fg) {
    a_t[f] <- a[f] + sum(prob[levels == f])
    b_t[f] <- b[f] + sum(1 - prob[levels == f])
  }
  
  # Update non-sparse coefficients ----------------------------------------------------
  if (m > 0) {
    Xbeta <- numeric(n) + 0
    for (u in 1:pL) {
      beta_u <- Reduce(`+`, xi[ancestors[[u]]])
      Xbeta[outcomes_units[[u]]] <- X[outcomes_units[[u]], ] %*% beta_u
    }
    wwT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = wwT, g_eta = g_eta, K = m)
    for (v in 1:p) {
      # Update Omega_v
      leaf_descendants <- outcomes_nodes[[v]]
      omega_t <- omega[levels[v]]
      Omega_inv[[v]] <- 2 * Reduce(`+`, wwT_g_eta[leaf_descendants]) +
        diag(1 / omega_t, nrow = m)
      Omega[[v]] <- solve(Omega_inv[[v]])
      Omega_det[v] <- det(Omega[[v]])
      # Update delta_v
      delta[[v]] <- delta[[v]] * 0
      for (u in leaf_descendants) {
        anc_u_mv <- setdiff(ancestors[[u]], v)
        units_u <- outcomes_units[[u]]
        theta_u_mv <- Reduce(`+`, delta[anc_u_mv])
        delta[[v]] <- delta[[v]] + crossprod(W[units_u, , drop = FALSE],
                     (y[units_u] / 2 - 2 * g_eta[units_u] *
                     (W[units_u, , drop = FALSE] %*% theta_u_mv + Xbeta[units_u])
                                             ) )
      }
      delta[[v]] <- Omega[[v]] %*% delta[[v]]
    }
  }
  
  # Return ---------------------------------------------------------------------------
  return(list(prob = prob, mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
              Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
              Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det,
              a_t = a_t, b_t = b_t))
}