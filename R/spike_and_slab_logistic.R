#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_logistic} performs group variable selection via a spike
#'   and slab prior. The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Vector of outcomes data.
#' @param X Array of design matrices for each variable group, with dimensions dim(X) = c(G,n,K),
#'   where G is the number of variable groups, n is the number of observations, and
#'   K is the number of variables.
#' @param tol Convergence tolerance for ELBO. Default = 1E-4.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples Add this later from test file.
#' @family spike and slab functions

spike_and_slab_logistic <- function(y, X, W, tol = 1E-4, max_iter = 1E5,
                                  update_hyper = T, update_hyper_freq = 10,
                                  hyperparams_init = list(omega = 100,
                                                          tau = 100,
                                                          rho = 0.5)) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- length(X)
  n <- length(y)
  K <- sapply(X, ncol)
  m <- ncol(W)
  # Initial hyperparameter values
  eta <- abs(rnorm(n))
  g_eta <- gfun(eta)
  hyperparams <- hyperparams_init
  if (m == 0) {
    hyperparams$omega <- 1
  }
  hyperparams$eta <- eta
  hyperparams$g_eta <- g_eta
  # Variational parameter initial values
  A_eta <- Matrix::Diagonal(n = n, g_eta)
  Sigma_inv <- sapply(X = X, 
                FUN = function(X, A, tau) 2 * Matrix::t(X) %*% A %*% X + 
                Matrix::Diagonal(ncol(X), 1 / tau),
                tau = hyperparams$tau,
                A = A_eta)
  Sigma <- sapply(Sigma_inv, Matrix::solve)
  Sigma_det <- sapply(Sigma, Matrix::det)
  mu <- sapply(K, rnorm, mean = 0 , sd = 10, simplify = F)
  mu <- sapply(mu, Matrix::Matrix, ncol = 1)
  prob <- rep(rho, G)
  tau_t <- rep(tau, G)
  delta <- Matrix::Matrix(rnorm(m, sd = 10), ncol = 1)
  Omega_inv <- 2 * Matrix::t(W) %*% A_eta %*% W + 
    Matrix::Diagonal(m, 1 / hyperparams$omega)
  if (m != 0) {
    Omega <- Matrix::solve(Omega_inv)
    Omega_det <- Matrix::det(Omega)
  } else {
    Omega <- Matrix::Matrix(nrow = 0, ncol = 0)
    Omega_det <- 1
  }
  # Put VI parameters in list
  vi_params <- list(mu = mu, prob = prob, Sigma = Sigma,
                    Sigma_inv = Sigma_inv, Sigma_det = Sigma_det,
                    tau_t = tau_t, delta = delta,
                    Omega = Omega, Omega_inv = Omega_inv,
                    Omega_det = Omega_det)
  # Compute initial ELBO
  hyperparams <-  update_hyperparams_logistic(X = X, W = W,
                                              y = y, n = n,
                                              K = K, G = G, m = m,
                                              prob = prob, mu = mu,
                                              Sigma = Sigma, Sigma_det = Sigma_det,
                                              tau_t = tau_t,
                                              delta = delta,
                                              Omega = Omega, Omega_det = Omega_det,
                                              eta = hyperparams$eta,
                                              g_eta = hyperparams$g_eta,
                                              omega = hyperparams$omega,
                                              tau = hyperparams$tau,
                                              rho = hyperparams$rho,
                                              update_hyper = F)
  ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track[1] <- hyperparams$ELBO
  ELBO_track2 <- numeric(max_iter + 1)
  ELBO_track2[1] <- hyperparams$ELBO
  omega_track <- numeric(max_iter %/% update_hyper_freq + 1)
  omega_track[1] <- hyperparams$omega
  rho_track <- numeric(max_iter %/% update_hyper_freq + 1)
  rho_track[1] <- hyperparams$rho
  tau_track <- numeric(max_iter %/% update_hyper_freq + 1)
  tau_track[1] <- hyperparams$tau
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    i <- i + 1
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    vi_params <- update_vi_params_logistic(X = X, W = W,
                                         y = y, n = n, K = K, G = G, m = m,
                                         prob = vi_params$prob, 
                                         mu = vi_params$mu, 
                                         Sigma = vi_params$Sigma, 
                                         Sigma_inv = vi_params$Sigma_inv, 
                                         Sigma_det = vi_params$Sigma_det, 
                                         tau_t = vi_params$tau_t,
                                         delta = vi_params$delta, 
                                         Omega = vi_params$Omega,
                                         Omega_inv = vi_params$Omega_inv, 
                                         Omega_det = vi_params$Omega_det,
                                         eta = hyperparams$eta,
                                         g_eta = hyperparams$g_eta,
                                         omega = hyperparams$omega, 
                                         rho = hyperparams$rho, 
                                         tau = hyperparams$tau)
    if (!update_hyper_i) {
      hyperparams <- update_hyperparams_logistic(X = X, W = W,
                                               y = y, n = n,
                                               K = K, G = G, m = m,
                                               prob = vi_params$prob, 
                                               mu = vi_params$mu,
                                               Sigma = vi_params$Sigma, 
                                               Sigma_det = vi_params$Sigma_det,
                                               tau_t = vi_params$tau_t,
                                               delta = vi_params$delta,
                                               Omega = vi_params$Omega, 
                                               Omega_det = vi_params$Omega_det,
                                               eta = hyperparams$eta,
                                               g_eta = hyperparams$g_eta,
                                               omega = hyperparams$omega,
                                               tau = hyperparams$tau,
                                               rho = hyperparams$rho,
                                               update_hyper = F)
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track2[i + 1] - ELBO_track2[i]) < tol) {
        # If we are not updating hyperparameters, we have converged
        if (!update_hyper) break
        # Otherwise, fill in results until next hyperparameter update
        i2 <- ceiling(i / update_hyper_freq) * update_hyper_freq
        ELBO_track2[(i + 2):(i2 + 1)] <- hyperparams$ELBO
        i <- i2
        update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
        update_hyper_im1 <- (i %% update_hyper_freq == 1)
      }
    }
    # Update hyperparameters
    if (update_hyper_i) {
      hyperparams <- update_hyperparams_logistic(X = X, W = W,
                                                 y = y, n = n,
                                                 K = K, G = G, m = m,
                                                 prob = vi_params$prob, 
                                                 mu = vi_params$mu,
                                                 Sigma = vi_params$Sigma, 
                                                 Sigma_det = vi_params$Sigma_det,
                                                 tau_t = vi_params$tau_t,
                                                 delta = vi_params$delta,
                                                 Omega = vi_params$Omega, 
                                                 Omega_det = vi_params$Omega_det,
                                                 eta = hyperparams$eta,
                                                 g_eta = hyperparams$g_eta,
                                                 omega = hyperparams$omega,
                                                 tau = hyperparams$tau,
                                                 rho = hyperparams$rho,
                                                 update_hyper = T)
      j <- i %/% update_hyper_freq
      ELBO_track[j + 1] <- hyperparams$ELBO
      ELBO_track2[i + 1] <- hyperparams$ELBO
      omega_track[j + 1] <- hyperparams$omega
      rho_track[j + 1] <- hyperparams$rho
      tau_track[j + 1] <- hyperparams$tau
      if (abs(ELBO_track[j + 1] - ELBO_track[j]) < tol) break
    } 
    if (i %% update_hyper_freq == 0) print(i)
    if (i == max_iter) {
      rlang::warn("Maximum number of iterations reached")
      break
    }
  }
  j <- i %/% update_hyper_freq
  return(list(vi_params = vi_params, hyperparams = hyperparams,
              ELBO_track = ELBO_track[1:(j + 1)], rho_track = rho_track[1:(j + 1)],
              tau_track = tau_track[1:(j + 1)], omega_track = omega_track[1:(j + 1)],
              ELBO_track2 = ELBO_track2[1:(i + 1)]))
}