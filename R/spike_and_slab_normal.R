#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_normal} performs group variable selection via a spike
#'   and slab prior. The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Vector of outcomes data.
#' @param X List of length g of design matrices for each variable group, each with dimension n x K
#'   where G is the number of variable groups, n is the number of observations, and
#'   K is the number of variables.
#' @param tol Convergence tolerance for ELBO. Default = 1E-4.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @param print_freq How often to print out iteration number. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples Add this later from test file.
#' @family spike and slab functions

spike_and_slab_normal <- function(y, X, W = Matrix::Matrix(nrow = length(y), ncol = 0), 
                                  tol = 1E-4, max_iter = 1E5,
                                  update_hyper = T, update_hyper_freq = 10,
                                  hyperparams_init = list(omega = 100,
                                                          tau = 100,
                                                          sigma2 = var(y),
                                                          rho = 0.5)) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- length(X)
  n <- length(y)
  K <- sapply(X, ncol)
  m <- ncol(W)
  # Computing XtX and WtW so we don't have to do this repeatedly
  XtX <- sapply(X, Matrix::crossprod)
  WtW <- Matrix::crossprod(W)
  # Initial hyperparameter values
  hyperparams <- hyperparams_init
  # Variational parameter initial values
  Sigma_inv <- sapply(XtX, FUN = function(XtX, tau, sigma2) XtX / sigma2 + 
                      Matrix::Diagonal(ncol(XtX), 1 / tau),
                      tau = hyperparams$tau, 
                      sigma2 = hyperparams$sigma2)  
  Sigma <- sapply(Sigma_inv, Matrix::solve)
  Sigma_det <- sapply(Sigma, Matrix::det)
  mu <- sapply(K, rnorm, mean = 0 , sd = 10, simplify = F)
  mu <- sapply(mu, Matrix::Matrix, ncol = 1, simplify = F)
  prob <- rep(rho, G)
  tau_t <- rep(tau, G) # this should not be changed; tau_t = tau according to algorithm
  delta <- Matrix::Matrix(rnorm(m, sd = 10), ncol = 1)
  Omega_inv <- WtW / hyperparams$sigma2 + Matrix::Diagonal(m, 1 / hyperparams$omega)
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
  hyperparams <-  update_hyperparams_normal(X = X, XtX = XtX, 
                                              W = W, WtW = WtW,
                                              y = y, n = n,
                                              K = K, G = G, m = m,
                                              prob = prob, mu = mu,
                                              Sigma = Sigma, Sigma_det = Sigma_det,
                                              tau_t = tau_t,
                                              delta = delta,
                                              Omega = Omega, Omega_det = Omega_det,
                                              omega = hyperparams$omega,
                                              sigma2 = hyperparams$sigma2,
                                              tau = hyperparams$tau,
                                              rho = hyperparams$rho,
                                              update_hyper = F)
  ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track[1] <- hyperparams$ELBO
  ELBO_track2 <- numeric(max_iter + 1)
  ELBO_track2[1] <- hyperparams$ELBO
  sigma2_track <- numeric(max_iter %/% update_hyper_freq + 1)
  sigma2_track[1] <- hyperparams$sigma2
  rho_track <- numeric(max_iter %/% update_hyper_freq + 1)
  rho_track[1] <- hyperparams$rho
  tau_track <- numeric(max_iter %/% update_hyper_freq + 1)
  tau_track[1] <- hyperparams$tau
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    i <- i + 1
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    update_hyper_im1 <- (i %% update_hyper_freq == 1) & update_hyper
    vi_params <- update_vi_params_normal(X = X, XtX = XtX, 
                                         W = W, WtW = WtW,
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
                                         omega = hyperparams$omega, 
                                         sigma2 = hyperparams$sigma2,  
                                         rho = hyperparams$rho, 
                                         tau = hyperparams$tau,
                                         update_hyper_last = update_hyper_im1)
    if (!update_hyper_i) {
      hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, 
                                               W = W, WtW = WtW,
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
                                               omega = hyperparams$omega,
                                               sigma2 = hyperparams$sigma2,
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
      hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, 
                                               W = W, WtW = WtW,
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
                                               omega = hyperparams$omega,
                                               sigma2 = hyperparams$sigma2,
                                               tau = hyperparams$tau,
                                               rho = hyperparams$rho,
                                               update_hyper = T)
      j <- i %/% update_hyper_freq
      ELBO_track[j + 1] <- hyperparams$ELBO
      ELBO_track2[i + 1] <- hyperparams$ELBO
      sigma2_track[j + 1] <- hyperparams$sigma2
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
              tau_track = tau_track[1:(j + 1)], sigma2_track = sigma2_track[1:(j + 1)],
              ELBO_track2 = ELBO_track2[1:(i + 1)]))
}