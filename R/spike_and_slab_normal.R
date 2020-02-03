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
#' @param X Array of design matrices for each variable group, with dimensions dim(X) = c(G,n,K),
#'   where G is the number of variable groups, n is the number of observations, and
#'   K is the number of variables.
#' @param tol Convergence tolerance for ELBO. Default = 1E-4.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @param print_freq How often to print out iteration number. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples Add this later from test file.
#' @family spike and slab functions

spike_and_slab_normal <- function(y, X, tol = 1E-4, max_iter = 1E5,
                                  update_hyper = T, update_hyper_freq = 10,
                                  hyperparams_init = list(tau = 100,
                                                          sigma2 = var(y),
                                                          rho = 0.5)) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- dim(X)[1]
  n <- dim(X)[2]
  K <- dim(X)[3]
  # Computing XtX so we don't have to do this repeatedly
  XtX <- plyr::aaply(X, 1, function(X) crossprod(X, X), .drop = F)
  attributes(XtX)$dimnames <- NULL
  # Initial hyperparameter values
  hyperparams <- hyperparams_init
  # Variational parameter initial values
  Sigma_inv <- plyr::aaply(.data = XtX, .margins = 1,
                           .fun = function(XtX, tau, K, sigma2) XtX / sigma2 + diag(1 / tau, K),
                           tau = hyperparams$tau, K = K, 
                           sigma2 = hyperparams$sigma2, .drop = F)
  attributes(Sigma_inv)$dimnames <- NULL
  Sigma <- plyr::aaply(.data = Sigma_inv, .margins = 1, .fun = solve, .drop = F)
  attributes(Sigma)$dimnames <- NULL
  if (K == 1) {
    Sigma_det <- Sigma[, 1, 1]
  } else {
    Sigma_det <- plyr::aaply(.data = Sigma, .margins = 1, .fun = det)
  }
  mu <- matrix(rnorm(G * K, sd = 10), nrow = G)
  prob <- rep(hyperparams$rho, G)
  tau_t <- rep(hyperparams$tau, G)
  # Put VI parameters in list
  vi_params <- list(mu = mu, prob = prob, Sigma = Sigma,
                    Sigma_inv = Sigma_inv, Sigma_det = Sigma_det,
                    tau_t = tau_t)
  # Compute initial ELBO
  hyperparams <-  update_hyperparams_normal(X = X, XtX = XtX, y = y, n = n,
                                              K = K, G = G, prob = prob, mu = mu,
                                              Sigma = Sigma, Sigma_det = Sigma_det,
                                              tau_t = tau_t,
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
    update_hyper_im1 <- (i %% update_hyper_freq == 1)
    vi_params <- update_vi_params_normal(X = X, XtX = XtX, y = y, n = n, K = K, G = G,
                                         prob = vi_params$prob, mu = vi_params$mu, 
                                         Sigma = vi_params$Sigma, Sigma_inv = vi_params$Sigma_inv, 
                                         Sigma_det = vi_params$Sigma_det, tau_t = vi_params$tau_t, 
                                         sigma2 = hyperparams$sigma2, rho = hyperparams$rho, 
                                         tau = hyperparams$tau,
                                         update_hyper_last = update_hyper_im1)
    if (!update_hyper_i) {
      hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, y = y,
                                               n = n, K = K, G = G,
                                               prob = vi_params$prob,
                                               mu = vi_params$mu,
                                               Sigma = vi_params$Sigma,
                                               Sigma_det = vi_params$Sigma_det, 
                                               tau_t = vi_params$tau_t,
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
      hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, y = y,
                                               n = n, K = K, G = G,
                                               prob = vi_params$prob,
                                               mu = vi_params$mu,
                                               Sigma = vi_params$Sigma,
                                               Sigma_det = vi_params$Sigma_det, 
                                               tau_t = vi_params$tau_t,
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