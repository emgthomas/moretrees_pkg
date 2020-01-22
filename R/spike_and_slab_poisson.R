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

spike_and_slab_poisson <- function(y, X, tol = 1E-4, max_iter = 1E5,
                                  update_hyper = T, update_hyper_freq = 10,
                                  hyperparams_init = NULL,
                                  vi_params_init = NULL) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- dim(X)[1]
  n <- dim(X)[2]
  K <- dim(X)[3]
  # Computing some constant values
  sum_log_y_fac <- sum(sapply(y, logfac))
  Xy <- matrix(nrow = G, ncol = K)
  for (g in 1:G) {
    Xy[g, ] <- t(X[g, , ]) %*% y
  }
  # Initial hyperparameter values
  if (is.null(hyperparams_init)) {
    hyperparams <- list(tau = 100,
         rho = 0.5,
         tau_alpha = 10)
  } else {
    hyperparams <- hyperparams_init
  }
  # Variational parameter initial values
  if (is.null(vi_params_init)) {
    u <- rep(0, G)
    mu <- matrix(rnorm(G * K, sd = 1), nrow = G)
    Sigma <- array(dim = c(G, K, K))
    for (g in 1:G) {
      Sigma[g, , ] <- diag(1, nrow = K)
    }
    mu_alpha <- rnorm(1, sd = 10)
    tau_t_alpha <- 1
    vi_params <- list(u = u, mu = mu, Sigma = Sigma,
                      mu_alpha = mu_alpha, tau_t_alpha = tau_t_alpha)
  } else {
    vi_params <- vi_params_init
  }
  # compute expA
  expA <- matrix(nrow = n, ncol = G)
  b <- numeric(n)
  for (g in 1:G) {
    for (i in 1:n) {
      b[i] <- t(X[g, i, ]) %*% Sigma[g, , ] %*% X[g, i, ]
    }
    expA[, g] <- exp(mu[g, ] %*% t(X[g, , ])  + b / 2)
  }
  vi_params$expA <- expA
  # Compute initial ELBO
  ELBO <-  elbo_poisson(X = X, y = y, n = n, K = K, G = G, sum_log_y_fac = sum_log_y_fac, 
                                      u = vi_params$u,
                                      mu = vi_params$mu, 
                                      Sigma = vi_params$Sigma, 
                                      mu_alpha = vi_params$mu_alpha, 
                                      tau_t_alpha = vi_params$tau_t_alpha, 
                                      expA = vi_params$expA,
                                      rho = hyperparams$rho,
                                      tau = hyperparams$tau,
                                      tau_alpha = hyperparams$tau_alpha)
  # ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track <- numeric(max_iter)
  ELBO_track[1] <- ELBO
  # ELBO_track2 <- numeric(max_iter + 1)
  # ELBO_track2[1] <- hyperparams$ELBO
  # sigma2_track <- numeric(max_iter %/% update_hyper_freq + 1)
  # sigma2_track[1] <- hyperparams$sigma2
  # rho_track <- numeric(max_iter %/% update_hyper_freq + 1)
  # rho_track[1] <- hyperparams$rho
  # tau_track <- numeric(max_iter %/% update_hyper_freq + 1)
  # tau_track[1] <- hyperparams$tau
  # Run algorithm -----------------------------------------------------------------
  i <- 1
  repeat {
    i <- i + 1
    # update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    # update_hyper_im1 <- (i %% update_hyper_freq == 1)
    vi_params <- update_vi_params_poisson(X = X, y = y, Xy = Xy, n = n, K = K, G = G, 
                                          u = vi_params$u,
                                          mu = vi_params$mu, 
                                          Sigma = vi_params$Sigma, 
                                          mu_alpha = vi_params$mu_alpha, 
                                          tau_t_alpha = vi_params$tau_t_alpha, 
                                          expA = vi_params$expA,
                                          rho = hyperparams$rho,
                                          tau = hyperparams$tau,
                                          tau_alpha = hyperparams$tau_alpha)
    ELBO_track[i] <- elbo_poisson(X = X, y = y, n = n, K = K, G = G, sum_log_y_fac = sum_log_y_fac, 
                                  u = vi_params$u,
                                  mu = vi_params$mu, 
                                  Sigma = vi_params$Sigma, 
                                  mu_alpha = vi_params$mu_alpha, 
                                  tau_t_alpha = vi_params$tau_t_alpha, 
                                  expA = vi_params$expA,
                                  rho = hyperparams$rho,
                                  tau = hyperparams$tau,
                                  tau_alpha = hyperparams$tau_alpha)
    # if (!update_hyper_i) {
    #   hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, y = y,
    #                                            n = n, K = K, G = G,
    #                                            prob = vi_params$prob,
    #                                            mu = vi_params$mu,
    #                                            Sigma = vi_params$Sigma,
    #                                            Sigma_det = vi_params$Sigma_det, 
    #                                            tau_t = vi_params$tau_t,
    #                                            sigma2 = hyperparams$sigma2,
    #                                            tau = hyperparams$tau,
    #                                            rho = hyperparams$rho,
    #                                            update_hyper = F)
    #   ELBO_track2[i + 1] <- hyperparams$ELBO
    #   if (abs(ELBO_track2[i + 1] - ELBO_track2[i]) < tol) {
    #     # If we are not updating hyperparameters, we have converged
    #     if (!update_hyper) break
    #     # Otherwise, fill in results until next hyperparameter update
    #     i2 <- ceiling(i / update_hyper_freq) * update_hyper_freq
    #     ELBO_track2[(i + 2):(i2 + 1)] <- hyperparams$ELBO
    #     i <- i2
    #     update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    #     update_hyper_im1 <- (i %% update_hyper_freq == 1)
    #   }
    # }
    # # Update hyperparameters
    # if (update_hyper_i) {
    #   hyperparams <- update_hyperparams_normal(X = X, XtX = XtX, y = y,
    #                                            n = n, K = K, G = G,
    #                                            prob = vi_params$prob,
    #                                            mu = vi_params$mu,
    #                                            Sigma = vi_params$Sigma,
    #                                            Sigma_det = vi_params$Sigma_det, 
    #                                            tau_t = vi_params$tau_t,
    #                                            sigma2 = hyperparams$sigma2,
    #                                            tau = hyperparams$tau,
    #                                            rho = hyperparams$rho,
    #                                            update_hyper = T)
    #   j <- i %/% update_hyper_freq
    #   ELBO_track[j + 1] <- hyperparams$ELBO
    #   ELBO_track2[i + 1] <- hyperparams$ELBO
    #   sigma2_track[j + 1] <- hyperparams$sigma2
    #   rho_track[j + 1] <- hyperparams$rho
    #   tau_track[j + 1] <- hyperparams$tau
    #   if (abs(ELBO_track[j + 1] - ELBO_track[j]) < tol) break
    # } 
    if (i %% update_hyper_freq == 0) print(i)
    if (abs(ELBO_track[i] - ELBO_track[i - 1]) < tol) break
    if (i == max_iter) {
      rlang::warn("Maximum number of iterations reached")
      break
    }
  }
  # j <- i %/% update_hyper_freq
  j <- i - 1
  return(list(vi_params = vi_params, hyperparams = hyperparams,
              ELBO_track = ELBO_track[1:(j + 1)]))
              # rho_track = rho_track[1:(j + 1)],
              # tau_track = tau_track[1:(j + 1)], sigma2_track = sigma2_track[1:(j + 1)],
              # ELBO_track2 = ELBO_track2[1:(i + 1)]))
}