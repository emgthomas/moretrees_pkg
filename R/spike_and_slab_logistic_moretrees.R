#' Group spike and slab variable selection with Gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_logistic_moretrees} performs group variable selection via a spike
#'   and slab prior for binary data.
#'   The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Integer vector of length n containing outcomes; 1 = success, 0 = failure.
#' @param X Matrix of dimension n x sum(K), where n is the number of units, and
#' K[g] is the number of variables in group g.
#' @param W Matrix of non-sparse regression coefficients of dimension n x m
#' @param model "ss" or "moretrees". If "ss", regular spike and slab variable selection is implemented.
#' If "moretrees", the multi-outcome moretrees model is fitted.
#' @param groups If model = "ss", groups is a list of length G (number of groups), where groups[[g]] 
#' is an integer vector specifying the columns of X that belong to group g. NULL if model = "moretrees".
#' @param outcomes_units: If model = "moretrees", outcomes_units is a list of length equal to the number 
#' of unique outcomes. Each element of the list is an integer vector indicating which units (entries 
#' of y_reord, rows of X_reord) correspond to each outcomes. NULL if model = "ss".
#' @param outcomes_nodes: If model = "moretrees", outcomes_nodes is a list of length equal to the number 
#' of unique nodes. Each element of the list is an integer vector indicating which outcomes/leaves
#' are descendants of each node. NULL if model = "ss".
#' @param ancestors If model = "moretrees", ancestors ia  list of length equal to the number of unique 
#' outcomes. Each element of the list is an integer vector indicating which nodes on the tree (including 
#' leaves) are ancestors of the corresponding outcome. NULL if model = "ss".
#' @param tol Convergence tolerance for ELBO.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. Default = every 10 iterations.
#' @return A list of variational parameters.
#' @examples
#' @family spike and slab functions

spike_and_slab_logistic_moretrees <- function(dsgn,
                                              tol, max_iter,
                                              update_hyper, 
                                              update_hyper_freq,
                                              hyper_fixed,
                                              print_freq,
                                              hyper_random_init,
                                              vi_random_init) {
  if (is.null(dsgn$W)) {
    W <- Matrix::Matrix(nrow = length(dsgn$y), ncol = 0)
  }
  # Prepare for running algorithm ---------------------------------------------------
  n <- length(dsgn$y)
  m <- ncol(dsgn$W)
  p <- length(unique(unlist(dsgn$ancestors)))
  K <- ncol(dsgn$X)
  # Initial hyperparameter values
  eta <- abs(rnorm(n, mean = 0, sd = vi_random_init$eta_sd))
  g_eta <- gfun(eta)
  if (update_hyper) {
    # If hyperparameters will be updated, randomly initialise them
    hyperparams <- list(omega = runif(1, 0, hyper_random_init$omega_max),
                        tau = runif(1, 0, hyper_random_init$tau_max),
                        rho = runif(1, 0, 1))
  } else {
    # Otherwise, use fixed values
    hyperparams <- hyper_fixed
  }
  
  if (m == 0) {
    hyperparams$omega <- 1
  }
  hyperparams$eta <- eta
  hyperparams$g_eta <- g_eta
  # Variational parameter initial values
  # A_eta <- Matrix::Diagonal(n = n, g_eta)
  xxT <- plyr::alply(X, 1, tcrossprod)
  xxT_g_eta <- mapply(`*`, xxT, g_eta, SIMPLIFY = F)
  Sigma_inv <- sapply(X = outcomes_nodes, 
                      FUN = function(outcomes, x, K, tau) 2 * Reduce(`+`, x[outcomes]) + 
                        diag(1 / tau, nrow = K),
                      x = xxT_g_eta,
                      K = K,
                      tau = hyperparams$tau,
                      simplify = F)
  Sigma <- sapply(Sigma_inv, solve, simplify = F)
  Sigma_det <- sapply(Sigma, det, simplify = T)
  mu <- sapply(X = 1:p, FUN = function(i) matrix(rnorm(K), ncol = 1),
               simplify = F)
  prob <- runif(p, 0 , 1)
  tau_t <- rep(hyperparams$tau, G)
  delta <- sapply(X = 1:p, FUN = function(i) matrix(rnorm(m), ncol = 1),
                  simplify = F)
  wwT <- plyr::alply(W, 1, tcrossprod)
  Omega_inv <- sapply(X = outcomes_nodes, 
                      FUN = function(outcomes, w, m, omega) 2 * Reduce(`+`, w[outcomes]) + 
                        diag(1 / omega, nrow = m),
                      w = xxT_g_eta,
                      m = m,
                      omega = hyperparams$omega,
                      simplify = F)
  if (m != 0) {
    Omega <- sapply(Omega_inv, solve, simplify = F)
    Omega_det <- sapply(Omega, det, simplify = T)
  } else {
    Omega <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    Omega_det <- rep(1, p)
  }
  # Put VI parameters in list
  vi_params <- list(mu = mu, prob = prob, Sigma = Sigma,
                    Sigma_inv = Sigma_inv, Sigma_det = Sigma_det,
                    tau_t = tau_t, delta = delta,
                    Omega = Omega, Omega_inv = Omega_inv,
                    Omega_det = Omega_det)
  # Compute initial ELBO
  hyperparams <-  update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                        W = dsgn$W,
                                                        y = dsgn$y, 
                                                        outcomes_units = dsgn$outcomes_units,
                                                        ancestors = dsgn$ancestors,
                                                        n = n, K = K, p = p, m = m,
                                                        prob = prob, mu = mu,
                                                        Sigma = Sigma, Sigma_det = Sigma_det,
                                                        tau_t = tau_t, delta = delta,
                                                        Omega = Omega, Omega_det = Omega_det,
                                                        eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                        omega = hyperparams$omega, tau = hyperparams$tau,
                                                        rho = hyperparams$rho, update_hyper = F)
  ELBO_track <- numeric(max_iter %/% update_hyper_freq + 1)
  ELBO_track[1] <- hyperparams$ELBO
  ELBO_track2 <- numeric(max_iter + 1)
  ELBO_track2[1] <- hyperparams$ELBO
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    if (i >= max_iter) {
      cat(paste("Iteration", i, "complete.\n"))
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
    i <- i + 1
    if (i %% print_freq == 0) cat("Iteration", i, "\n")
    update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
    vi_params <- update_vi_params_logistic_moretrees(y = dsgn$y, X = dsgn$X,
                                                     W = dsgn$W,
                                                     outcomes_nodes = dsgn$outcomes_nodes,
                                                     n = n, K = K, p = p, m = m,
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
      hyperparams <-  update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                            W = dsgn$W,
                                                            y = dsgn$y, 
                                                            outcomes_units = dsgn$outcomes_units,
                                                            ancestors = dsgn$ancestors,
                                                            n = n, K = K, p = p, m = m,
                                                            prob = prob, mu = mu,
                                                            Sigma = Sigma, Sigma_det = Sigma_det,
                                                            tau_t = tau_t, delta = delta,
                                                            Omega = Omega, Omega_det = Omega_det,
                                                            eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                            omega = hyperparams$omega, tau = hyperparams$tau,
                                                            rho = hyperparams$rho,
                                                            update_hyper = F)
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track2[i + 1] - ELBO_track2[i]) < tol) {
        # If we are not updating hyperparameters, we have converged
        if (!update_hyper) break
        # Otherwise, fill in results until next hyperparameter update
        i2 <- ceiling(i / update_hyper_freq) * update_hyper_freq
        if (i2 >= max_iter) {
          ELBO_track2[(i + 2):max_iter] <- hyperparams$ELBO
          i <- max_iter
          cat("Iteration", i, "complete.\n")
          cat("\nWarning: Maximum number of iterations reached!\n")
          break
        }
        ELBO_track2[(i + 2):(i2 + 1)] <- hyperparams$ELBO
        i <- i2
        update_hyper_i <- (i %% update_hyper_freq == 0) & update_hyper
        update_hyper_im1 <- (i %% update_hyper_freq == 1)
      }
    }
    # Update hyperparameters
    if (update_hyper_i) {
      hyperparams <-   hyperparams <-  update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                                             W = dsgn$W,
                                                                             y = dsgn$y, 
                                                                             outcomes_units = dsgn$outcomes_units,
                                                                             ancestors = dsgn$ancestors,
                                                                             n = n, K = K, p = p, m = m,
                                                                             prob = prob, mu = mu,
                                                                             Sigma = Sigma, Sigma_det = Sigma_det,
                                                                             tau_t = tau_t, delta = delta,
                                                                             Omega = Omega, Omega_det = Omega_det,
                                                                             eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                                             omega = hyperparams$omega, tau = hyperparams$tau,
                                                                             rho = hyperparams$rho,
                                                                             update_hyper = T)
      j <- i %/% update_hyper_freq
      ELBO_track[j + 1] <- hyperparams$ELBO
      ELBO_track2[i + 1] <- hyperparams$ELBO
      if (abs(ELBO_track[j + 1] - ELBO_track[j]) < tol) break
    }
  }
  return(list(vi_params = vi_params, hyperparams = hyperparams,
              ELBO_track = ELBO_track2[1:(i + 1)]))
}