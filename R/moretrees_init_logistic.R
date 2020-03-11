# --------------------------------------------------------------------------------- #
# -------------------- moretrees initial values function -------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees_init_rand} Randomly generates starting values for moretrees
#'   models. Not recommended if the model is converging slowly!
#' 
#' @export
#' @useDynLib moretrees
#' 
#' @section Model Description:
#' Describe MOReTreeS model and all parameters here.
#' 
#' @param dsgn Design list generated by moretrees_design_tree()
#' @param xxT Computed from exposure matrix X
#' @param wwT Computed from covariate matrix W
#' @param hyper_fixed Fixed values of hyperparameters to use if update_hyper = FALSE.
#' If family = "bernoulli", this should be a list including the following elements:
#' tau (prior variance for sparse node coefficients)
#' rho (prior node selection probability for sparse node coefficients)
#' omega (prior variance for non-sparse node coefficients)
#' If family = "gaussian", in addition to the above, the list should also include:
#' sigma2 (variance of residuals)
#' @param hyper_random_init If update_hyper = TRUE, this is a list containing the 
#' maximum values of the hyperparameters. Each hyperparameter will be initialised
#' uniformly at random between 0 and the maximum values given by the list elements
#' below. If multiple random restarts are being used, it is recommended
#' to use a large range for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' tau_max (maxmimum of prior sparse node variance)
#' omega_max (maximum of prior non-sparse node variance)
#' sigma2_max (maximum of residual error variance--- for gaussian data only)
#' @param vi_random_init A list with parameters that determine the distributions from
#' which the initial VI parameters will be randomly chosen. All parameters will be randomly
#' selected from independent normal distributions with the standard deviations given by
#' the list elements below. If multiple random restarts are being used, it is recommended
#' to use large standard deviations for these initial values so that the parameter space
#' can be more effectively explored. The list contains the following elements:
#' mu_sd (standard deviation for posterior means of sparse node coefficients)
#' delta_sd (standard deviation for posterior means of non-sparse node coefficients)
#' xi_sd (standard deviation for auxilliary parameters xi--- for bernoulli data only)
#' @return A list containing starting values
#' @examples 
#' @family MOReTreeS functions

moretrees_init_logistic <- function(X, W, y, A,
                                outcomes_units,
                                outcomes_nodes,
                                ancestors,
                                levels,
                                xxT, wwT,
                                hyper_fixed) {
  
  n <- length(y)
  m <- ncol(W)
  p <- length(unique(unlist(ancestors)))
  pL <- length(ancestors)
  K <- ncol(X)
  Fg <- max(levels)
  vi_params <- list()
  hyperparams <- list()
  
  # Get coefficient estimates from maximum likelihood ----------------------------------
  beta_ml <- matrix(0, nrow = p, ncol = K)
  theta_ml <- matrix(0, nrow = p, ncol = m)
  for (v in 1:p) {
    u <- outcomes_nodes[[v]]
    units <- unlist(outcomes_units[u])
    suppressWarnings(suppressMessages(
      if (m > 0){
        mod <- glm(y[units] == 1 ~ 0 + X[units,  , drop = F] 
                   + W[units,  , drop = F],
                   family = "binomial")
      } else {
        mod <- glm(y[units] == 1 ~ 0 + X[units,  , drop = F],
                   family = "binomial")
      }
    ))
    beta_ml[v, ] <- mod$coefficients[1:K]
    if (m > 0) {
      theta_ml[v, ] <- mod$coefficients[(K+1):(K + m)]
    }
  }
  # replace any NA vals with zero
  beta_ml[is.na(beta_ml)] <- 0
  theta_ml[is.na(theta_ml)] <- 0
  # transform to get initial values of mu and delta
  A_inv <- solve(A)
  mu <- A_inv %*% beta_ml
  delta <- A_inv %*% theta_ml
  vi_params$mu <- lapply(1:p, function(v, mu) matrix(mu[v, ], ncol = 1),
                                        mu = mu)
  vi_params$delta <- lapply(1:p, function(v, delta) matrix(delta[v, ], ncol = 1),
                                           delta = delta)
  
  # Choose fixed hyperparameters ------------------------------------------------------------
  if (is.null(hyper_fixed)) {
    hyper_fixed$a_tau <- sapply(1:max(levels), function(l) sum(levels == l)) * K / 2
    hyper_fixed$b_tau <- sapply(1:max(levels), function(l) sum(mu[levels == l, ] ^ 2)) / 2
    if (m > 0) {
      hyper_fixed$a_omega <- sapply(1:max(levels), function(l) sum(levels == l)) * m / 2
      hyper_fixed$b_omega <- sapply(1:max(levels), function(l) sum(delta[levels == l, ] ^ 2)) / 2
    }
  }
  
  # Set initial values for hyperpriors ------------------------------------------------------
  vi_params$prob <- rep(0.95, p)
  vi_params$a <- numeric(Fg)
  vi_params$b <- numeric(Fg)
  vi_params$a_t_tau <- numeric(Fg)
  vi_params$b_t_tau <- numeric(Fg)
  for (f in 1:Fg) {
    # need to initialise these parameters using VI updates
    # so that terms cancel in ELBO.
    vi_params$a[f] <- 1 + sum(vi_params$prob[levels == f]) 
    vi_params$b[f] <- 1 + sum(levels == f) - sum(vi_params$prob[levels == f]) 
    vi_params$a_t_tau[f] <- sum(levels == f) * K / 2 + hyper_fixed$a_tau[f]
    vi_params$b_t_tau[f] <- sum(mu[levels == f, ] ^ 2) / 2 + hyper_fixed$b_tau[f]
  }
  
  if (m > 0) {
    vi_params$a_t_omega <- numeric(Fg)
    vi_params$b_t_omega <- numeric(Fg)
    for (f in 1:Fg) {
      vi_params$a_t_omega[f] <- sum(levels == f) * m / 2 + hyper_fixed$a_omega[f]
      vi_params$b_t_omega[f] <- sum(delta[levels == f, ] ^ 2) / 2 + hyper_fixed$b_omega[f]
    }
  }
  
  # Get starting values for eta --------------------------------------------------------
  # Use expected linear predictor squared 
  # (this is close to the real update for eta)
  xi <- mapply(FUN = function(prob, mu) prob * mu,
               prob = vi_params$prob, mu = vi_params$mu, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:pL) {
    beta_v <- Reduce(`+`, xi[ancestors[[v]]])
    theta_v <- Reduce(`+`, vi_params$delta[ancestors[[v]]])
    lp[outcomes_units[[v]]] <- X[outcomes_units[[v]], , drop = F] %*% beta_v +
      W[outcomes_units[[v]], , drop = F ] %*% theta_v
  }
  hyperparams$eta <- abs(lp)
  hyperparams$g_eta <- gfun(hyperparams$eta)

  # Sigma and Omega initial values ---------------------------------------------------
  xxT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                      xxT = xxT, g_eta = hyperparams$g_eta, K = K)
  vi_params$tau_t <- (vi_params$b_t_tau / vi_params$a_t_tau)[levels]
  vi_params$Sigma_inv <- lapply(X = 1:length(outcomes_nodes), 
                      FUN = function(v, outcomes, x, K, tau_t) 2 * 
                        Reduce(`+`, x[outcomes[[v]]]) + 
                        diag(1 / tau_t[v], nrow = K),
                      outcomes = outcomes_nodes,
                      x = xxT_g_eta,
                      K = K,
                      tau_t = vi_params$tau_t)
  vi_params$Sigma <- lapply(vi_params$Sigma_inv, solve)
  vi_params$Sigma_det <- sapply(vi_params$Sigma, det)
  if (m > 0) {
    wwT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = wwT, g_eta = hyperparams$g_eta, K = m)
    omega_t <- (vi_params$a_t_omega / vi_params$b_t_omega)[levels]
    vi_params$Omega_inv <- lapply(X = 1:length(outcomes_nodes), 
                                  FUN = function(v, outcomes, w, m, omega_t) 2 * 
                                    Reduce(`+`, w[outcomes[[v]]]) + 
                                    diag(1 / omega_t[v], nrow = m),
                                  outcomes = outcomes_nodes,
                                  w = wwT_g_eta,
                                  m = m,
                                  omega_t = omega_t)
    vi_params$Omega <- sapply(vi_params$Omega_inv, solve, simplify = F)
    vi_params$Omega_det <- sapply(vi_params$Omega, det, simplify = T)
  } else {
    vi_params$Omega <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_inv <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_det <- rep(1, p)
  }
  
  # Make up ELBO
  hyperparams$ELBO <- 1E-16
  return(list(vi_params = vi_params, 
              hyperparams = hyperparams,
              hyper_fixed = hyper_fixed))
}