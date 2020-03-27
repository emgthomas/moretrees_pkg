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
#' @param vi_params_init A list with any starting values supplied by the user for
#' VI parameters. NULL values will be filled in.
#' @param hyperparams_init A list with any starting values supplied by the user for
#' hyperparameters. NULL values will be filled in.
#' @return A list containing starting values for both VI and hyper parameters
#' @examples 
#' @family MOReTreeS functions

moretrees_init_logistic <- function(X, W, y, A,
                                    outcomes_units,
                                    outcomes_nodes,
                                    ancestors,
                                    levels,
                                    xxT, wwT,
                                    vi_params,
                                    hyperparams,
                                    hyper_fixed,
                                    random_init,
                                    random_init_vals) {
  
  n <- length(y)
  m <- ncol(W)
  p <- length(unique(unlist(ancestors)))
  pL <- length(ancestors)
  K <- ncol(X)
  Fg <- max(levels)
  
  # Get coefficient estimates from maximum likelihood ----------------------------------
  if (is.null(vi_params[["mu"]]) | is.null(vi_params[["delta"]])) {
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
        theta_ml[v, ] <- mod$coefficients[(K + 1):(K + m)]
      }
    }
    # replace any NA vals with zero
    beta_ml[is.na(beta_ml)] <- 0
    theta_ml[is.na(theta_ml)] <- 0
    # transform to get initial values of mu and delta
    A_inv <- solve(A)
    mu <- A_inv %*% beta_ml
    delta <- A_inv %*% theta_ml
  }
  if (is.null(vi_params[["mu"]])) {
    vi_params$mu <- lapply(1:p, function(v, mu) matrix(mu[v, ], ncol = 1),
                           mu = mu)
  } else {
    check <- is.list(vi_params$mu) &&
      sapply(vi_params$mu, function(x) all.equal(dim(x), c(K, 1)))
    if (!check) stop("Incompatible initial value supplied for mu")

  }
  if (random_init) {
    vi_params$mu <- lapply(vi_params$mu,
                          function(mu) mu + rnorm(nrow(mu), 
                          sd = abs(mu) * random_init_vals$mu_sd_frac))
  }
  if (is.null(vi_params[["delta"]])) {
    vi_params$delta <- lapply(1:p, function(v, delta) matrix(delta[v, ], ncol = 1),
                              delta = delta)
  } else {
    check <- is.list(vi_params$delta) &&
      sapply(vi_params$delta, function(x) all.equal(dim(x), c(m, 1)))
    if (!check) stop("Incompatible initial value supplied for delta")
  }
  if (random_init) {
    vi_params$delta <- lapply(vi_params$delta,
       function(delta) delta + rnorm(nrow(delta), 
       sd = abs(delta) * random_init_vals$delta_sd_frac))
  }
  
  # Initial values for hyperparms to be updated via EB --------------------------------------
  if (is.null(hyperparams[["tau"]])) {
    hyperparams$tau <- sapply(1:Fg, function(l) mean(unlist(vi_params$mu[levels == l]) ^ 2))
  } else {
    check <- is.numeric(hyperparams$tau) &&
      length(hyperparams$tau) == Fg
    if (!check) stop("Incompatible initial value supplied for tau")
  }
  if (random_init) {
    hyperparams$tau <- sapply(hyperparams$tau, 
       function(tau) runif(1, min = tau * random_init_vals$tau_lims[1],
          max = tau * random_init_vals$tau_lims[2]))
  }
  if (m > 0) {
    if (is.null(hyperparams[["omega"]])) {
      hyperparams$omega <- sapply(1:Fg, function(l) mean(unlist(vi_params$delta[levels == l]) ^ 2))
    } else {
      check <- is.numeric(hyperparams$omega) &&
        length(hyperparams$omega) == Fg
      if (!check) stop("Incompatible initial value supplied for omega")
    }
    if (random_init) {
      hyperparams$omega <- sapply(hyperparams$omega, 
             function(omega) runif(1, min = omega * random_init_vals$omega_lims[1],
             max = omega * random_init_vals$omega_lims[2]))
    }
  } else {
    hyperparams$omega <- rep(1 , Fg)
  }
  
  # Set initial values for hyperpriors ------------------------------------------------------
  if (is.null(vi_params[["prob"]])) {
    vi_params$prob <- rep(0.95, p)
  } else {
    check <- is.numeric(vi_params$prob) &&
      length(vi_params$prob) == p &&
      sum(vi_params$prob >= 0) == p &&
      sum(vi_params$prob <= 1) == p
    if (!check) stop("Incompatible initial value supplied for prob")
  }
  if (random_init) {
    prob <- vi_params$prob
    prob[prob > 0.99] <- 0.99
    prob[prob < 0.01] <- 0.01
    u <- log(prob / (1 - prob))
    u <- u + rnorm(p) * random_init_vals$u_sd_frac * abs(u)
    vi_params$prob <- moretrees::expit(u)
  }
  if (is.null(vi_params[["a_t"]])) {
    vi_params$a_t <- numeric(Fg)
    for (f in 1:Fg) {
      # initialise these parameters using VI updates
      vi_params$a_t[f] <- hyper_fixed$a[f] + sum(vi_params$prob[levels == f]) 
    }
  } else {
    check <- is.numeric(vi_params$a_t) &&
      length(vi_params$a_t) == Fg &&
      sum(vi_params$a_t > 0) == Fg
    if (!check) stop("Incompatible initial value supplied for a_t")
  }
  if (is.null(vi_params[["b_t"]])) {
    vi_params$b_t <- numeric(Fg)
    for (f in 1:Fg) {
      # initialise these parameters using VI updates
      vi_params$b_t[f] <- hyper_fixed$b[f] + sum(1 - vi_params$prob[levels == f]) 
    }
  } else {
    check <- is.numeric(vi_params$b_t) &&
      length(vi_params$b_t) == Fg &&
      sum(vi_params$b_t > 0) == Fg
    if (!check) stop("Incompatible initial value supplied for b_t")
  }
  
  # Get starting values for eta --------------------------------------------------------
  # Use expected linear predictor squared 
  # (this is close to the real update for eta)
  if (is.null(vi_params[["eta"]])) {
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
  } else {
    check <- is.numeric(hyperparams$eta) &&
      length(hyperparams$eta) == n
    check <- is.numeric(vi_params$a_t) &&
      length(vi_params$a_t) == Fg &&
      sum(vi_params$a_t > 0) == Fg
    if (!check) stop("Incompatible initial value supplied for eta")
  }
  if (random_init) {
    hyperparams$eta <- abs(hyperparams$eta * 
        (1 + rnorm(length(hyperparams$eta)) * random_init_vals$eta_sd_frac))
  }
  hyperparams$g_eta <- gfun(hyperparams$eta)
  
  # Sigma and Omega initial values ---------------------------------------------------
  if (is.null(vi_params[["tau_t"]])) {
    vi_params$tau_t <- hyperparams$tau[levels]
  } else {
    check <- is.numeric(vi_params$tau_t) &&
      length(vi_params$tau_t) == p &&
      sum(vi_params$tau_t > 0) == p
    if (!check) stop("Incompatible initial value supplied for tau_t")
  }
  if (is.null(vi_params[["Sigma_inv"]])) {
    xxT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                        xxT = xxT, g_eta = hyperparams$g_eta, K = K)
    vi_params$Sigma_inv <- lapply(X = 1:length(outcomes_nodes), 
                                  FUN = function(v, outcomes, x, K, tau_t) 2 * 
                                    Reduce(`+`, x[outcomes[[v]]]) + 
                                    diag(1 / tau_t[v], nrow = K),
                                  outcomes = outcomes_nodes,
                                  x = xxT_g_eta,
                                  K = K,
                                  tau_t = vi_params$tau_t)
  } else {
    check <- is.list(vi_params$Sigma_inv) &&
      sum(sapply(vi_params$Sigma_inv, is.matrix)) == p &&
      sum(sapply(vi_params$Sigma_inv, function(x) all.equal(dim(x), c(K, K)))) == p &&
      sum(sapply(vi_params$Sigma_inv, function(x) all.equal(t(x), x))) == p
    if (!check) stop("Incompatible initial value supplied for Sigma_inv")
  }
  if (is.null(vi_params[["Sigma"]])) {
    vi_params$Sigma <- lapply(vi_params$Sigma_inv, solve)
  } else {
    check <- is.list(vi_params$Sigma) &&
      sum(sapply(vi_params$Sigma, is.matrix)) == p &&
      sum(sapply(vi_params$Sigma, function(x) all.equal(dim(x), c(K, K)))) == p &&
      sum(sapply(vi_params$Sigma, function(x) all.equal(t(x), x))) == p
    if (!check) stop("Incompatible initial value supplied for Sigma")
  }
  if (is.null(vi_params[["Sigma_det"]])) {
    vi_params$Sigma_det <- sapply(vi_params$Sigma, det)
  } else {
    check <- is.numeric(vi_params$Sigma_det) &&
      length(vi_params$Sigma_det) == p &&
      sum(vi_params$Sigma_det > 0) == p
    if (!check) stop("Incompatible initial value supplied for Sigma_det")
  }
  if (m > 0) {
    if (is.null(vi_params[["Omega_inv"]])) {
      wwT_g_eta <- lapply(X = outcomes_units, FUN = xxT_g_eta_fun,
                          xxT = wwT, g_eta = hyperparams$g_eta, K = m)
      omega_t <- hyperparams$omega[levels]
      vi_params$Omega_inv <- lapply(X = 1:length(outcomes_nodes), 
                                    FUN = function(v, outcomes, w, m, omega_t) 2 * 
                                      Reduce(`+`, w[outcomes[[v]]]) + 
                                      diag(1 / omega_t[v], nrow = m),
                                    outcomes = outcomes_nodes,
                                    w = wwT_g_eta,
                                    m = m,
                                    omega_t = omega_t)
    } else {
      check <- is.list(vi_params$Omega_inv) &&
        sum(sapply(vi_params$Omega_inv, is.matrix)) == p &&
        sum(sapply(vi_params$Omega_inv, function(x) all.equal(dim(x), c(K, K)))) == p &&
        sum(sapply(vi_params$Omega_inv, function(x) all.equal(t(x), x))) == p
      if (!check) stop("Incompatible initial value supplied for Omega_inv")
    }
    if (is.null(vi_params[["Omega"]])) {
      vi_params$Omega <- sapply(vi_params$Omega_inv, solve, simplify = F)
    } else {
      check <- is.list(vi_params$Omega) &&
        sum(sapply(vi_params$Omega, is.matrix)) == p &&
        sum(sapply(vi_params$Omega, function(x) all.equal(dim(x), c(K, K)))) == p &&
        sum(sapply(vi_params$Omega, function(x) all.equal(t(x), x))) == p
      if (!check) stop("Incompatible initial value supplied for Omega")
    }
    if (is.null(vi_params[["Omega_det"]])) {
      vi_params$Omega_det <- sapply(vi_params$Omega, det, simplify = T)
    } else {
      check <- is.numeric(vi_params$Omega_det) &&
        length(vi_params$Omega_det) == p &&
        sum(vi_params$Omega_det > 0) == p
      if (!check) stop("Incompatible initial value supplied for Omega_det")
    }
  } else {
    vi_params$Omega <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_inv <- rep(list(matrix(nrow = 0, ncol = 0)), p)
    vi_params$Omega_det <- rep(1, p)
  }
  
  # Make up ELBO
  hyperparams$ELBO <- 1E-16
  return(list(vi_params = vi_params, 
              hyperparams = hyperparams))
}