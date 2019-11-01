#' ssMOReTreeS with a Gaussian outcome,
#' with model fitting via variational inference.
#' 
#' @param y Vector of outcomes data.
#' @param X Array of design matrices for each variable group.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)

moretrees_normal <- function(y, X, tol = 1E-4, 
                             update_hyper = T,
                             update_hyper_freq = 10,
                             print_freq=10) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- dim(X)[1]
  n <- dim(X)[2]
  K <- dim(X)[3]
  # Computing XtX so we don't have to do this repeatedly
  XtX <- plyr::aaply(X, 1, function(X) crossprod(X, X), .drop = F)
  attributes(XtX)$dimnames <- NULL
  # Hyperparameter initial values
  rho <- 0.5
  tau <- 1
  sigma2 <- var(y)
  # Variational parameter initial values
  Sigma_inv <- plyr::aaply(.data = XtX, .margins = 1,
                     .fun = function(XtX, tau, K, sigma2) XtX / sigma2 + diag(1 / tau, K),
                     tau = tau, K = K, sigma2 = sigma2, .drop = F)
  attributes(Sigma_inv)$dimnames <- NULL
  Sigma <- plyr::aaply(.data = Sigma_inv, .margins = 1, .fun = solve, .drop = F)
  attributes(Sigma)$dimnames <- NULL
  if (K == 1) {
    Sigma_det <- Sigma[, 1, 1]
  } else {
    Sigma_det <- plyr::aaply(.data = Sigma, .margins = 1, .fun = det)
  }
  mu <- matrix(rnorm(G * K, sd = 10), nrow = G)
  prob <- runif(G)
  tau_t <- tau
  # Put arguments in list
  args <- list(mu = mu, prob = prob, Sigma = Sigma,
                 Sigma_inv = Sigma_inv, Sigma_det = Sigma_det, # variational params
                 tau_t = tau_t, rho = rho, sigma2 = sigma2, tau = tau, # hyperparams
                 update_hyper = FALSE, update_hyper_last = FALSE) 
  # Run algorithm -----------------------------------------------------------------
  ELBO_track <- c()
  rho_track <- c()
  tau_track <- c()
  sigma2_track <- c()
  ELBO_old <- -1E16
  ELBO_new <- 1E16
  i <- 0
  while (abs(ELBO_new - ELBO_old) > tol) {
    ELBO_old <- ELBO_new
    i <- i + 1
    params$update_hyper <- (i %% update_hyper_freq == 0) & update_hyper
    params$update_hyper_last <- (i %% update_hyper_freq == 1)
    params <- do.call(update_params_normal, params) # update all parameters
    ELBO_new <- params$ELBO
    ELBO_track <- c(ELBO_track, ELBO_new)
    rho_track <- c(rho_track, params$rho)
    tau_track <- c(tau_track, params$tau)
    sigma2_track <- c(sigma2_track, params$sigma2)
    if (i %% print_freq == 0) print(i)
  }
  return(list(params = params, ELBO = ELBO_track, rho_track = rho_track,
              tau_track = tau_track, sigma2_track = sigma2_track))
}