# --------------------------------------------------------------------------------- #
# --------------------- Runs VI algorithm for Gaussian outcome  ------------------- #
# --------------------------------------------------------------------------------- #

moretrees_normal <- function(y, X, tol = 1E-4, update.hyper.freq = 10, 
                             print.freq=10) {
  # Prepare for running algorithm ---------------------------------------------------
  G <- dim(X)[1]
  n <- dim(X)[2]
  K <- dim(X)[3]
  # Computing XtX so we don't have to do this repeatedly 
  XtX <- plyr::aaply(X, 1, function(X) crossprod(X, X), .drop = F)
  attributes(XtX)$dimnames <- NULL
  # Hyperparameter initial values 
  rho <- 0.5
  tau <- var(as.numeric(beta[beta!=0]))
  sigma2 <- var(y)
  # Variational parameter initial values 
  Sigma_inv <- plyr::aaply(.data = XtX, .margins = 1, 
                     .fun = function(XtX, tau, K, sigma2) XtX / sigma2 + diag(1 / tau, K),
                     tau = tau, K = K, sigma2 = sigma2, .drop = F)
  attributes(Sigma_inv)$dimnames <- NULL
  Sigma <- plyr::aaply(.data = Sigma_inv, .margins = 1, .fun = solve, .drop = F)
  attributes(Sigma)$dimnames <- NULL
  if (K == 1) {
    Sigma_det <- Sigma[,1,1]
  } else {
    Sigma_det <- aaply(.data = Sigma, .margins = 1, .fun = det)
  }
  mu <- matrix(rnorm(G*K,sd=10),nrow=G)
  prob <- runif(G)
  tau_t <- tau
  # Put parameters in list
  params <- list(mu = mu, prob = prob, Sigma = Sigma,
                 Sigma_inv = Sigma_inv, Sigma_det = Sigma_det, # variational params
                 tau_t = tau_t, rho = rho, sigma2 = sigma2, tau = tau) # hyperparams
  # Run algorithm -----------------------------------------------------------------
  ELBO.track <- c()
  rho.track <- c()
  tau.track <- c()
  sigma2.track <- c()
  ELBO.old <- -1E16
  ELBO.new <- 1E16
  i <- 0
  while(abs(ELBO.new-ELBO.old)>tol){
    ELBO.old <- ELBO.new
    i <- i+1
    update.hyper <- (i %% update.hyper.freq == 0)
    update.hyper.last <- (i %% update.hyper.freq == 1)
    # update.hyper <- FALSE
    # update.hyper.last <- FALSE
    params <- update_params_normal(X = X, XtX = XtX, y = y, n = n, K = K, G = G,
                                   prob = params$prob, mu = params$mu, Sigma = params$Sigma,
                                   Sigma_inv = params$Sigma_inv, Sigma_det = params$Sigma_det,
                                   tau_t = params$tau_t, sigma2 = params$sigma2, rho = params$rho,
                                   tau = params$tau,
                                   update.hyper = update.hyper, update.hyper.last = update.hyper.last)
    ELBO.new <- params$ELBO
    ELBO.track <- c(ELBO.track, ELBO.new)
    rho.track <- c(rho.track, params$rho)
    tau.track <- c(tau.track, params$tau)
    sigma2.track <- c(sigma2.track, params$sigma2)
    if (i %% print.freq == 0) print(i)
  }
  return(list(params = params, ELBO = ELBO.track, rho.track = rho.track,
              tau.track = tau.track, sigma2.track = sigma2.track))
}






