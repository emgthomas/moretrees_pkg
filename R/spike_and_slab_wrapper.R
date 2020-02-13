# --------------------------------------------------------------------------------- #
# ----------- spike & slab variable selection wrapper function -------------------- #
# --------------------------------------------------------------------------------- #

#' Group spike and slab variable selection with bernoulli or gaussian outcome
#' 
#' Here's a brief description.
#'   \code{spike_and_slab_normal} performs group variable selection via a spike
#'   and slab prior. The posterior is approximated via variational inference.
#'   This function returns coefficient estimates and 95% credible intervals.
#' 
#' All the details go here!
#' 
#' @section Model Description:
#'   Describe group spike and slab prior and all parameters here.
#' 
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X List of length G of design matrices for each variable group.
#'   Each element of the list has dimension n x K_g
#'   where n is the number of observations, and K_g is the number of variables in group g.
#' @param W Matrix of data with non-sparse regression coefficients of dimension n x m
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95% credible interval.
#' @param tol Convergence tolerance for ELBO. Default = 1E-8.
#' @param update_hyper Update hyperparameters? Default = TRUE.
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param print_freq How often to print out iteration number. 
#' Default = every 50 iterations.
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples Add this later from test file.
#' @family spike and slab functions


spike_and_slab <- function(y, X, W = NULL, 
                                  family = "bernoulli",
                                  ci_level = 0.95,
                                  tol = 1E-4, max_iter = 1E5,
                                  update_hyper = T, update_hyper_freq = 10,
                                  hyperparams_init = NULL) {
  # Get correct function
  if (!(family %in% c("bernoulli", "gaussian"))) {
    stop("family must be a string (\"bernoulli\" or \"gaussian\")")
  }
  if (family == "bernoulli") {
    ss_fun <- spike_and_slab_logistic
    y[y == 0] <- -1
  }
  if (family == "gaussian") ss_fun <- spike_and_slab_normal
  
  # Run algorithm
  mod <- ss_fun(y = y, X = X, W = W, 
             tol = tol, 
             max_iter = max_iter,
             update_hyper = update_hyper,
             update_hyper_freq = update_hyper_freq,
             hyperparams_init = hyperparams_init)
  
  # Compute estimates and credible intervals
  beta_est <- compute_betas(mod, ci_level)
  theta_est <- compute_thetas(mod, ci_level)
  
  # Return
  return(list(sparse_est = beta_est, nonsparse_est = theta_est, mod = mod))
}
