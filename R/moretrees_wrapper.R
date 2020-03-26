# --------------------------------------------------------------------------------- #
# ------------------------- moretrees wrapper function ---------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees} Fits Multi-Outcome Regression with Tree-structured Shrinkage
#'   (MOReTreeS) model to normally-distributed or binary outcome data.
#'   The posterior is approximated via variational inference.
#'   Returns coefficient estimates and 95% credible intervals.
#' 
#' All the details go here!
#' 
#' @export
#' @useDynLib moretrees
#' 
#' @section Model Description:
#' Describe MOReTreeS model and all parameters here.
#' 
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X An n x K matrix of exposure data, where K is the dimension of the exposure.
#' Grouping of the outcomes will be based on their relationships with the variables in X.
#' @param W Matrix of covariates of dimension n x m.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' @param outcomes Character vector of length n. outcomes[i] is a string indicating the 
#' outcome experienced by unit i.
#' @param tr A directed igraph object. This is a tree representing the relationships
#' among the outcomes. The leaves represent individual outcomes, and internal nodes
#' represent outcome categories consisting of their leaf descendants. All nodes
#' of tr must have unique names as given by names(V(tr)). The names of the leaves must 
#' be equal to the unique elements of outcomes. The vertices of tr, V(tr), may have 
#' an attribute 'levels' containing integer values from 1 to max(V(tr)$levels). 
#' In this case, the levels attribute specifies groups of nodes that share common 
#' hyperparameters rho[f], tau[f], and omega[f]. If V(tr)$levels is NULL, 
#' levels will be automatically chosen based on distance from the root node of the tree.
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval. 
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval
#' @param get_ml If TRUE, moretrees will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model. The default is FALSE.
#' @param tol Convergence tolerance for ELBO. Default = 1E-8.
#' @param tol_hyper If hyper_method = "EB", the convergence tolerance for ELBO
#' between subsequent hyperparmeter updates. Typically a more generous
#' tolerance than tol. Default = 1E-6.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' @param hyper_fixed Fixed values of hyperprior parameters for rho.
#' If family = "bernoulli", this should be a list including the following elements:
#' a_rho, b_rho (parameters of beta prior on rho for each level)
#' If family = "gaussian", in addition to the above, the list should also include:
#' sigma2 (variance of residuals)
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param random_init If TRUE, initial values will be randomly permuted.
#' @param random_init_vals If random_init = TRUE, 
#' this is a list containing parameters for randomly permuting the inital values.
#' The list contains the following elements:
#' @param print_freq How often to print out iteration number and current value of epsilon
#' (the difference in objective function value for the two most recent iterations). 
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that 
#' gives the highest ELBO will be returned. It is recommended to choose nrestarts > 1.
#' The default is 3.
#' @param keep_restarts If TRUE, the results from all random restarts will be returned.
#' If FALSE, only the restart with the highest ELBO is returned.
#' @param parallel If TRUE, the random restarts will be run in parallel. It is recommended
#' to first set the number of cores using doParallel::registerDoParallel(). Otherwise,
#' the default number of cores specified by the doParallel package will be used.
#' @param log_restarts If TRUE, progress of each random restart will be logged to a text
#' file in log_dir.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples 
#' @family MOReTreeS functions

moretrees <- function(X, W = NULL, y, outcomes, tr,
                      initial_values = NULL,
                      ci_level = 0.95,
                      get_ml = FALSE,
                      update_hyper_freq = 50,
                      print_freq = update_hyper_freq,
                      hyper_fixed = NULL,
                      tol = 1E-8, 
                      tol_hyper = 1E-6,
                      max_iter = 5000,
                      nrestarts = 3,
                      keep_restarts = nrestarts > 1,
                      parallel = nrestarts > 1,
                      log_restarts = nrestarts > 1,
                      log_dir = getwd(),
                      random_init = nrestarts > 1,
                      random_init_vals = list(omega_lims = c(0.5, 1.5),
                                              tau_lims = c(0.5, 1.5),
                                              eta_sd_frac = 0.2,
                                              mu_sd_frac = 0.2,
                                              delta_sd_frac = 0.2)) {
  
  if (!is.matrix(X)) stop("X must be a matrix")
  if (!is.null(W) & !(is.matrix(W))) stop("If W is not NULL, must be a matrix")
  if (!(length(get_ml) == 1 & is.logical(get_ml))) stop("get_ml must be either TRUE or FALSE")
  
  # Get MOReTreeS design elements
  dsgn <- moretrees_design_tree(X = X, W = W, y = y, outcomes = outcomes, tr = tr)
  
  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }
  
  # Run algorithm
  args <- sapply(as.list(stackoverflow::match.call.defaults()), eval)
  mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
    if (log_restarts) {
      sink(file = paste0(log_dir, "restart_", i, "_log.txt"))
      cat("Initialising random restart", i, "...\n\n")
    }
    mod <- R.utils::doCall("spike_and_slab_logistic_moretrees", 
                           dsgn = dsgn,
                           args = args)
    if (log_restarts) {
      cat("\nRestart", i, "complete.")
      sink()
    }
    mod
  }
  
  # Select random restart that gave the highest ELBO
  ELBO_restarts <- sapply(mod_restarts, FUN = function(mod) mod$ELBO_track[length(mod$ELBO_track)])
  best_restart <- which.max(ELBO_restarts)
  mod <- mod_restarts[[best_restart]]
  if (keep_restarts) {
    mod_restarts <- mod_restarts[- best_restart]
  } else {
    rm(mod_restarts)
    mod_restarts <- NULL
  }
  
  # Compute MOReTreeS exposure coefficient estimates from model output
  betas <- moretrees_compute_betas(mod = mod, ci_level = ci_level,
                                   outcomes = outcomes,
                                   A_leaves = dsgn$A[names(igraph::V(tr))[igraph::V(tr)$leaf], ])
  
  # Compute MOReTreeS covariate coefficient estimates from model output
  if (!is.null(W)) {
    theta_est <- moretrees_compute_thetas(mod = mod, ci_level = ci_level, 
                                          m = ncol(W), A_leaves = dsgn$A[names(igraph::V(tr))[igraph::V(tr)$leaf], ])
  } else {
    theta_est <- NULL
  }
  
  # Get maximum likelihood estimates by group for comparison
  if (get_ml) {
    beta_ml <- ml_by_group(X = X, W = W, y = y, outcomes = outcomes,
                           outcome_groups = betas$beta_moretrees$outcomes,
                           ci_level = ci_level,
                           family = "binomial")
  } else {
    beta_ml <- NULL
  }
  
  # Return results
  return(list(beta_est = betas$beta_est,
              beta_moretrees = betas$beta_moretrees,
              beta_ml = beta_ml, 
              theta_est = theta_est,
              mod = mod,
              mod_restarts = mod_restarts))
}