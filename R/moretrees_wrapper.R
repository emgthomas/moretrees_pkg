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
#' @param outcomes Character vector specifying which outcomes the elements of y/rows of 
#' X and W correspond to. 
#' @param W_method = "shared" if information about the effect of variables in W wil be shared
#' across the outcomes according to the tree structure. If W_method = "individual", the effect of
#' W will be estimated separately for each outcome (no infromation sharing).
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95% credible interval.
#' @param get_ml If TRUE, moretrees will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model. The default is FALSE.
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

moretrees <- function(X, W = NULL, y, outcomes, tr,
                      W_method = "shared",
                      family = "bernoulli",
                      ci_level = 0.95,
                      get_ml = FALSE,
                      update_hyper = T, update_hyper_freq = 50,
                      tol = 1E-8, max_iter = 5000,
                      hyperparams_init = NULL) {

  if (!(family %in% c("bernoulli", "gaussian"))) {
    stop("family must be a string (\"bernoulli\" or \"gaussian\")")
  }
  if (family == "bernoulli") ss_fun <- spike_and_slab_logistic
  if (family == "gaussian") ss_fun <- spike_and_slab_normal
  if (!(length(get_ml) == 1 & is.logical(get_ml))) stop("get_ml must be either TRUE or FALSE")
  
  # Get MOReTreeS design matrices
  dsgn <- moretrees_design_matrix(X = X, W = W, y = y,
                                  outcomes = outcomes, tr = tr,
                                  W_method = W_method)

  # Run algorithm
  mod <- ss_fun(y = dsgn$y_reord, X = dsgn$Xstar, W = dsgn$Wstar,
                update_hyper = update_hyper, 
                update_hyper_freq = update_hyper_freq,
                tol = tol,
                max_iter = max_iter,
                hyperparams_init = hyperparams_init)
  
  # Compute MOReTreeS exposure coefficient estimates from model output
  betas <- moretrees_compute_betas(mod = mod, ci_level = ci_level,
                                   A_leaves = dsgn$A[names(V(tr))[V(tr)$leaf], ])
  
  # Compute MOReTreeS covariate coefficient estimates from model output
  theta_est <- moretrees_compute_thetas(mod = mod, ci_level = ci_level, 
                                        m = ncol(W),
                                        A_leaves = dsgn$A[names(V(tr))[V(tr)$leaf], ],
                                        W_method = W_method)
  
  # Get maximum likelihood estimates by group for comparison
  if (get_ml) {
    beta_ml <- ml_by_group(X = X, W = W, y = y, outcomes = outcomes,
                           outcome_groups = betas$beta_moretrees$outcomes,
                           ci_level = ci_level,
                           family = family)
  } else {
    beta_ml <- NULL
  }
  
  # Return results
  return(list(beta_est = betas$beta_est,
              beta_moretrees = betas$beta_moretrees,
              beta_ml = beta_ml, 
              theta_est = theta_est,
              mod = mod))
}