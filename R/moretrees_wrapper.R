# --------------------------------------------------------------------------------- #
# ------------------------- moretrees wrapper function ---------------------------- #
# --------------------------------------------------------------------------------- #

#' Here's a brief description.
#'   \code{moretrees} Fits Multi-Outcome Regression with Tree-structured Shrinkage
#'   (MOReTreeS) model to matched case-control or case-crossover data.
#'   The posterior is approximated via variational inference.
#'   Returns estimated outcome groups and group-specific coefficient 
#'   estimates with credible intervals.
#' 
#' All the details go here!
#' 
#' @export
#' @useDynLib moretrees
#' 
#' @section Model Description:
#' Describe MOReTreeS model and all parameters here.
#' 
#' @param Xcase An n x K matrix of exposure data for cases, where K is the dimension of the exposure.
#' Grouping of the outcomes is based on their associations with variables in Xcase.
#' Rows of Xcase correspond to inividual cases, columns correspond to variables.
#' @param Xcontrol An n x K matrix of exposure data for controls; row i in Xcontrol is the matched
#' control for case i.
#' @param Wcase An n x m matrix of covariate data for cases, where m is the dimension of the exposure.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' Rows of Wcase correspond to inividual cases, columns correspond to variables.
#' @param Wcontrol An n x m matrix of covariate data for controls; row i in Wcontrol is the matched
#' control for case i.
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
#' the default is two levels of hyperparameters: one for all leaf nodes, and one
#' for all internal nodes.
#' @param ci_level A number between 0 and 1 giving the desired credible interval. 
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval
#' @param get_ml If TRUE, moretrees will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model. 
#' Default is TRUE.
#' @param tol Convergence tolerance for ELBO. 
#' Default is 1E-8.
#' @param tol_hyper If hyper_method = "EB", the convergence tolerance for ELBO
#' between subsequent hyperparmeter updates. Typically a more generous
#' tolerance than tol. 
#' Default is 1E-4.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' Default is 5000.
#' @param hyper_fixed Fixed values of hyperprior parameters for rho.
#' This should be a list with two elements: 
#' a and b, both numeric vectors of length L, representing the 
#' parameters of the beta prior on rho for each level, where L is the
#' number of levels.
#' Default is list(a = rep(1, L), b = rep(1, L)).
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param random_init The initial values for the MOReTreeS model are selected based
#' on maximum likelihood estimates of the effect size at every level of the tree.
#' If random_init = TRUE, some randomness will be added to these initial values.
#' The default is FALSE, unless nrestarts > 1, in which case random_init will be set
#' to TRUE and a warning message will be printed.
#' @param random_init_vals If random_init = TRUE, 
#' this is a list containing the following parameters for randomly permuting 
#' the inital values:
#' 1. tau_lims: a vector of length 2, where tau_lims[1] is between 0 and 1,
#' and tau_lims[2] > 1. The initial values for the hyperparameter tau will
#' be chosen uniformly at random in the range (tau_init \* tau_lims[1], tau_init \* tau_lims[2]),
#' where tau_init is the initial value for tau either supplied in initial_values or guessed
#' using moretrees_init_logistic().
#' 1. omega_lims: a vector of length 2, where omega_lims[1] is between 0 and 1,
#' and omega_lims[2] > 1. The initial values for the hyperparameter omega will
#' be chosen uniformly at random in the range (omega_init \* omega_lims[1], omega_init \* omega_lims[2]),
#' where omega_init is the initial value for omega either supplied in initial_values or guessed
#' using moretrees_init_logistic().
#' 1. eta_sd_frac: a value between 0 and 1. The initial values for the auxilliary parameters
#' eta will have a normal random variate added to them with standard deviation equal to 
#' eta_sd_frac multiplied by the initial value for eta either supplied in initial_values or guessed
#' using moretrees_init_logistic(). Absolute values are then taken for any 
#' values of eta that are < 0.
#' 1. mu_sd_frac: a value between 0 and 1. The initial values for the exposure effects
#' mu will have a normal random variate added to them with standard deviation equal to 
#' mu_sd_frac multiplied by the absolute value of the initial value for mu either supplied in 
#' initial_values or guessed using moretrees_init_logistic().
#' 1. delta_sd_frac: a value between 0 and 1. The initial values for the exposure effects
#' delta will have a normal random variate added to them with standard deviation equal to 
#' delta_sd_frac multiplied by the absolute value of the initial value for delta either supplied in 
#' initial_values or guessed using moretrees_init_logistic().
#' @param print_freq How often to print out iteration number and current value of epsilon
#' (the difference in objective function value for the two most recent iterations). 
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that 
#' gives the highest ELBO will be returned. It is recommended to choose nrestarts > 1.
#' The default is 3.
#' @param keep_restarts If TRUE, the results from all random restarts will be returned.
#' If FALSE, only the restart with the highest ELBO is returned. Default is TRUE.
#' @param parallel If TRUE, the random restarts will be run in parallel.
#' It is recommended to first set the number of cores using doParallel::registerDoParallel(). 
#' Otherwise, the default number of cores specified by the doParallel package will be used.
#' Default is TRUE.
#' @param log_restarts If TRUE, when nrestarts > 1 progress of each random restart will be 
#' logged to a text file in log_dir. If FALSE and nrestarts > 1, progress will not be shown.
#' If nrestarts = 1, progress will always be printed to the console.
#' Default is FALSE.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples 
#' @family MOReTreeS functions

moretrees <- function(Xcase, Xcontrol, 
                      Wcase = NULL, Wcontrol = NULL,
                      outcomes, 
                      tr,
                      vi_params_init = list(),
                      hyperparams_init = list(),
                      ci_level = 0.95,
                      get_ml = TRUE,
                      update_hyper_freq = 50,
                      print_freq = 50,
                      hyper_fixed = NULL,
                      tol = 1E-8, 
                      tol_hyper = 1E-4,
                      max_iter = 5000,
                      nrestarts = 3,
                      keep_restarts = TRUE,
                      parallel = TRUE,
                      log_restarts = FALSE,
                      log_dir = ".",
                      random_init = FALSE,
                      random_init_vals = list(omega_lims = c(0.5, 1.5),
                                              tau_lims = c(0.5, 1.5),
                                              eta_sd_frac = 0.2,
                                              mu_sd_frac = 0.2,
                                              delta_sd_frac = 0.2)) {
  
  if (!(is.matrix(Xcase) & is.matrix(Xcontrol))) stop("Xcase and Xcontrol must be matrices")
  if (!identical(dim(Wcase), dim(Wcontrol))) stop("Xcase and Xcontrol must have same dimension")
  if (!is.null(Wcase)) {
    if (!(is.matrix(Wcase) & is.matrix(Wcontrol))) stop("If not NULL, Wcase & Wcontrol must be matrices")
    if (!identical(dim(Wcase), dim(Wcontrol))) stop("If not NULL, Wcase and Wcontrol must have same dimension")
  }
  if (!(length(get_ml) == 1 & is.logical(get_ml))) stop("get_ml must be either TRUE or FALSE")
  log_dir <- sub("/$", "", log_dir)
  if (log_restarts) message("Algorithm progress for restart i will be printed to ",
                        log_dir, "/restart_i_log.txt\n", sep = "")
  
  # Fill in some arguments
  if (nrestarts > 1 & !random_init) {
    message("Setting random_init = TRUE since nrestarts > 1\n")
    random_init <- TRUE
  }
  if (nrestarts == 1) parallel <- FALSE
  
  # Get MOReTreeS design elements
  X <- Xcase - Xcontrol
  if (!is.null(Wcase)) {
    W <- Wcase - Wcontrol
  } else {
    W <- NULL
  }
  y <- rep(1, nrow(X))
  dsgn <- moretrees_design_tree(X = X, W = W, y = y, outcomes = outcomes, tr = tr)
  
  # Get hyper_fixed if not supplied
  if (is.null(hyper_fixed)) {
    L <- max(dsgn$levels)
    hyper_fixed <- list(a = rep(1, L), b = rep(1, L))
  }
  
  # Setting up parallelization
  if (parallel) {
    `%doRestarts%` <- foreach::`%dopar%`
  } else {
    `%doRestarts%` <- foreach::`%do%`
  }
  
  # Run algorithm
  mod_restarts <- foreach::foreach(i = 1:nrestarts) %doRestarts% {
    if (log_restarts) {
      sink(file = paste0(log_dir, "/restart_", i, "_log.txt"))
    }
    cat("\nInitialising restart", i, "...\n\n")
    mod <- moretrees:::spike_and_slab_logistic_moretrees(dsgn = dsgn,
                                             vi_params_init = vi_params_init,
                                             hyperparams_init = hyperparams_init,
                                             random_init = random_init,
                                             random_init_vals = random_init_vals,
                                             tol = tol,
                                             tol_hyper = tol_hyper,
                                             max_iter = max_iter,
                                             print_freq = print_freq,
                                             update_hyper_freq = update_hyper_freq,
                                             hyper_fixed = hyper_fixed)
    cat("\nRestart", i, "complete.\n")
    if (log_restarts) {
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
                                   A_leaves = dsgn$A_leaves)
  
  # Compute MOReTreeS covariate coefficient estimates from model output
  if (!is.null(W)) {
    theta_est <- moretrees_compute_thetas(mod = mod, ci_level = ci_level, 
                                          m = ncol(W), 
                                          dsgn$A_leaves)
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