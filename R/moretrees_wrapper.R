# --------------------------------------------------------------------------------- #
# ------------------------- moretrees wrapper function ---------------------------- #
# --------------------------------------------------------------------------------- #

#' Description
#' 
#' @import Rcpp
#' @import stats
#' @export
#' @useDynLib moretrees
#' 
#' @param Xcase An \code{n} x \code{K} matrix of exposure data for cases, where \code{K} is the dimension of the exposure.
#' Grouping of the outcomes is based on their associations with variables in \code{Xcase}.
#' Rows of \code{Xcase} correspond to inividual cases, columns correspond to variables.
#' @param Xcontrol An \code{n} x \code{K} matrix of exposure data for controls; row \code{i} in \code{Xcontrol} is the matched
#' control for case \code{i}.
#' @param Wcase An \code{n} x \code{m} matrix of covariate data for cases, where m is the dimension of the exposure.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' Rows of \code{Wcase} correspond to inividual cases, columns correspond to variables.
#' @param Wcontrol An \code{n} x \code{m} matrix of covariate data for controls; row \code{i} in \code{Wcontrol} is the matched
#' control for case \code{i}.
#' @param outcomes Character vector of length \code{n}. \code{outcomes[i]} is a string indicating the 
#' outcome experienced by unit \code{i}.
#' @param tr A directed \code{igraph} object. This is a tree representing the relationships
#' among the outcomes. The leaves represent individual outcomes, and internal nodes
#' represent outcome categories consisting of their leaf descendants. All nodes
#' of tr must have unique names as given by \code{names(V(tr))}. The names of the leaves must 
#' be equal to the unique elements of outcomes. The vertices of \code{tr}, \code{V(tr)}, may have 
#' an attribute \code{levels} containing integer values from 1 to \code{max(V(tr)$levels)}. 
#' In this case, the levels attribute specifies groups of nodes that share common 
#' hyperparameters \code{rho[f]}, \code{tau[f]}, and \code{omega[f]}. If \code{V(tr)$levels} is \code{NULL}, 
#' the default is two levels of hyperparameters: one for all leaf nodes, and one
#' for all internal nodes.
#' @param ci_level A number between 0 and 1 giving the desired credible interval. 
#' For example, \code{ci_level = 0.95} (the default) returns a 95\% credible interval
#' @param get_ml If \code{TRUE}, moretrees will also return the maximum likelihood estimates of the
#' coefficients for each outcome group discovered by the model. 
#' Default is \code{TRUE}.
#' @param update_hyper_freq How frequently to update hyperparameters. 
#' Default = every 50 iterations.
#' @param print_freq How often to print out iteration number and current value of epsilon
#' (the difference in objective function value for the two most recent iterations). 
#' @param hyper_fixed Fixed values of hyperprior parameters for rho.
#' This should be a list with two elements: 
#' a and b, both numeric vectors of length \code{L}, representing the 
#' parameters of the beta prior on rho for each level, where \code{L} is the
#' number of levels.
#' Default is \code{list(a = rep(1, L), b = rep(1, L))} (uniform hyperprior)
#' @param tol Convergence tolerance for the objective function.
#' Default is \code{1E-8}.
#' @param tol_hyper The convergence tolerance for the objective function between
#' between subsequent hyperparmeter updates. Typically a more generous
#' tolerance than \code{tol}. 
#' Default is \code{1E-4}.
#' @param max_iter Maximum number of iterations of the VI algorithm.
#' Default is 5000.
#' @param nrestarts Number of random re-starts of the VI algorithm. The result that 
#' gives the highest value of the objective function will be returned.
#' It is recommended to choose \code{nrestarts > 1}.
#' The default is 3.
#' @param keep_restarts If \code{TRUE}, the results from all random restarts will be returned.
#' If \code{FALSE}, only the restart with the highest objective function is returned. '
#' Default is \code{TRUE}.
#' @param parallel If \code{TRUE}, the random restarts will be run in parallel.
#' It is recommended to first set the number of cores using \code{doParallel::registerDoParallel()}. 
#' Otherwise, the default number of cores specified by the \code{doParallel} package will be used.
#' Default is \code{TRUE}.
#' @param log_restarts If \code{TRUE}, when \code{nrestarts > 1} progress of each random restart will be 
#' logged to a text file in \code{log_dir}. If \code{FALSE} and \code{nrestarts > 1}, 
#' progress will not be shown.
#' If \code{nrestarts = 1}, progress will always be printed to the console.
#' Default is \code{FALSE}.
#' @param log_dir Directory for logging progress of random restarts.
#' Default is the working directory.
#' @param vi_params_init,hyperparams_init Named lists containing initial values for the 
#' variational parameters and hyperparameters. Supplying good initial values can be challenging,
#' and \code{moretrees()} provides a way to guess initial values based on transformations
#' of conditional logistic regression estimates of the effect sizes
#' for each individual outcome (see \code{moretrees_init_logistic()}).
#' The most common use for \code{vi_params_init} and \code{hyperparams_init} is to supply starting
#' values based on previous output from \code{moretrees()}; 
#' see the \code{vignette('moretrees')} for examples.
#' The user can provide initial values for all parameters or a subset. 
#' When initial values for one or more parameters are not
#' supplied, the missing values will be filled in by \code{moretrees_init_logistic()}.
#' @param random_init 
#' If \code{TRUE}, some random variability will be added to the initial values.
#' The default is \code{FALSE}, unless \code{nrestarts > 1}, in which case 
#' \code{random_init} will be set to \code{TRUE} and a warning message will be printed.
#' The amount of variability is determined by \code{random_init_vals}.
#' @param random_init_vals If \code{random_init = TRUE}, 
#' this is a list containing the following parameters for randomly permuting 
#' the inital values:
#' \describe{
#' \item{\code{tau_lims}}{a vector of length 2, where \code{tau_lims[1]} is between 0 and 1,
#' and \code{tau_lims[2] > 1}. The initial values for the hyperparameter \code{tau} will
#' be chosen uniformly at random in the range \code{(tau_init * tau_lims[1], tau_init * tau_lims[2])},
#' where \code{tau_init} is the initial value for \code{tau} either supplied in \code{hyperparams_init} 
#' or guessed using \code{moretrees_init_logistic()}.}
#' \item{\code{omega_lims}}{a vector of length 2, where \code{omega_lims[1]} is between 0 and 1,
#' and \code{omega_lims[2] > 1}. The initial values for the hyperparameter omega will
#' be chosen uniformly at random in the range \code{(omega_init * omega_lims[1], omega_init * omega_lims[2])},
#' where omega_init is the initial value for omega either supplied in \code{hyperparams_init} or guessed
#' using \code{moretrees_init_logistic()}.}
#' \item{\code{eta_sd_frac}}{a value between 0 and 1. The initial values for the auxilliary parameters
#' \code{eta} will have a normal random variate added to them with standard deviation equal to 
#' \code{eta_sd_frac} multiplied by the initial value for eta either supplied in \code{hyperparams_init} or guessed
#' using \code{moretrees_init_logistic()}. Absolute values are then taken for any 
#' values of \code{eta} that are \code{< 0}.}
#' \item{\code{mu_sd_frac}}{a value between 0 and 1. The initial values for
#' \code{mu} will have a normal random variate added to them with standard deviation equal to 
#' \code{mu_sd_frac} multiplied by the absolute value of the initial value for \code{mu} either supplied in 
#' \code{vi_params_init} or guessed using \code{moretrees_init_logistic()}.}
#' \item{\code{delta_sd_frac}}{a value between 0 and 1. The initial values for
#' \code{delta} will have a normal random variate added to them with standard deviation equal to 
#' \code{delta_sd_frac} multiplied by the absolute value of the initial value for delta either supplied in 
#' \code{vi_params_init} or guessed using \code{moretrees_init_logistic()}.}
#' \item{\code{u_sd_frac}}{a value between 0 and 1. The initial value for the node inclusion probabilities
#' will first be transformed to the log odds scale to obtain \code{u}. A normal random variate will be
#' added to \code{u} with standard deviation eqaul to u_sd_frac multiplied by the absolute value of the
#' initial value for \code{u} either supplied in \code{vi_params_init} or guessed using \code{moretrees_init_logistic()}.
#' \code{u} will then be transformed back to the probability scale.}
#' }
#' @return A list containing the following elements:
#' \describe{
#' \item{\code{beta_est}}{estimated exposure coefficients and credible intervals for each outcome.
#' This is a data frame where the variables \code{est1, cil1, ciu1} correspond to the estimated
#' coefficient and lower and upper credible interval bounds for the variable in first column
#' of \code{Xcase}/\code{Xcontrol}. \code{est2, cil2, ciu2}, correspond to the second column in 
#' \code{Xcase}/\code{Xcontrol},
#' and so on. The variable group indicates to which estimated group each outcome belongs. } 
#' \item{\code{beta_moretrees}}{estimated exposure coefficients and credible intervals for each outcome group.
#' This is the same information in beta_est but presented by group.  \code{Outcomes} is a list of the
#' outcomes in each group and \code{n_obs} is the number of 
#' matched pairs corresponding to those outcomes.} 
#' \item{\code{theta_est}}{estimated covariate coefficients and credible intervals for each outcome.
#' This is a matrix where the columns \code{est1, cil1, ciu1} correspond to the estimated
#' coefficient and lower and upper credible interval bounds for the variable in first column
#' of \code{Wcase}/\code{Wcontrol}. \code{est2, cil2, ciu2}, correspond to the second column in 
#' \code{Wcase}/\code{Wcontrol}, and so on.}
#' \item{\code{beta_ml, theta_ml}}{Results from running separate, classic conditional logisitic 
#' regression models on the data from observations corresponding to each outcome group
#' shown in \code{beta_moretrees}.}
#' \item{\code{mod}}{outputs from variational inference algorithm}
#' \item{\code{mod_restarts}}{outputs from other random restarts of the algorithm, if
#' \code{nrestarts > 1}}
#' }
#' @examples vignette('moretrees')
#' @family MOReTreeS functions

moretrees <- function(Xcase, Xcontrol, 
                      Wcase = NULL, Wcontrol = NULL,
                      outcomes, 
                      tr,
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
                      vi_params_init = list(),
                      hyperparams_init = list(),
                      random_init = FALSE,
                      random_init_vals = list(omega_lims = c(0.5, 1.5),
                                              tau_lims = c(0.5, 1.5),
                                              eta_sd_frac = 0.2,
                                              mu_sd_frac = 0.2,
                                              delta_sd_frac = 0.2,
                                              u_sd_frac = 0.2)) {
  
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
    mod <- spike_and_slab_logistic_moretrees(dsgn = dsgn,
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
    ml <- ml_by_group(X = X, W = W, y = y, outcomes = outcomes,
                      outcome_groups = betas$beta_moretrees$outcomes,
                      ci_level = ci_level,
                      family = "binomial",
                      return_theta = TRUE)
    beta_ml <- ml$beta_ml
    theta_ml <- ml$theta_ml
  } else {
    beta_ml <- NULL
    theta_ml <- NULL
  }
  
  # Return results
  return(list(beta_est = betas$beta_est,
              beta_moretrees = betas$beta_moretrees,
              beta_ml = beta_ml, 
              theta_est = theta_est,
              theta_ml = theta_ml,
              mod = mod,
              mod_restarts = mod_restarts))
}