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

spike_and_slab_logistic_moretrees <- function(dsgn, initial_values,
                                              random_init,
                                              tol, max_iter,
                                              update_hyper, 
                                              update_hyper_freq,
                                              hyper_fixed,
                                              print_freq,
                                              hyper_random_init,
                                              vi_random_init) {
  if (is.null(dsgn$W)) {
    dsgn$W <- matrix(nrow = length(dsgn$y), ncol = 0)
  }
  
  # Prepare for running algorithm ---------------------------------------------------
  dsgn$n <- length(dsgn$y)
  dsgn$m <- ncol(dsgn$W)
  dsgn$p <- length(unique(unlist(dsgn$ancestors)))
  dsgn$pL <- length(dsgn$ancestors)
  dsgn$K <- ncol(dsgn$X)
  dsgn$Fg <- max(dsgn$levels)
  if (dsgn$K == 1) {
    dsgn$xxT <- dsgn$X ^ 2
  } else {
    dsgn$xxT <- rowOuterProds(dsgn$X)
  }
  if (dsgn$m > 0) {
    if (dsgn$m == 1) {
      dsgn$wwT <- dsgn$W ^ 2
    } else {
      dsgn$wwT <- rowOuterProds(dsgn$W)
    }
  } else {
    dsgn$wwT <- NULL
  }
  # Initial values
  if (is.null(initial_values)) {
    if (random_init) {
      initial_values <- R.utils::doCall(moretrees_init_rand, vi_random_init = vi_random_init,
                                        hyper_random_init = hyper_random_init,
                                        hyper_fixed = hyper_fixed,
                                        args = dsgn)
    } else {
      initial_values <- R.utils::doCall(moretrees_init_logistic, 
                                        hyper_fixed = hyper_fixed,
                                        args = dsgn)
    }
  }
  
  # else {
  #   # In some cases, we may supply starting values for mu but not delta
  #   if (m > 0 & nrow(initial_values$vi_params$delta[[1]]) != m) {
  #     initial_values <- R.utils::doCall(moretrees_init_W_logistic,
  #                                        initial_values = initial_values,
  #                                        hyper_fixed = hyper_fixed,
  #                                        args = dsgn)
  #   } # otherwise, all initial values should already be supplied
  # }
  vi_params <- initial_values$vi_params
  hyperparams <- initial_values$hyperparams
  hyper_fixed <- initial_values$hyper_fixed
  
  # Do first updates so that ELBO calculation will be correct
  # (some cancellation in formula relies on recent updates)
  vi_params <- R.utils::doCall(update_vi_params_logistic_moretrees, 
                               args = c(dsgn, vi_params, hyperparams, hyper_fixed))
  hyperparams <-  R.utils::doCall(update_hyperparams_logistic_moretrees,
                                  args = c(dsgn, vi_params, hyperparams, hyper_fixed))
  
  # Initialise ELBO
  ELBO_track <- numeric(max_iter %/% update_hyper_freq)
  
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    # check if max_iter reached
    if (i > max_iter) {
      cat(paste("Iteration", i, "complete.\n"))
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
  
    # iterate i
    i <- i + 1
    
    # update vi params
    vi_params <- R.utils::doCall(update_vi_params_logistic_moretrees, 
                      args = c(dsgn, vi_params, hyperparams, hyper_fixed))
    
    # compute ELBO and update eta
    hyperparams <-  R.utils::doCall(update_hyperparams_logistic_moretrees,
                      args = c(dsgn, vi_params, hyperparams, hyper_fixed))
    ELBO_track[i] <- hyperparams$ELBO
    
    # check tolerance
    if (i > 2 && abs(ELBO_track[i] - ELBO_track[i - 1]) < tol) break
    
    # print progress
    if (i %% print_freq == 0 & i > 2) {
      cat("Iteration", i, "; epsilon =", ELBO_track[i] - ELBO_track[i - 1] , "\n")
    }
  }
  
  # return results
  return(list(vi_params = vi_params, 
              hyperparams = hyperparams,
              hyper_fixed = hyper_fixed,
              ELBO_track = ELBO_track[1:i]))
}