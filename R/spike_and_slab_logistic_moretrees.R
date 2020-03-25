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
#' @param dsgn A list containing all data elements needed for the algorithm.
#' @param tol Convergence tolerance for ELBO.
#' @param maxiter Maximum number of iterations of the VI algorithm.
#' @param print_freq How often to print out iteration number.
#' @param hyper_fixed Fixed values of hyperprior parameters.
#' @param update_hyper_freq How frequently to update hyperparameters.
#' @return A list of variational parameters.
#' @examples
#' @family spike and slab functions

spike_and_slab_logistic_moretrees <- function(dsgn, 
                                              initial_values,
                                              random_init,
                                              tol,
                                              tol_hyper,
                                              max_iter,
                                              print_freq,
                                              update_hyper_freq,
                                              hyper_fixed,
                                              hyper_random_init,
                                              vi_random_init) {
  
  if (is.null(dsgn$W)) {
    dsgn$W <- matrix(nrow = length(dsgn$y), ncol = 0)
  }
  
  # Add some data elements to dsgn ---------------------------------------------------
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
  
  # Get initial values ------------------------------------------------------------
  if (is.null(initial_values)) {
    if (random_init) {
      initial_values <- R.utils::doCall(moretrees_init_rand, 
                                        vi_random_init = vi_random_init,
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
                               hyper_method = hyper_method,
                               args = c(dsgn, vi_params, hyperparams, hyper_fixed))
  hyperparams <-  R.utils::doCall(update_hyperparams_logistic_moretrees,
                                  hyper_method = hyper_method,
                                  update_hyper = F,
                                  args = c(dsgn, vi_params, hyperparams, hyper_fixed))
  
  # Initialise ELBO
  ELBO_track <- numeric(max_iter)
  
  # Run algorithm -----------------------------------------------------------------
  i <- 0
  repeat {
    
    # iterate i
    i <- i + 1
    
    # check if max_iter reached
    if (i > max_iter) {
      i <- max_iter
      cat(paste("Iteration", i, "complete.\n"))
      cat("\nWarning: Maximum number of iterations reached!\n")
      break
    }
    
    # update vi params
    vi_params <- R.utils::doCall(update_vi_params_logistic_moretrees, 
                                 hyper_method = hyper_method,
                                 args = c(dsgn, vi_params, hyperparams, hyper_fixed))
    
    # compute ELBO and update eta
    update_hyper <- i %% update_hyper_freq == 0
    hyperparams <-  R.utils::doCall(update_hyperparams_logistic_moretrees, 
                                    hyper_method = hyper_method,
                                    update_hyper = update_hyper,
                                    args = c(dsgn, vi_params, hyperparams, hyper_fixed))
    ELBO_track[i] <- hyperparams$ELBO
    
    # print progress
    if (i %% print_freq == 0 & i > 3) {
      cat("Iteration", i, "; epsilon =", ELBO_track[i] - ELBO_track[i - 1] , "\n")
    }
    
    # check tolerance
    if (update_hyper & i >= 2 * update_hyper_freq) {
      # if we just updated hyperparameters, check for convergence of hyperparameters
      criterion1 <- abs(ELBO_track[i] - ELBO_track[i - update_hyper_freq]) < tol_hyper
      if (criterion1) {
        # did last VI update reach convergence?
        criterion2 <- abs(ELBO_track[i - 1] - ELBO_track[i - 2]) < tol
        # if yes, both have converged. if not, continue.
        if (criterion2) break else next
      } else next
    } else {
      criterion3 <- (i > 2) && (abs(ELBO_track[i] - ELBO_track[i - 1]) < tol)
      # if criterion 3, fill in results until just before the 
      # next hyperparameter update (or max_iter, whichever comes first)
      if (criterion3) {
        i2 <- min(ceiling(i / update_hyper_freq) * update_hyper_freq - 1, 
                  max_iter)
        ELBO_track[(i + 1):i2] <- hyperparams$ELBO 
        i <- i2
      }
    }
    
  } # close repeat loop
  
  # return results
  return(list(vi_params = vi_params, 
              hyperparams = hyperparams,
              hyper_fixed = hyper_fixed,
              ELBO_track = ELBO_track[1:i]))
}