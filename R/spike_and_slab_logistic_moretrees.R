#' Description
#'
#'   \code{spike_and_slab_logistic_moretrees()} fits MOReTreeS models for
#'   binary data.
#'   The posterior is approximated via variational inference.
#'   This function returns the parameters of the variational approximation.
#'   
#' @param dsgn list of outputs from \code{moretrees_design_tree()} 
#' @param vi_params_init,hyperparams_init,random_init,random_init_vals,tol,tol_hyper,max_iter,print_freq,update_hyper_freq,hyper_fixed
#' see documentation for \code{moretrees()}
#' @return A named list containing the following entried:
#' \describe{
#' \item{\code{vi_params}}{named list of final variational parameter estimates}
#' \item{\code{hyperparams}}{named list of final hyperparameter estimates}
#' \item{\code{hyper_fixed}}{named list of fixed hyperparameters}
#' \item{\code{ELBO_track}}{numeric vector containing the values of the objective function
#' (ELBO) at the end of every iteration}
#' }
#' @family Internal VI functions

spike_and_slab_logistic_moretrees <- function(dsgn, 
                                              vi_params_init,
                                              hyperparams_init,
                                              random_init,
                                              random_init_vals,
                                              tol,
                                              tol_hyper,
                                              max_iter,
                                              print_freq,
                                              update_hyper_freq,
                                              hyper_fixed) {
  
  # Add some data elements to dsgn ---------------------------------------------------
  if (is.null(dsgn$W)) {
    dsgn$W <- matrix(nrow = length(dsgn$y), ncol = 0)
  }
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
  init <- R.utils::doCall(moretrees_init_logistic, 
                          vi_params = vi_params_init,
                          hyperparams = hyperparams_init,
                          hyper_fixed = hyper_fixed,
                          random_init = random_init,
                          random_init_vals = random_init_vals,
                          args = dsgn)
  vi_params <- init$vi_params
  hyperparams <- init$hyperparams
  
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
      warning("Maximum number of iterations reached! Consider increasing max_iter.")
      break
    }
    
    # update vi params
    vi_params <- R.utils::doCall(update_vi_params_logistic_moretrees, 
                                 args = c(dsgn, vi_params, hyperparams, hyper_fixed))
    
    # compute ELBO and update eta
    update_hyper <- i %% update_hyper_freq == 0
    hyperparams <-  R.utils::doCall(update_hyperparams_logistic_moretrees, 
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