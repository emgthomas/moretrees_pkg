# --------------------------------------------------------------------------------- #
# ------------------------- moretrees wrapper function ---------------------------- #
# --------------------------------------------------------------------------------- #

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
#' @param X An n x K matrix of exposure data, where K is the dimension of the exposure.
#' Grouping of the outcomes will be based on their relationships with the variables in X.
#' @param W Matrix of covariates of dimension n x m.
#' Coefficients for these variables do not affect grouping of the outcomes.
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
#' 
# X is an n x K matrix
# y is a vector of outcomes (either continuous or binary)
# outcomes is a character vector of length n, where entry i
# tells us which outcome is represented by unit i
# tr is an igraph tree, where the leaves represent outcomes

moretrees <- function(X, y, outcomes, tr,
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
  
  # Get design matrix
  dsgn <- moretrees_design_matrix(X = X, y = y, outcomes = outcomes, tr = tr)
  
  # Run algorithm
  mod <- ss_fun(dsgn$y, dsgn$Xstar, W,
                update_hyper = update_hyper, 
                update_hyper_freq = update_hyper_freq,
                tol = tol,
                max_iter = max_iter,
                hyperparams_init = hyperparams_init)
  
  # Get betas from xis
  p <- length(mod$vi_params$prob)
  xi_est <- t(sapply(1:p, 
                   function(v) as.numeric(mod$vi_params$mu[[v]] * (mod$vi_params$prob[v] >= 0.5)),
                   simplify = T))
  beta_est <- as.matrix(dsgn$A[names(V(tr))[V(tr)$leaf], ] %*% xi_est) %>%
    as.data.frame
  K <- ncol(beta_est)
  beta_names <- sapply(1:K, function(i) paste0("est",i))
  colnames(beta_est) <- beta_names
  
  # Compute credible intervals
  xi_var_est <- sapply(1:p, 
                         function(v) diag(as.matrix(mod$vi_params$Sigma[[v]])) * 
                                                 (mod$vi_params$prob[v] >= 0.5),
                         simplify = T) %>% t
  beta_sd_est <- as.matrix(dsgn$A[names(V(tr))[V(tr)$leaf], ] %*% xi_var_est) %>%
                             sqrt
  z <- qnorm(ci_level + (1 - ci_level) / 2)
  beta_ci_l <- beta_est - z * beta_sd_est
  names(beta_ci_l) <- sapply(1:K, function(i) paste0("cil",i))
  beta_ci_u <- beta_est + z * beta_sd_est
  names(beta_ci_u) <- sapply(1:K, function(i) paste0("ciu",i))
  beta_est <- cbind(beta_est, beta_ci_l, beta_ci_u)
  beta_est$group <- as.numeric(beta_est$est1) %>%
    as.factor %>% as.integer

  # Get estimated coefficients and CIs by group
  beta_moretrees <- beta_est[!duplicated(beta_est), ]
  beta_moretrees <- beta_moretrees[order(beta_moretrees$group), ]
  row.names(beta_moretrees) <- NULL
  beta_moretrees$outcomes <- vector(mode = "list", length = nrow(beta_moretrees))
  G <- nrow(beta_moretrees)
  for (i in 1:G) {
    beta_moretrees$outcomes[[i]] <- row.names(beta_est)[beta_est$group == i]
  }
  # re-order columns for readability
  cols <- c("est", "cil", "ciu")
  cols <- sapply(1:K, function(i) paste0(cols, i), simplify = T) %>%
    as.vector
  beta_moretrees <- beta_moretrees[ , c("group", cols, "outcomes")]
  
  # Get maximum likelihood estimates by group for comparison
  if (get_ml) {
    beta_ml <- matrix(nrow = G, ncol = K * 3 + 1) %>%
      as.data.frame
    names(beta_ml) <- c("group", cols)
    beta_ml$group <- 1:G
    beta_ml$outcomes <- beta_moretrees$outcomes
    if (family == "bernoulli") family <- "binomial"
    for (g in 1:G) {
      which_i <- outcomes %in% beta_moretrees$outcomes[[g]]
      mod_ml <- glm(y[which_i] ~ 0 + X[which_i, ], 
                    family = family)
      suppressMessages(beta_ml_ci <- confint(mod_ml, level = ci_level))
      beta_ml[g, paste0("est", 1:K)] <- mod_ml$coefficients
      beta_ml[g, paste0("cil", 1:K)] <- beta_ml_ci[ , 1]
      beta_ml[g, paste0("ciu", 1:K)] <- beta_ml_ci[ , 2]
    }
  } else {
    beta_ml <- NULL
  }
  
  # Return results
  return(list(beta_est = beta_est,
              beta_moretrees = beta_moretrees,
              beta_ml = beta_ml, 
              mod = mod))
}