# --------------------------------------------------------------------------------- #
# ------ computing maximum likelihood coefficient estimates by outcome group ------ #
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
#' @param outcomes Character vector specifying which outcomes the elements of y 
#' (rows of X and W) correspond to. 
#' @param outcome_groups A list of length G, where each element is a character vector
#' indiciating which outcomes belong to that group.
#' @param family A string specifying the distribution of the outcomes: 
#' either "bernoulli" (for classification) or "gaussian" (for regression)
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95% credible interval.
#' @param return_theta Return ML estimates for theta? Default = FALSE
#' @return A list containing the following elements:
#' 1. estimated coefficients and credible intervals; 
#' 2. outputs from variational inference algorithm
#' @examples
#' @family spike and slab functions

ml_by_group <- function(X, W = NULL, y, outcomes, outcome_groups,
                        ci_level, family, return_theta = FALSE) {
  G <- length(outcome_groups)
  K <- ncol(X)
  beta_ml <- matrix(nrow = G, ncol = K * 3 + 1) %>%
    as.data.frame
  cols <- sapply(1:K, 
    function(i) paste0(c("est", "cil", "ciu"), i), simplify = T) %>%
    as.vector
  names(beta_ml) <- c("group", cols)
  beta_ml$group <- 1:G
  beta_ml$outcomes <- outcome_groups
  if (!is.null(W)) {
    m <- ncol(W)
    if (return_theta) {
      theta_ml <- matrix(nrow = G, ncol = m * 3 + 1) %>%
        as.data.frame
    }
    cols <- sapply(1:m, 
      function(i) paste0(c("est", "cil", "ciu"), i), simplify = T) %>%
      as.vector
    names(theta_ml) <- c("group", cols)
    theta_ml$group <- 1:G
    theta_ml$outcomes <- outcome_groups
  }
  for (g in 1:G) {
    which_i <- outcomes %in% outcome_groups[[g]]
    if (!is.null(W)) {
      mod_ml <- glm(y[which_i] ~ 0 + as.matrix(X[which_i, ]) + 
                      as.matrix(W[which_i, ]), 
                    family = family)
    } else {
      mod_ml <- glm(y[which_i] ~ 0 + as.matrix(X[which_i, ]), 
                    family = family)
    }
    suppressWarnings(suppressMessages(beta_ml_ci <- confint(mod_ml, level = ci_level)))
    if (K == 1 & is.null(W)) beta_ml_ci <- matrix(beta_ml_ci, nrow = K)
    beta_ml[g, paste0("est", 1:K)] <- mod_ml$coefficients[1:K]
    beta_ml[g, paste0("cil", 1:K)] <- beta_ml_ci[1:K , 1]
    beta_ml[g, paste0("ciu", 1:K)] <- beta_ml_ci[1:K , 2]
    if (return_theta) {
      theta_ml[g, paste0("est", 1:m)] <- mod_ml$coefficients[(K + 1):(K + m)]
      theta_ml[g, paste0("cil", 1:m)] <- beta_ml_ci[(K + 1):(K + m) , 1]
      theta_ml[g, paste0("ciu", 1:m)] <- beta_ml_ci[(K + 1):(K + m) , 2]
    }
  }
  if (return_theta) {
    return(list(beta_ml = beta_ml, theta_ml = theta_ml))
  } else {
    return(beta_ml)
  }
}
