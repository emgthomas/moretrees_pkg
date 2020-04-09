# --------------------------------------------------------------------------------- #
# ----- computing MOReTreeS estimates for covariates effects from model output ---- #
# --------------------------------------------------------------------------------- #

#' \code{moretrees_compute_thetas} Computes MOReTreeS estimates and credible
#'   intervals for covariate effects from model output
#' 
#' @param mod List containing outputs from spike and slab VI algorithm
#' @param ci_level A number between 0 and 1 giving the desired credible interval.
#' For example, ci_level = 0.95 (the default) returns a 95\% credible interval.
#' @param m Integer number of variables in covariate matrix W
#' @param A_leaves pL x p sparse ancestor Matrix where rows correspond to leaves
#' of tree (outcomes) and columns correspond to nodes on tree. Results in mod
#' must have same ordering as columns of A_leaves.
#' @param outcomes Character vector of length n. outcomes[i] is a string indicating the 
#' outcome experienced by unit i.
#' @return A matrix containing the estimates and confidence intervals.
#' @family Processing model output

moretrees_compute_thetas <- function(mod, ci_level,
                                     m, A_leaves = NULL) {
  
  # Get thetas from xis
  pL <- nrow(A_leaves)
  p <- ncol(A_leaves)
  theta_est <- matrix(nrow = pL, ncol = m)
  for (j in 1:m) {
    theta_est[ , j] <- as.numeric(A_leaves %*%
             sapply(mod$vi_params$delta, function(delta) delta[j , ]))
  }
  theta_names <- sapply(1:m, function(i) paste0("est",i))
  colnames(theta_est) <- theta_names
  
  # Compute credible intervals
  theta_sd_est <- matrix(nrow = pL, ncol = m)
  for (j in 1:m) {
    theta_sd_est[ , j] <- as.numeric(A_leaves %*%
         sapply(mod$vi_params$Omega, function(Omega) diag(Omega)[j]))
  }
  z <- qnorm(ci_level + (1 - ci_level) / 2)
  theta_ci_l <- theta_est - z * theta_sd_est
  colnames(theta_ci_l) <- sapply(1:m, function(i) paste0("cil",i))
  theta_ci_u <- theta_est + z * theta_sd_est
  colnames(theta_ci_u) <- sapply(1:m, function(i) paste0("ciu",i))
  theta_est <- cbind(theta_est, theta_ci_l, theta_ci_u)
  rownames(theta_est) <- rownames(A_leaves)
  
  # re-order columns for readability
  cols <- c("est", "cil", "ciu")
  cols <- sapply(1:m, function(i) paste0(cols, i), simplify = T) %>%
    as.vector
  theta_est <- theta_est[ , cols]
  
  # Return results
  return(theta_est)
  
}