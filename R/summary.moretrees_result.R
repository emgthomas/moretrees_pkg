
#' \code{summary.moretrees_result} summarizes the results from \code{moretrees()}.
#' 
#' @export
#' @param object Output from \code{moretrees()}. 
#' An object of class "moretrees_result".
#' @param transform Either "exp", "identity", or "er" (excess rate)
#' @param coeff_type Either "moretrees" or "clr"
#' @param compact If \code{TRUE}, a more compact summary of results is printed.
#' Only works well when the dimension of the exposure variable is low
#' (say, < 4); otherwise the table takes up too much horizontal space.
#' Default is \code{FALSE}.
#' @param ... Not used.
#' @return Summary showing, for each outcome group discovered by MOReTreeS,
#' the exposure coefficients, 95\% credible intervals, number of outcomes
#' per group, and number of matched pairs.
#' @examples 
#' # See vignette
#' vignette("moretrees")
#' @family MOReTrees results

summary.moretrees_result <- function(object,
                                     transform = "exp",
                                     coeff_type = "moretrees",
                                     compact = F,
                                     ...) {
  
  # Get estimates
  k <- length(object$mod$vi_params$mu[[1]])
  cols <- unlist(lapply(1:k, function(k) paste0(c("est", "cil", "ciu"), k)))
  if (coeff_type == "clr") {
    coeff_type2 <- "ml"
  } else {
    coeff_type2 <- coeff_type
  }
  beta_type <- paste0("beta_", coeff_type2)
  est <- object[beta_type][[1]][ , cols]
  if (transform %in% c("exp", "et")) est <- exp(est)
  if (transform == "er") est <- (est - 1) * 100
  est$n_outcomes <- sapply(object$beta_moretrees$outcomes, length)
  est$outcomes <- sapply(object$beta_moretrees$outcomes, paste0, collapse = ", ")
  est$n_obs <- object$beta_moretrees$n_obs
  est$group <- 1:nrow(est)
  
  if (compact) {
    rslt <- est[ , c("group","n_outcomes", "n_obs", 
                     paste0(c("est", "cil", "ciu"), rep(1:k, each = 3)),
                     "outcomes")]
    rslt <- list(est = rslt, transform = transform,
                 coeff_type = coeff_type)
    class(rslt) <- "summary.moretrees_compact"
    return(rslt)
  }
  
  # make separate data.frames per group
  grps <- list()
  for (g in 1:max(est$group)) {
    est_g <- est[est$group == g, ]
    grps[[g]] <- list(n_outcomes = est_g$n_outcomes,
                      n_obs = est_g$n_obs,
                      outcomes = est_g$outcomes)
    grps[[g]]$est <- reshape(est_g[ , cols],
                             varying = list(paste0("est", 1:k),
                                            paste0("cil", 1:k),
                                            paste0("ciu", 1:k)),
                             v.names = c("est", "cil", "ciu"),
                             timevar = "dim",
                             direction = "long")
    grps[[g]]$est$id <- NULL
    
  }
  
  names(grps) <- paste0("Group", 1:length(grps))
  grps$transform <- transform
  grps$coeff_type <- coeff_type
  
  class(grps) <- "summary.moretrees_long"
  
  # Return
  return(grps)
}
