#' Long-form printing of MOReTreeS model fits
#' 
#' \code{print.summary.moretrees_long} is a print method for class 
#' \code{summary.moretrees_long}.
#' 
#' @export
#' @param x output from \code{summary.moretrees_result} with \code{compact = FALSE}.
#' @param print_outcomes If \code{TRUE}, for each discovered group the full list
#' of outcomes will be printed. Set this to \code{FALSE} if these outcome lists
#' make output difficult to read.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#' @return Summary showing, for each outcome group discovered by MOReTreeS,
#' the exposure coefficients, 95\% credible intervals, number of outcomes
#' per group, and number of matched pairs.
#' @examples 
#' # See vignette
#' vignette("moretrees")
#' @family MOReTrees results

print.summary.moretrees_long <- function(x,
                                         print_outcomes = T,
                                         digits = max(3L, getOption("digits") - 3L),
                                         ...) {
  
  transform <- x$transform
  x$transform <- NULL
  
  cat("Group-specific exposure effect estimates for the", length(x), "groups discovered by MOReTreeS\n")
  if (x$coeff_type == "clr") {
    cat("Showing conditional logistic regression estimates for discovered groups.\n\n")
  }
  x$coeff_type <- NULL
  
  for (g in 1:length(x)) {
    
    est <- x[[g]]
    cat("------------------ Group", g, "------------------\n\n")
    
    cat("Number of outcomes:", est$n_outcomes, "\n")
    cat("Number of matched pairs:", est$n_obs, "\n")
    if (print_outcomes) {
      cat("List of outcomes:\n")
      cat(est$outcomes, "\n\n")
    } else {
      cat("\n\n")
    }
    
    if (transform == "exp") cat("Odds ratio estimate(s):")
    if (transform == "er") cat("Excess rate estimate(s):")
    if (transform == "identity") cat("Coefficient estimate(s):")
    print(knitr::kable(est$est, digits = digits))
    
    cat("\n")
  }
  
  cat("\nIf this is hard to read, try print(x, compact = TRUE)\n\n")
  
  # Return
  return(invisible(x))
}
