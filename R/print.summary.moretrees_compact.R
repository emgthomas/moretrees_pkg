#' Compact printing of MOReTreeS model fits
#' 
#' \code{print.summary.moretrees_compact} is a print method for class 
#' \code{summary.moretrees_compact}.
#' 
#' @export
#' @param x output from \code{summary.moretrees_result} with \code{compact = TRUE}.
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

print.summary.moretrees_compact <- function(x,
                                            print_outcomes = T,
                                            digits = max(3L, getOption("digits") - 3L),
                                            ...) {
  
  transform <- x$transform
  if (!print_outcomes) x$est$outcomes <- NULL
  
  if (x$coeff_type == "clr") {
    cat("Showing conditional logistic regression estimates for discovered groups.\n\n")
  }
  
  if (transform == "exp") cat("Group-specific odds ratio estimate(s):")
  if (transform == "er") cat("Group-specific excess rate estimate(s):")
  if (transform == "identity") cat("Group-specific coefficient estimate(s):")
  print(knitr::kable(x$est, digits = digits))
  
  cat("\n")
  
  # Return
  return(invisible(x))
}
