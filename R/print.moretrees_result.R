
#' \code{print.moretrees_result} summarizes the results from \code{moretrees()}.
#' 
#' @export
#' @param x Output from \code{moretrees()}.
#' @param ... Arguments passed to summary and printing methods.
#' @return Summary showing, for each outcome group discovered by MOReTreeS,
#' the exposure coefficients, 95\% credible intervals, number of outcomes
#' per group, and number of matched pairs.
#' @examples 
#' # See vignette
#' vignette("moretrees")
#' @family MOReTrees results

print.moretrees_result <- function(x, ...) {
  
  print(summary(x, ...), ...)
  
  # Return
  return(invisible(x))
}
