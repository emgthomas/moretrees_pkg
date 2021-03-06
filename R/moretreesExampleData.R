#' Simulated multi-outcome matched case-control data.
#'
#' Simulated data for testing the \code{moretrees()} function.
#' Used in the package vignette. The dataset contains
#' 1000 matched case-control pairs. The outcomes
#' are cardiovascular diseases in category 7.4 
#' (diseases of arteries, arterioles, and capillaries) 
#' in the multilevel Clinical Classifications Software (CCS)
#' hierarchical disease classification system.
#'
#' @format A named list with the following elements:
#' \describe{
#'   \item{Xcase}{A matrix with 1000 rows and 2 columns.
#'   \code{Xcase[i, ]} represents exposure data for the case in 
#'   matched-pair \code{i}.}
#'   \item{Xcontrol}{A matrix with 1000 rows and 2 columns.
#'   Xcontrol[i, ] represents exposure data for the control in 
#'   matched-pair i.}
#'   \item{Wcase}{A matrix with 1000 rows and 2 columns.
#'   \code{Wcase[i, ]} represents covariate data for the case in 
#'   matched-pair \code{i}.}
#'   \item{Wcontrol}{A matrix with 1000 rows and 2 columns.
#'   Wcontrol[i, ] represents covariate data for the control in 
#'   matched-pair i.}
#'   \item{outcomes}{A character vector of length 1000.
#'   \code{outcomes[i]} indicates which outcome was experienced by
#'   the case in matched pair \code{i}. The entries are CCS codes
#'   belonging to category 7.4.}
#' }
"moretreesExampleData"