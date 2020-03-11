% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moretrees_design_tree.R
\name{moretrees_design_tree}
\alias{moretrees_design_tree}
\title{Here's a brief description.
  \code{moretrees_design_tree} converts outcome, exposure, and covariate data
  into format suitable for analysis using MOReTreeS.}
\usage{
moretrees_design_tree(y, X, W = NULL, outcomes, tr)
}
\arguments{
\item{y}{Vector of length n containing outcomes data.
If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
If family = "gaussian", y must be a numeric vector containing continuous data.}

\item{X}{An n x K matrix of exposure data, where K is the dimension of the exposure.
Grouping of the outcomes will be based on their relationships with the variables in X.}

\item{W}{Matrix of covariates of dimension n x m.
Coefficients for these variables do not affect grouping of the outcomes.
Default is NULL (no covariates).}

\item{outcomes}{is a character vector of length n, where entry i}

\item{tr}{is an igraph tree, where the leaves represent outcomes}
}
\value{
A list containing the following elements:
y: Re-ordered outcome vector.
X: Re-ordered exposure matrix.
W: Re-ordered covariate matrix.
outcomes_units: list of length equal to the number of unique outcomes. Each element of
the list is an integer vector indicating which units (entries of y_reord, rows of X_reord)
correspond to each outcomes.
#' outcomes_nodes: list of length equal to the number of unique nodes. Each element of
the list is an integer vector indicating which outcomes are descendants of each node.
ancestors: list of length equal to the number of unique outcomes. Each element of the 
list is an integer vector indicating which nodes on the tree (including leaves) are ancestors
of the corresponding outcome.
}
\description{
All the details go here!
}
\section{Model Description}{

  Describe group spike and slab prior and all parameters here.
}

\seealso{
Other MOReTreeS functions: \code{\link{moretrees_compute_betas}},
  \code{\link{moretrees_compute_thetas}},
  \code{\link{moretrees_design_matrix}},
  \code{\link{moretrees_init_W_logistic}},
  \code{\link{moretrees_init_logistic}},
  \code{\link{moretrees_init_rand}},
  \code{\link{moretrees}}
}
\concept{MOReTreeS functions}