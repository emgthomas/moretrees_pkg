# --------------------------------------------------------------------------------- #
# -------- Code for converting design matrix + tree into components --------------- #
# -------- needed to fit MOReTreeS model ------------------------------------------ #
# --------------------------------------------------------------------------------- #

#'   \code{moretrees_design_tree} converts outcome, exposure, and covariate data
#'   into format suitable for analysis using MOReTreeS.
#' 
#' 
#' @param y Vector of length n containing outcomes data.
#' If family = "bernoulli", y must be an integer vector where 1 = success, 0 = failure.
#' If family = "gaussian", y must be a numeric vector containing continuous data.
#' @param X An n x K matrix of exposure data, where K is the dimension of the exposure.
#' Grouping of the outcomes will be based on their relationships with the variables in X.
#' @param W Matrix of covariates of dimension n x m.
#' Coefficients for these variables do not affect grouping of the outcomes.
#' Default is NULL (no covariates).
#' @param outcomes is a character vector of length n, where entry i
#  tells us which outcome is represented by unit i
#' @param tr is an igraph tree, where the leaves represent outcomes.
#' The vertices of tr, V(tr), may have an attribute 'levels' containing integer
#' values from 1 to max(V(tr)$levels). In this case, the levels attribute specifies
#' groups of nodes that share common hyperparameters rho[f], tau[f], and omega[f].
#' If V(tr)$levels is NULL, moretrees_design_tree will specify levels based on
#' distance from the root node of the tree.
#' @return A list containing the following elements:
#' y: Re-ordered outcome vector.
#' X: Re-ordered exposure matrix.
#' W: Re-ordered covariate matrix.
#' A: Matrix for transforming node-level parameters to outcome-level parameters.
#' outcomes_units: list of length equal to the number of unique outcomes. Each element of
#' the list is an integer vector indicating which units (entries of y_reord, rows of X_reord)
#' correspond to each outcomes.
#' outcomes_nodes: list of length equal to the number of unique nodes. Each element of
#' the list is an integer vector indicating which outcomes are descendants of each node.
#' ancestors: list of length equal to the number of unique outcomes. Each element of the 
#' list is an integer vector indicating which nodes on the tree (including leaves) are ancestors
#' of the corresponding outcome.
#' levels: an integer vector of length equal to the number of unique nodes. Contains
#' integer values from 1 to max(levels), where the values indicate groups of nodes
#' that share common hyperparameters rho[f], tau[f], and omega[f], where f is the value
#' in levels.
#' @family Data processing functions

moretrees_design_tree <- function(y, X, W = NULL, outcomes, tr) {
  # Some checks
  if (!is.character(outcomes)) stop("outcomes is not a character object")
  if (!igraph::is.igraph(tr)) stop("tr is not a graph object")
  if (!igraph::is.directed(tr)) stop
  if (is.integer(y) & !(sum(y %in% c(0, 1)) == length(y))) 
    stop("y contains values other than zero or one")
  
  nodes <- names(igraph::V(tr))
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  if(!setequal(unique(outcomes), leaves)) {
    stop("Not all outcomes are leaves of tree")
  }
  # Re-order nodes to have internal nodes first, then leaves
  nodes <- c(nodes[!(nodes %in% leaves)], leaves)
  # Get levels for specifying groups of hyperparameters
  if (is.null(igraph::V(tr)$levels)) {
    # The default is to have two levels: one for leaf nodes, one for internal nodes
    levels <- rep(1, length(nodes))
    levels[nodes %in% leaves] <- 2
  } else {
    # Otherwise, use levels supplied
    levels <- igraph::V(tr)[nodes]$levels
    levels <- levels[match(names(igraph::V(tr)), nodes)]
  }
  if (sum(table(levels) < 5) > 0) {
    warning("Some levels contain fewer than five nodes: this may lead to problems
            with estimating variance parameters. Recommend increasing the number
            of nodes per level.")
  }

  # Extract relevant parameters
  p <- length(nodes)
  pL <- length(leaves)
  K <- ncol(X)
  n <- nrow(X)
  
  # Get ancestor matrix
  A <- igraph::as_adjacency_matrix(tr, sparse = T)
  A <- A[nodes, nodes] # re-order rows/columns to mirror nodes
  A <- Matrix::expm(Matrix::t(A))
  A[A > 0 ] <- 1 
  A <- Matrix::Matrix(A, sparse = T)
  
  # Sort by outcomes, where order is specified by ordering in tr
  ord <- order(ordered(outcomes, levels = leaves))
  X <- X[ord, , drop = F]
  X <- as.matrix(X)
  if (!is.null(W)){
    W <- as.matrix(W)
    W <- W[ord, , drop = F]
  } 
  y <- y[ord]
  outcomes <- outcomes[ord]
  
  # Replace y = 0 with y = -1 for compatibility with moretrees algorithm
  if (is.integer(y)) y[y == 0] <- -1
  
  # get lists of ancestors for each outcome
  d <- igraph::diameter(tr)
  ancestors <- igraph::ego(tr, order = d + 1, nodes = leaves, mode = "in")
  ancestors <- sapply(ancestors, names, simplify = F)
  ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a), nodes = nodes,
                      simplify = F)
  names(ancestors) <- leaves
  
  # get lists of which units correspond to each outcome
  outcomes_units <- sapply(leaves, function(v) which(outcomes == v), simplify = F)
  names(outcomes_units) <- leaves
  
  # get lists of outcomes are descendants of each node
  descendants <- igraph::ego(tr, order = d + 1, nodes = nodes, mode = "out")
  descendants <- sapply(descendants, names)
  outcomes_nodes <- sapply(descendants, function(d, leaves) which(leaves %in% d), leaves = leaves,
                           simplify = F)
  names(outcomes_nodes) <- nodes
  
  # Change outcomes to integers
  outcomes <- sapply(outcomes, function(o) which(leaves == o))
  
  # return
  return(list(y = y, 
              X = X, 
              W = W, 
              A = A,
              A_leaves = A[nodes %in% leaves, ],
              outcomes = outcomes,
              outcomes_units = outcomes_units, 
              outcomes_nodes = outcomes_nodes,
              ancestors = ancestors,
              levels = levels))
  
}