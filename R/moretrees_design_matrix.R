# --------------------------------------------------------------------------------- #
# -------- Code for converting design matrix + tree into MOReTreeS ---------------- #
# -------- Design Matrix ---------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

# X is an n x K matrix
# tr is an igraph tree, where the leaves represent outcomes
# outcomes is a character vector of length n, where entry i
# tells us which outcome is represented by unit i

moretrees_design_matrix <- function(X, tr, outcomes) {
  # Check
  if (!is.character(outcomes)) stop("outcomes is not a character object")
  if (!is.igraph(tr)) stop("tr is not a graph object")
  if (!is.directed(tr)) stop("tr is not a directed graph object")
  
  nodes <- names(V(tr))
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  if(!setequal(unique(outcomes), leaves)) stop("Not all outcomes are leaves of tree")
  
  # Extract relevant parameters
  p <- length(nodes)
  pL <- length(leaves)
  K <- ncol(X)
  n <- nrow(X)
  A <- igraph::as_adjacency_matrix(tr, sparse = T)
  A <- expm(Matrix::t(A))
  A[A > 0 ] <- 1 
  
  # Get list of MOReTreeS design matrices for each node
  Xstar <- rep(list(Matrix::Matrix(0, nrow = n, ncol = K)), p)
  names(Xstar) <- nodes
  for (k in 1:K) {
    # Get design matrix for variable k
    Xsplt_k <- sapply(leaves, function(v) X[outcomes == v, k], simplify = F)
    Xmat_k <- Matrix(0, nrow = n, ncol = p)
    Xmat_k[ , nodes %in% leaves] <- Matrix::bdiag(Xsplt_k)
    rm(Xsplt_k)
    Xstar_k <- Xmat_k %*% A
    rm(Xmat_k)
    # Split design matrix among nodes as appropriate
    for (v in 1:v) {
      Xstar[[v]][ , k] <- Xstar_k[ , v]
    }
    rm(Xstar_k)
  }
  
  return(list(Xstar = Xstar, A = A))
}