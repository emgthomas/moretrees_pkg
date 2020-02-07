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
  
  p <- length(nodes)
  pL <- length(leaves)
  K <- ncol(X)
  n <- nrow(X)
  A <- igraph::as_adjacency_matrix(tr, sparse = T)
  A <- expm(Matrix::t(A))
  A[A > 0 ] <- 1 
  
  Xstar <- list()
  for (k in 1:K) {
    X_splt_k <- sapply(leaves, function(v) X[outcomes == v, k], simplify = F)
    Xmat_k <- Matrix(0, nrow = n, ncol = p)
    Xmat_k[ , nodes %in% leaves] <- Matrix::bdiag(X_splt_k)
    colnames(Xmat_k) <- nodes
    Xstar[[k]] <- Xmat_k %*% A
  }
  
  # Generate data
  X_list <- rep(list(), p)
  for (v in 1:p) {
    X_splt <- sapply(leaves, function(v) X[outcomes == v, ], simplify = F)
    X_splt <- Matrix::bdiag(X_splt)
    Xmat_v <- Matrix::Matrix(0, nrow = n, ncol = K)
    for (k in 1:K) {
      Xmat_v[ , k] <- Xstar[[k]][ , v]
    }
    X_list[[v]] <- Xmat_v
  }
  
  return(list(X = X_list, A = A))
}