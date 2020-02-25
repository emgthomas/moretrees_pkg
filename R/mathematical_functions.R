# --------------------------------------------------------------------------------- #
# ---------- Useful mathematical functions not implemented in base R  ------------- #
# --------------------------------------------------------------------------------- #

trace_prod <- function(A, B) {
  # computes trace(A %*% B) for square matrices A, B
  sum(Matrix::t(B) * A)
}

log1p_exp <- function(x) {
  # computes log(1 + exp(x))
  if (x > 20) {
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

log1p.exp.vec <- function(x){
  # computes log(1 + exp(x)) for x a vector
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}

logexpit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

expit <- function(x) {
  # computes 1/(1+exp(-x)) for x a vector
  exp(logexpit(x))
}

gfun <- function(x) {
  # computes function g described in manuscript;
  # needed for normal approx to logistic likelihood
  (expit(x) - 1 / 2) / (2 * x)
}

# quadFormByRow <- function(Sigma, X) Matrix::rowSums(Matrix::tcrossprod(X, Sigma) * X)
