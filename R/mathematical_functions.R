# --------------------------------------------------------------------------------- #
# ---------- Useful mathematical functions not implemented in base R  ------------- #
# --------------------------------------------------------------------------------- #

trace_prod <- function(A, B) {
  # computes trace(A %*% B) for square matrices A, B
  sum(t(B) * A)
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

loglogit <- function(x) {
  # computes -log(1+exp(-x)) for x a vector
  -log1p.exp.vec(-x)
}

logfac <- function(x) {
  if(x == 0) return(0)
  sum(log(1:x))
}