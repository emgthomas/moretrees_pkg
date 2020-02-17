##### Generate some test data #####
devtools::load_all() # Sources all files in R/

# Input parameters -------------------------------------------------------------------
group <- "7"
tr <- ccs_tree(group)$tr
leaves <- names(igraph::V(tr)[igraph::V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- Matrix::expm(Matrix::t(A))
A[A > 0 ] <- 1
G <- length(igraph::V(tr))
n <- 1000
K_g <- 2 # number of variables
K <- rep(K_g, G)
m <- 2
tau <- 3
rho <- 0.5
omega <- 2

# Generate some data -----------------------------------------------------------------
X <- Matrix::Matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
if (m > 0) {
  W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n, ncol = m)
} else {
  W <- NULL
}

dsgn <- moretrees_design_matrix(y = rep(1, n), X = X, W = W, outcomes = outcomes, tr = tr)
X <- dsgn$Xstar
W <- dsgn$Wstar
groups <- dsgn$groups

# Generate fake variational params ---------------------------------------------------
G <- length(groups)
n <- nrow(X)
K <- sapply(groups, length)
m <- ncol(W)
eta <- abs(rnorm(n, mean = 0, sd = 2))
g_eta <- gfun(eta)
A_eta <- Matrix::Diagonal(n = n, g_eta)
Sigma_inv <- sapply(X = groups, 
                    FUN = function(cols, Xg, A, tau) 2 * 
                      Matrix::crossprod(Xg[ , cols, drop = F], A) %*% Xg[ , cols, drop = F] + 
                      Matrix::Diagonal(length(cols), 1 / tau),
                    Xg = X,
                    tau = tau,
                    A = A_eta)
Sigma <- sapply(Sigma_inv, Matrix::solve)
Sigma_det <- sapply(Sigma, Matrix::det)
mu <- sapply(K, rnorm, mean = 0 , sd = 10, simplify = F)
mu <- sapply(mu, Matrix::Matrix, ncol = 1)
prob <- runif(G, 0 , 1)
tau_t <- rep(tau, G)
delta <- Matrix::Matrix(rnorm(m, sd = 10), ncol = 1)
Omega_inv <- 2 * Matrix::t(W) %*% A_eta %*% W + 
  Matrix::Diagonal(m, 1 / omega)
if (m != 0) {
  Omega <- Matrix::solve(Omega_inv)
  Omega_det <- Matrix::det(Omega)
} else {
  Omega <- Matrix::Matrix(nrow = 0, ncol = 0)
  Omega_det <- 1
}

# Computing expected linear predictor squared -----------------------------------------------------

# Without using Matrix::crossprod
f0 <- function(X) {
  lp2 <- numeric(n) + 0
  for (g in 1:G) {
    lp2 <- lp2 + prob[g] * apply(X[ , groups[[g]], drop = F], MARGIN = 1,
                                 FUN = function(x, Sigma) as.numeric(Matrix::t(x) %*% Sigma %*% x),
                                 Sigma = Sigma[[g]] + (1 - prob[g]) * mu[[g]] %*% Matrix::t(mu[[g]]))
  }
  lp2
}

# Current approach
f1 <- function(X) {
  lp2 <- numeric(n) + 0
  for (g in 1:G) {
    lp2 <- lp2 + prob[g] * apply(X[ , groups[[g]], drop = F], MARGIN = 1,
                                 FUN = function(x, Sigma) as.numeric(Matrix::crossprod(x, Sigma) %*% x),
                                 Sigma = Sigma[[g]] + (1 - prob[g]) * Matrix::tcrossprod(mu[[g]]))
  }
  lp2
}


# Defining function outside apply
xT_Sigma_x <- function(x, Sigma) {
  as.numeric(Matrix::crossprod(x, Sigma) %*% x)
}
f2 <- function(X) {
  lp2 <- numeric(n) + 0
  for (g in 1:G) {
    lp2 <- lp2 + prob[g] * apply(X[ , groups[[g]], drop = F], MARGIN = 1,
                                 FUN = xT_Sigma_x,
                                 Sigma = Sigma[[g]] + (1 - prob[g]) * Matrix::tcrossprod(mu[[g]]))
  }
  lp2
}

# Benchmarking ---------------------------------------------------------------

# Compare f1 & f0 
require(microbenchmark)
microbenchmark(f1(X), f0(X), times = 10)
# f1 slightly better

# Compare f1 & f2
microbenchmark(f1(X), f2(X), times = 100)
# f1 still faster! weird?

# Compare f1 and f2 on sparse matrices
microbenchmark(f1(X), f2(X), times = 10)
# f2 is slight winner!
