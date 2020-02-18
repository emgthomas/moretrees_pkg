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
n <- 1E6
K_g <- 2 # number of variables
K <- rep(K_g, G)
m <- 2
tau <- 3
rho <- 0.5
omega <- 2

# Generate some data -----------------------------------------------------------------
X0 <- Matrix::Matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
if (m > 0) {
  W0 <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n, ncol = m)
} else {
  W0 <- NULL
}

dsgn <- moretrees_design_matrix(y = rep(1, n), X = X0, W = W0, outcomes = outcomes, tr = tr)
X <- dsgn$Xstar
W <- dsgn$Wstar
groups <- dsgn$groups
X1 <- as.matrix(X0[dsgn$ord, ])
W1 <- as.matrix(W0[dsgn$ord, ])
outcomes1 <- outcomes[dsgn$ord]
nodes <- colnames(dsgn$A)

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
delta2 <- Matrix::Matrix(as.numeric(delta), nrow = G, byrow = F)
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


# Defining function outside apply
xT_Sigma_x <- function(x, Sigma) {
  as.numeric(Matrix::crossprod(x, Sigma) %*% x)
}
f2 <- function(X) {
  lp2 <- numeric(n) + 0
  for (g in 1:G) {
    lp2 <- lp2 + apply(X[ , groups[[g]], drop = F], MARGIN = 1,
                       FUN = xT_Sigma_x,
                       Sigma = prob[g] * (Sigma[[g]] + (1 - prob[g]) * Matrix::tcrossprod(mu[[g]])))
  }
  lp2
}

# Current approach
f1 <- function(X) {
  lp2 <- numeric(n) + 0
  for (g in 1:G) {
    lp2 <- lp2 + apply(X[ , groups[[g]], drop = F], MARGIN = 1,
                       FUN = function(x, Sigma) as.numeric(Matrix::crossprod(x, Sigma) %*% x),
                       Sigma = prob[g] * (Sigma[[g]] + (1 - prob[g]) * Matrix::tcrossprod(mu[[g]])))
  }
  lp2
}

# Doing it the moretrees way
ancestors <- igraph::ego(tr, order = 100, nodes = leaves, mode = "in")
ancestors <- sapply(ancestors, names)
ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a), nodes = nodes)
names(ancestors) <- leaves
outcomes_list <- sapply(leaves, function(v) which(outcomes1 == v), simplify = F)
names(outcomes_list) <- leaves
Sigma2 <- sapply(Sigma, as.matrix, simplify = F)
mu2 <- sapply(mu, as.matrix, simplify = F)
f3 <- function(X1) {
  lp2 <- numeric(n) + 0
  Sigma_u <- mapply(FUN = function(prob, Sigma, mu) prob * 
                      (Sigma + (1 - prob) * tcrossprod(mu)),
                    prob = prob, Sigma = Sigma2, mu = mu2,
                    SIMPLIFY = F)
  for (v in 1:length(ancestors)) {
    Sigma_v <- Reduce(`+`, Sigma_u[ancestors[[v]]])
    lp2[outcomes_list[[v]]] <- apply(X1[outcomes_list[[v]], ], 1,
                                     FUN = function(x, Sigma) crossprod(x, Sigma) %*% x,
                                     Sigma = Sigma_v)
  }
  lp2
}

# Benchmarking ---------------------------------------------------------------

# Compare sparse to non-sparse matrix 
require(microbenchmark)
X2 <- Matrix::Matrix(X, sparse = F)
microbenchmark(f0(X), f0(X2), times = 1)
# f1 slightly better

# Compare f0 to f1
microbenchmark(f0(X), f1(X), times = 10)
# f1 slightly better

# Compare f1 & f2
microbenchmark(f1(X), f2(X), times = 10)
# f1 still faster! weird?

# Compare f1 & f2
lp2_1 <- f1(X)
lp2_3 <- f3(X1)
all.equal(lp2_1, lp2_3)
microbenchmark(f1(X), f3(X1), times = 1)
# f3 (moretrees way) is about 10 times faster

X2 <- Matrix::Matrix(X1, sparse = F)
W2 <- Matrix::Matrix(W1, sparse = F)
lp2_3 <- f3(X2)
lp2_4 <- f3(X1)
microbenchmark(f3(X1), f3(X2), times = 10)
# moretrees way is slightly faster with non-sparse matrix

# Doing it the moretrees way for expected linear predictor -------------------------------

# Expected linear predictor
f1 <- function(X, W) {
  lp <- W %*% delta
  for (g in 1:G) {
    lp <- lp + prob[g] *  X[ , groups[[g]], drop = F] %*% mu[[g]]
  }
  as.numeric(lp)
}

# Expected linear predictor moretrees way
f2 <- function(X1, W1) {
  xi_u <- mapply(FUN = function(prob, mu) prob * mu,
                 prob = prob, mu = mu2, SIMPLIFY = F)
  lp <- numeric(n) + 0
  for (v in 1:length(ancestors)) {
    beta_v <- Reduce(`+`, xi_u[ancestors[[v]]])
    theta_v <- apply(delta2[ancestors[[v]], ], 2, sum)
    lp[outcomes_list[[v]]] <- X1[outcomes_list[[v]], ] %*% beta_v + W1[outcomes_list[[v]], ] %*% theta_v
  }
  lp
}

# Expected linear predictor


# benchmark
microbenchmark(f1(X, W), f2(X1, W1), times = 10)
# moretrees way is about 30 times faster

f2(X1, W1)
f1(X, W)

X2 <- Matrix::Matrix(X1, sparse = F)
W2 <- Matrix::Matrix(W1, sparse = F)
microbenchmark(f2(X1, W1), f2(X2, W2), times = 10)
# about 20 times faster with non-sparse matrix (but this won't have a huge overall impact)


