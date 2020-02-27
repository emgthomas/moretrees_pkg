### Sparse Matrix Multiplication Optimization ###

rm(list = ls())
devtools::load_all()

# 1. Generating test data -----------------------------------------------------------------------

# Sample size
n <- 1E6 # NOTE: ultimately I would like this to work for sample size about 20 * 1E6

# Input parameters -------------------------------------------------------------------
# No need to pay too much attention to the following
{family <- "bernoulli"
group <- "7"
tr <- ccs_tree(group)$tr
leaves <- names(igraph::V(tr)[igraph::V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- Matrix::expm(Matrix::t(A))
A[A > 0 ] <- 1
G <- length(igraph::V(tr))
p <- G
pL <- sum(igraph::V(tr)$leaf)
K_g <- 3 # number of variables
K <- rep(K_g, G)
m <- 2
tau <- 3
rho1 <- 0.6 # rho for internal nodes
rho2 <- 0.05 # rho for leaf nodes
rho <- sum(1 + 0.8 * (p - pL - 1) + 0.05 * pL) / p # overall rho
omega <- 2
sigma2 <- 2

# Generate randomly grouped beta (groups follow tree)
set.seed(4563457)
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
s_true <- c(1, rbinom(n = p - pL - 1, size = 1, prob = rho1), 
            rbinom(n = pL, size = 1, prob = rho2))
xi <- t(mapply(function(gamma, s) matrix(gamma * s, nrow = 1),
             gamma = gamma_true, s = s_true,
             SIMPLIFY = T))
if (K_g == 1) xi <- t(xi)
beta <- A[leaves, ] %*% xi
zeta <- matrix(rnorm(m * G, mean = 0, sd = sqrt(omega)), nrow = G, ncol = m)
theta <- A[leaves, ] %*% zeta

# Generate data --------------------------------------------------------------------------
X <- Matrix::Matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- sample(leaves, size = n, replace = T)
if (m > 0) {
  W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n, ncol = m)
} else {
  W <- NULL
}

# Get linear predictor
lp <- numeric(n)
for (v in leaves) {
  which_v <- outcomes == v
  lp[which_v] <- lp[which_v] + X[which_v, , drop = F] %*% Matrix::t(beta[v, , drop = F])
  if (m > 0) {
    lp[which_v] <- lp[which_v] + W[which_v, , drop = F] %*% Matrix::t(theta[v, , drop = F])
  }
}

# Simulate outcomes
p_success <- expit(lp)
y <- runif(n)
y <- as.integer(y <= p_success)
}

# Get large sparse matrix --------------------------------------------------------------------
# next line is a bit slow
dsgn <- moretrees_design_matrix(y, X, W, outcomes, tr)

# These are the sparse matrices I'm working with. You'll see they do have a certain structure.
# It's not block diagonal as I originally said 
X <- dsgn$Xstar
W <- dsgn$Wstar
groups <- dsgn$groups

# Prepare for running algorithm ---------------------------------------------------
G <- length(groups)
n <- length(y)
K <- sapply(groups, length)
m <- ncol(W)
# Initial hyperparameter values
eta <- abs(rnorm(n, mean = 0, sd = 10))
g_eta <- gfun(eta)
tau <- 5
omega <- 5
# Initial VI values
A_eta <- Matrix::Diagonal(n = n, g_eta)
Sigma_inv <- lapply(X = groups, 
                    FUN = function(cols, Xg, A, tau) 2 * 
                      Matrix::crossprod(Xg[ , cols, drop = F], A) %*% Xg[ , cols, drop = F] + 
                      Matrix::Diagonal(length(cols), 1 / tau),
                    Xg = X,
                    tau = tau,
                    A = A_eta)
Sigma <- lapply(Sigma_inv, Matrix::solve)
Sigma_det <- sapply(Sigma, Matrix::det)
mu <- sapply(K, rnorm, mean = 0 , sd = 10, simplify = F)
mu <- sapply(mu, Matrix::Matrix, ncol = 1)
prob <- runif(G, 0 , 1)
tau_t <- rep(tau, G)
delta <- Matrix::Matrix(rnorm(m, sd = 10), ncol = 1)
Omega_inv <- 2 * Matrix::t(W) %*% A_eta %*% W + 
  Matrix::Diagonal(m, 1 / 10)
if (m != 0) {
  Omega <- Matrix::solve(Omega_inv)
  Omega_det <- Matrix::det(Omega)
} else {
  Omega <- Matrix::Matrix(nrow = 0, ncol = 0)
  Omega_det <- 1
}

require(gdata)
keep(W, delta, X, Sigma, Omega, groups, xi, y, g_eta, prob, mu, sure = T)
save(W, delta, X, Sigma, Omega, groups, xi, y, g_eta, prob, mu, file = "/Users/emt380/Documents/PhD_Papers/MOReTreeS_Second_paper/R_code/tinkering/test_data.Rdata")

# Lines ------------------------------------------------------------------------------------

rm(list = ls())
load("/Users/emt380/Documents/PhD_Papers/MOReTreeS_Second_paper/R_code/tinkering/test_data.Rdata")
Wdelta <- W %*% delta

# First bottleneck ------------------------------------------------
# This will be done repeatedly for all values of g in 1:length(Sigma)
# Must be done in correct order so the calculations for each g can't be parallelized
g <- 2
mu[[g]] <- Sigma[[g]] %*% Matrix::crossprod( X[ , groups[[g]], drop = F],
        y / 2 - 2 * g_eta * (Wdelta + X[ , -groups[[g]], drop = F] %*% xi[-groups[[g]]]))

# Second bottleneck -----------------------------------------------
require(RcppEigen)
quadFormByRowCpp <- "using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::VectorXd;
const MappedSparseMatrix<double> X(as<MappedSparseMatrix<double> >(XX));
const MappedSparseMatrix<double> S(as<MappedSparseMatrix<double> >(SS));
return wrap(((X * S).cwiseProduct(X)) * VectorXd::Ones(X.cols()));"

quadFormByRow <- inline::cxxfunction(signature(XX = "dgCMatrix", SS = "dgCmatrix"), 
                                     quadFormByRowCpp, plugin="RcppEigen")


lp2 <- quadFormByRow(SS = Omega, XX = W) # this line is slowest
for (g in 1:length(Sigma)) { # this for loop could be parallelized
  Sigma_g <- Sigma[[g]] + (1 - prob[g]) * Matrix::tcrossprod(mu[[g]])
  lp2 <- lp2 + prob[g] * quadFormByRow(SS = as(Sigma_g, "dgCMatrix"), 
       XX = X[, groups[[g]], drop = F]) # this line is also somewhat slow
}




