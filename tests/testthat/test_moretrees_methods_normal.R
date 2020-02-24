# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian and Bernoulli outcomes  ----------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

rm(list = ls())
devtools::load_all() # Sources all files in R/

# Chose one --------------------------------------------------------------------------
family <- "gaussian"
# family <- "bernoulli"

# Input parameters -------------------------------------------------------------------
group <- "7"
tr <- ccs_tree(group)$tr
leaves <- names(igraph::V(tr)[igraph::V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- Matrix::expm(Matrix::t(A))
A[A > 0 ] <- 1
G <- length(igraph::V(tr))
p <- G
pL <- sum(igraph::V(tr)$leaf)
n <- 500
K_g <- 2 # number of variables
K <- rep(K_g, G)
m <- 2
tau <- 3
rho1 <- 0.6 # rho for internal nodes
rho2 <- 0.05 # rho for leaf nodes
rho <- sum(1 + 0.8 * (p - pL - 1) + 0.05 * pL) / p # overall rho
omega <- 2
sigma2 <- 1

# Generate randomly grouped beta (groups follow tree)
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
s_true <- c(1, rbinom(n = p - pL - 1, size = 1, prob = rho1), 
            rbinom(n = pL, size = 1, prob = rho2))
xi <- mapply(function(gamma, s) matrix(gamma * s, nrow = 1),
             gamma = gamma_true, s = s_true,
             SIMPLIFY = T) %>% t
if (K_g == 1) xi <- t(xi)
beta <- as.matrix(A[leaves, ] %*% xi)
zeta <- matrix(rnorm(m * G, mean = 0, sd = sqrt(omega)), nrow = G, ncol = m)
theta <- as.matrix(A[leaves, ] %*% zeta)
groups_true <- as.integer(as.factor(as.numeric(beta[ , 1])))
table(groups_true)

# Generate some data -----------------------------------------------------------------
X <- matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
if (m > 0) {
  W <- matrix(rnorm(m * n, sd = 0.5), nrow = n, ncol = m)
} else {
  W <- NULL
}

# Get linear predictor
lp <- numeric(n)
for (v in leaves) {
  which_v <- outcomes == v
  lp[which_v] <- lp[which_v] + X[which_v, , drop = F] %*% t(beta[v, , drop = F])
  if (m > 0) {
    lp[which_v] <- lp[which_v] + W[which_v, , drop = F] %*% t(theta[v, , drop = F])
  }
}

# Simulate outcomes
if (family == "gaussian") {
  y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
} else {
  p_success <- expit(lp)
  y <- sapply(p_success, rbinom, n = 1, size = 1)
}

# Design matrices -------------------------------------------------------------------

dsgn_mat <- moretrees_design_matrix(y = y, X = X, W = W, outcomes = outcomes,
                                    tr = tr, W_method = "shared")
dsgn_tr <- moretrees_design_tree(y = y, X = X, W = W, outcomes = outcomes,
                                    tr = tr)

# starting values ----------------------------------------------
hyper_random_init <- list(omega_max = 100,
                         tau_max = 100,
                         sigma2_max = 100)
vi_random_init <- list(eta_sd = 10,
                      mu_sd = 10,
                      delta_sd = 10)

# Prepare for running algorithm ---------------------------------------------------
set.seed(3457548)
G <- length(dsgn_mat$groups)
n <- length(dsgn_mat$y)
K <- sapply(dsgn_mat$groups, length)
m <- ncol(dsgn_mat$W)
# Initial hyperparameter values
hyperparams <- list(omega = runif(1, 0, hyper_random_init$omega_max),
                      tau = runif(1, 0, hyper_random_init$tau_max),
                      rho = runif(1, 0, 1),
                      sigma2 = runif(1, 0, hyper_random_init$sigma2_max))
# Variational parameter initial values
XtX <- lapply(dsgn_mat$groups, function(cols) crossprod(dsgn_mat$X[ , cols]))
Sigma_inv_mat <- lapply(XtX, FUN = function(XtX, tau, sigma2) XtX / sigma2 + 
                      diag(x = 1 / tau, nrow = ncol(XtX)),
                    tau = hyperparams$tau, 
                    sigma2 = hyperparams$sigma2)  
Sigma_mat <- lapply(Sigma_inv_mat, solve)
Sigma_det_mat <- sapply(Sigma_mat, det)
mu_mat <- lapply(K, rnorm, mean = 0 , sd = vi_random_init$mu_sd)
mu_mat <- lapply(mu_mat, matrix, ncol = 1)
prob <- runif(G, 0, 1)
tau_t <- rep(hyperparams$tau, G) # this should not be changed; tau_t = tau according to algorithm
delta_mat <- matrix(rnorm(m, sd = vi_random_init$delta_sd), ncol = 1)
WtW <- crossprod(dsgn_mat$W)
Omega_inv_mat <- WtW / hyperparams$sigma2 + diag(x = 1 / hyperparams$omega, nrow = m)
Omega_mat <- solve(Omega_inv_mat)
for (v in 1:nrow(Omega_mat)) {
  for (u in 1:ncol(Omega_mat)) {
    if (v == u | v == u + p | u == v + p) next
    Omega_mat[v, u] <- 0
  }
}
Omega_inv_mat <- solve(Omega_mat)
Omega_det_mat <- det(Omega_mat)

# matrix ---------------------------------
X <- dsgn_mat$X
XtX <- lapply(dsgn_mat$groups, function(cols) crossprod(X[ , cols]))
groups <- dsgn_mat$groups
W <- dsgn_mat$W
WtW <- crossprod(W)
y <- dsgn_mat$y
A <- dsgn_mat$A
delta <- delta_mat
Sigma <- Sigma_mat
Sigma_inv <- Sigma_inv_mat
Sigma_det <- Sigma_det_mat
mu <- mu_mat
m <- ncol(W)
K <- sapply(groups, length, simplify = T)
Omega <- Omega_mat
Omega_inv <- Omega_inv_mat
Omega_det <- Omega_det_mat

vi_mat <- update_vi_params_normal(X = X, groups = groups, W = W, y = y,
                                   XtX = XtX, WtW = WtW,
                                     n = n, K = K, G = G, m = m, prob = prob,
                                     mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
                                     Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
                                     Omega = Omega, Omega_inv = Omega_inv, Omega_det = Omega_det,
                                     omega = omega, rho = rho,
                                     tau = tau, sigma2 = sigma2, update_hyper_last = T)

# tr ---------------------------------
X <- dsgn_tr$X
XtX <- lapply(dsgn_tr$outcomes_units, function(units) crossprod(dsgn_tr$X[units , ]))
XtX <- lapply(dsgn_tr$outcomes_nodes, function(outcomes) Reduce(`+`, XtX[outcomes]))
W <- dsgn_tr$W
WtW <- lapply(dsgn_tr$outcomes_units, function(units) crossprod(dsgn_tr$W[units , ]))
WtW <- lapply(dsgn_tr$outcomes_nodes, function(outcomes) Reduce(`+`, WtW[outcomes]))
y <- dsgn_tr$y
A <- dsgn_tr$A
outcomes_units <- dsgn_tr$outcomes_units
outcomes_nodes <- dsgn_tr$outcomes_nodes
ancestors <- dsgn_tr$ancestors
pL <- length(ancestors)
p <- length(unique(unlist(ancestors)))
delta <- sapply(1:p, function(v) matrix(delta_mat[c(1, p + 1) + (v - 1)], ncol = 1),
                simplify = F)
Sigma <- Sigma_mat
Sigma_inv <- Sigma_inv_mat
Sigma_det <- Sigma_det_mat
mu <- mu_mat
K <- ncol(X)
m <- ncol(W)
Omega <- sapply(c(1:(p - 1), 0), function(v, Omega) as.matrix(Omega[1:(m*p) %% p == v,  1:(m*p) %% p == v]), 
                Omega = Omega_mat, simplify = F)
Omega_inv <- sapply(Omega, solve, simplify = F)
Omega_det <- sapply(Omega, det, simplify = T)

vi_tr <- update_vi_params_normal_moretrees(X = X, W = W, 
                                              y = y, WtW = WtW, XtX = XtX,
                                              outcomes_nodes = outcomes_nodes,
                                              outcomes_units = outcomes_units,
                                              ancestors = ancestors,
                                     n = n, K = K, p = p, pL = pL, m = m, prob = prob,
                                     mu = mu, Sigma = Sigma, Sigma_inv = Sigma_inv,
                                     Sigma_det = Sigma_det, tau_t = tau_t, delta = delta,
                                     Omega = Omega, Omega_inv = Omega_inv, 
                                     Omega_det = Omega_det, omega = omega, rho = rho,
                                     tau = tau, sigma2 = sigma2, update_hyper_last = T)
hyper_tr <-  update_hyperparams_logistic_moretrees(X = dsgn$X, 
                                                      W = dsgn$W,
                                                      y = dsgn$y, 
                                                      outcomes_units = dsgn$outcomes_units,
                                                      ancestors = dsgn$ancestors,
                                                      n = n, K = K, p = p, m = m,
                                                      prob = prob, mu = mu,
                                                      Sigma = Sigma, Sigma_det = Sigma_det,
                                                      tau_t = tau_t, delta = delta,
                                                      Omega = Omega, Omega_det = Omega_det,
                                                      eta = hyperparams$eta, g_eta = hyperparams$g_eta,
                                                      omega = hyperparams$omega, tau = hyperparams$tau,
                                                      rho = hyperparams$rho, update_hyper = F)

# compare ---------------------------------

check_sigma <- 0
for (v in 1:p) {
  Sigma_mat_v <- vi_mat$Sigma[[v]]
  Sigma_tr_v <- vi_tr$Sigma[[v]]
  check_v <- all.equal(Sigma_mat_v, Sigma_tr_v, check.attributes = F,
                       check.names = F)
  check_sigma <- check_sigma + (check_v != TRUE)
}

check_mu <- 0
for (v in 1:p) {
  mu_mat_v <- vi_mat$mu[[v]]
  mu_tr_v <- vi_tr$mu[[v]]
  check_v <- all.equal(mu_mat_v, mu_tr_v, check.attributes = F,
                       check.names = F)
  check_mu <- check_mu + (check_v != TRUE)
}

check_mu == 0
check_sigma == 0
all.equal(vi_mat$prob, vi_tr$prob)

delta_vi_mat <- vi_mat$delta
delta_vi_mat <- lapply(c(1:(p-1), 0), 
          FUN = function(v) delta_vi_mat[1:(m * p) %% p == v])

plot(unlist(delta_vi_mat), unlist(vi_tr$delta))

# plot(vi_mat$prob, vi_tr$prob)
# abline(a = 0, b = 1, col = "red")
# 
# plot(unlist(vi_mat$mu), unlist(vi_tr$mu))
# abline(a = 0, b = 1, col = "red")
# 
# plot(unlist(vi_mat$Sigma), unlist(vi_tr$Sigma))
# abline(a = 0, b = 1, col = "red")

vi_mat$delta[c(1, p + 1)]
vi_tr$delta[[1]]
