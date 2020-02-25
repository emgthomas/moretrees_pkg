# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian and Bernoulli outcomes  ----------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

rm(list = ls())
devtools::load_all() # Sources all files in R/

# Chose one --------------------------------------------------------------------------
# family <- "gaussian"
family <- "bernoulli"

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
n <- 1E6
K_g <- 3 # number of variables
K <- rep(K_g, G)
m <- 2
tau <- 3
rho1 <- 0.6 # rho for internal nodes
rho2 <- 0.05 # rho for leaf nodes
rho <- sum(1 + 0.8 * (p - pL - 1) + 0.05 * pL) / p # overall rho
omega <- 2
sigma2 <- 2
hyper_fixed <- list(tau = tau, rho = rho, omega = omega)
if (family == "gaussian") hyper_fixed$sigma2 <- sigma2
nrestarts <- 1
doParallel::registerDoParallel(cores = nrestarts)

# Generate randomly grouped beta (groups follow tree)
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
s_true <- c(1, rbinom(n = p - pL - 1, size = 1, prob = rho1), 
            rbinom(n = pL, size = 1, prob = rho2))
xi <- mapply(function(gamma, s) matrix(gamma * s, nrow = 1),
             gamma = gamma_true, s = s_true,
             SIMPLIFY = T) %>% t
if (K_g == 1) xi <- t(xi)
beta <- A[leaves, ] %*% xi
zeta <- matrix(rnorm(m * G, mean = 0, sd = sqrt(omega)), nrow = G, ncol = m)
theta <- A[leaves, ] %*% zeta
groups_true <- as.integer(as.factor(as.numeric(beta[ , 1])))
table(groups_true)

# Generate some data -----------------------------------------------------------------
X <- Matrix::Matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
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
if (family == "gaussian") {
  y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
} else {
  p_success <- expit(lp)
  y <- runif(n)
  y <- as.integer(y <= p_success)
}

require(gdata)
keep(X, W, y, outcomes, tr, family, hyper_fixed, nrestarts, sure = T,
     s_true, beta, theta, groups_true)

# Run algorithm ----------------------------------------------------------------------
require(profvis)
profvis(
  mod <- moretrees(X = X, W = W, y = y, outcomes = outcomes,
                   W_method = "shared",
                   tr = tr, family = family,
                   update_hyper = T, update_hyper_freq = 1,
                   hyper_fixed = hyper_fixed,
                   tol = 1E-8, max_iter = 4,
                   print_freq = 1,
                   nrestarts = nrestarts,
                   get_ml = F,
                   log_dir = "./tests/")
)
beta_est <- mod$beta_est
beta_moretrees <- mod$beta_moretrees
beta_ml <- mod$beta_ml
theta_est <- mod$theta_est
mod_restarts <- mod$mod_restarts
mod1 <- mod$mod

# Compare ELBOs for random restarts --------------------------------------------------
c(mod1$ELBO_track[length(mod1$ELBO_track)],
  sapply(mod_restarts, function(mod) mod$ELBO_track[length(mod$ELBO_track)]))

# Plot results -----------------------------------------------------------------------

# Check if the ELBO decreases
ELBO_track <- mod1$ELBO_track
if(min(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)]) < 0) {
  print("ELBO decreases at these time points:")
  which(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)] < 0)
} else {
  print("ELBO always increases!")
}

# ELBO at every time step
plot_start <- 1
plot_end <- length(ELBO_track)
# plot_end <- 240
plot(plot_start:plot_end,
     ELBO_track[plot_start:plot_end],
     type = "l")

# Compare sparse effect estimates to truth --------------------------------------------------
# compare estimated probabilities of variable inclusion to true variable inclusion indicators
tapply(mod1$vi_params$prob, s_true, summary)
plot(mod1$vi_params$prob, s_true)
# compare estimated coefficients to true coefficients
k <- 1
clmn <- paste0("est", k)
plot(beta_est[ , clmn <- paste0("est", k)], beta[ , k])
abline(a = 0, b = 1, col = "red")
j <- 1
clmn <- paste0("est", j)
plot(theta_est[, clmn], theta[ , j])
abline(a = 0, b = 1, col = "red")

# Compare estimated groups to truth ---------------------------------------------------------
groups_est <- beta_est$group
table(groups_est, groups_true)

# Compare moretrees estimates to maximum likelihood -----------------------------------------
k <- 1
clmn <- paste0("est", k)
plot(beta_ml[, clmn], beta_moretrees[ , clmn])
abline(a = 0, b = 1, col = "red")
cbind(beta_ml[, clmn], beta_moretrees[ , clmn])

# Compare hyperparameter estimates to truth -------------------------------------------------
if (family == "gaussian") {
  cbind(mod1$hyperparams[2:5], c(hyper_fixed$omega, hyper_fixed$sigma2, hyper_fixed$tau, hyper_fixed$rho))
} else {
  cbind(mod1$hyperparams[2:4], c(hyper_fixed$omega, hyper_fixed$tau, hyper_fixed$rho))
}

