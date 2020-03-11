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
group <- "7.4"
tr <- ccs_tree(group)$tr
leaves <- names(igraph::V(tr)[igraph::V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- Matrix::expm(Matrix::t(A))
A[A > 0 ] <- 1
G <- length(igraph::V(tr))
p <- G
pL <- sum(igraph::V(tr)$leaf)
n <- 1E3
K_g <- 1 # number of variables
K <- rep(K_g, G)
m <- 2
# mdim <- 3
tau <- 3
rho1 <- 0.4 # rho for internal nodes
rho2 <- 0.1 # rho for leaf nodes
rho <- sum(1 + 0.8 * (p - pL - 1) + 0.05 * pL) / p # overall rho
omega <- 2
sigma2 <- 2
hyper_fixed <- list(a_tau = c(0.001, 0.001, 0.001, 0.001), 
                    b_tau = c(0.001, 0.001, 0.001, 0.001),
                    a_omega = c(0.001, 0.001, 0.001, 0.001), 
                    b_omega = c(0.001, 0.001, 0.001, 0.001))
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
X <- matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- c(leaves, sample(leaves, size = n - pL, replace = T))

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
  lp[which_v] <- lp[which_v] + X[which_v, , drop = F] %*% Matrix::t(beta[v, , drop = F])
  if (m > 0) {
    lp[which_v] <- lp[which_v] + W[which_v, , drop = F] %*% Matrix::t(theta[v, , drop = F])
    # counter <- 0
    # for (j in 1:m) {
    #   for (l in 1:mdim) {
    #     counter <- counter + 1
    #     lp[which_v] <- lp[which_v] + W[which_v, j , drop = F] ^ (l - 1) %*% 
    #       Matrix::t(theta[v, counter, drop = F])
    #   }
    # }
  }
}

# v <- leaves[3]
# plot(W[outcomes == v, 2], lp[outcomes == v])

# Simulate outcomes
if (family == "gaussian") {
  y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
} else {
  p_success <- expit(lp)
  y <- runif(n)
  y <- as.integer(y <= p_success)
}

# require(splines)
# df <- 3
# Wspl <- matrix(nrow = n, ncol = 0)
# for (j in 1:m) {
#   Wspl <- cbind(Wspl, ns(W[ , j], df = df))
# }
# W <- Wspl


# Run algorithm ----------------------------------------------------------------------
require(gdata)
keep(X, W, y, outcomes, tr, family, hyper_fixed, nrestarts, hyper_fixed,
     s_true, groups_true, beta, theta, sure = T)
# require(profvis)
# profvis(
# Run model without W
  mod_start <- moretrees(X = X, W = W, y = y, outcomes = outcomes,
                   random_init = F,
                   method = "tree",
                   W_method = "shared",
                   tr = tr, family = family,
                   update_hyper = T, update_hyper_freq = 50,
                   hyper_fixed = hyper_fixed,
                   tol = 1E-8, max_iter = 1E5,
                   print_freq = 10,
                   nrestarts = nrestarts,
                   get_ml = T,
                   log_dir = "./tests/")
  mod_end <- mod_start
  # mod_end <- moretrees(X = X, W = W, y = y, outcomes = outcomes,
  #                  initial_values = mod_start$mod,
  #                  method = "tree",
  #                  W_method = "shared",
  #                  tr = tr, family = family,
  #                  update_hyper = T, update_hyper_freq = 20,
  #                  hyper_fixed = hyper_fixed,
  #                  tol = 1E-6, max_iter = 1E4,
  #                  print_freq = 10,
  #                  nrestarts = nrestarts,
  #                  get_ml = T,
  #                  log_dir = "./tests/")
# )
beta_est <- mod_end$beta_est
beta_moretrees <- mod_end$beta_moretrees
beta_ml <- mod_end$beta_ml
theta_est <- mod_end$theta_est
mod_restarts <- mod_end$mod_restarts
mod1 <- mod_end$mod

# Compare ELBOs for random restarts --------------------------------------------------
c(mod1$ELBO_track[length(mod1$ELBO_track)],
  sapply(mod_restarts, function(mod) mod$ELBO_track[length(mod$ELBO_track)]))

# Plot results -----------------------------------------------------------------------

# Check if the ELBO decreases
ELBO_track <- mod1$ELBO_track
# ELBO_track <- c(mod_start$mod$ELBO_track, mod_end$mod$ELBO_track[2:length(mod_end$mod$ELBO_track)])
if(min(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)]) < 0) {
  print("ELBO decreases at these time points:")
  which(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)] < 0)
} else {
  print("ELBO always increases!")
}

# ELBO at every time step
plot_start <- 1
plot_end <- length(ELBO_track)
# plot_end <- 4020
plot(ELBO_track[plot_start:plot_end],
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
  cbind(mod1$hyperparams[2:3], c(hyper_fixed$omega, hyper_fixed$sigma2, hyper_fixed$tau))
} else {
  cbind(mod1$hyperparams[2:3], c(hyper_fixed$omega, hyper_fixed$tau))
}

