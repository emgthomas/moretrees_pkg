# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian and Bernoulli outcomes  ----------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #
devtools::load_all() # Sources all files in R/

require(igraph)
require(Matrix)

# Chose one --------------------------------------------------------------------------
# family <- "gaussian"
family <- "bernoulli"

# Input parameters -------------------------------------------------------------------
group <- "7.3"
tr <- ccs_tree(group)$tr
leaves <- names(V(tr)[V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- expm(Matrix::t(A))
A[A > 0 ] <- 1 
G <- length(V(tr))
n <- 500
K_g <- 2 # number of variables
K <- rep(K_g, G)
m <- 0
tau <- 3
rho <- 0.5
omega <- 2
sigma2 <- 2

# Generate randomly grouped beta (groups follow tree)
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
# s_true <- rbinom(n = G, size = 1, prob = rho)
# s_true[1] <- 1
s_true <- rep(0, G)
s_true[c(1, 2, 5)] <- 1
xi <- sapply(1:G, function(g) Matrix::Matrix(gamma_true[[g]] * s_true[[g]],
                                             ncol = 1))
xi1 <- sapply(xi, function(xi) xi[1 , 1])
xi2 <- sapply(xi, function(xi) xi[2 , 1])
xi <- cbind(xi1, xi2)
beta <- A[leaves, ] %*% xi
beta1 <- beta[ , 1]
beta2 <- beta[ , 2]
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
groups_true <- as.integer(as.factor(as.numeric(beta1)))
table(groups_true)

# Generate some data -----------------------------------------------------------------
X <- matrix(rnorm(n * K_g), nrow = n)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)

# Get linear predictor
lp <- W %*% theta
for (v in leaves) {
  which_v <- outcomes == v
  lp[which_v] <- lp[which_v] + X[which_v, ] %*% beta[v , ]
}
lp <- as.numeric(lp)

# Simulate outcomes
if (family == "gaussian") {
  y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
} else {
  p_success <- expit(lp)
  y <- sapply(p_success, rbinom, n = 1, size = 1)
}


# Run algorithm ----------------------------------------------------------------------
# Create MORETreeS design matrix
mod <- moretrees(X = X, y = y, outcomes = outcomes, 
                  tr = tr, family = family,
                  update_hyper = T, update_hyper_freq = 10,
                  tol = 1E-8, max_iter = 1E5,
                  get_ml = T)
beta_est <- mod$beta_est
beta_moretrees <- mod$beta_moretrees
beta_ml <- mod$beta_ml
mod1 <- mod$mod

# Plot results -----------------------------------------------------------------------

# Check if the ELBO decreases
ELBO_track <- mod1$ELBO_track2
if(min(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)]) < 0) {
  print("ELBO decreases at these time points:")
  which(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)] < 0)
} else {
  print("ELBO always increases!")
}

# ELBO at every time step
plot_start <- 100
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
plot(beta_est[ , 1], beta1)
abline(a = 0, b = 1, col = "red")
plot(beta_est[ , 2], beta2)
abline(a = 0, b = 1, col = "red")

# Compare estimated groups to truth ---------------------------------------------------------
groups_est <- beta_est$group
table(groups_est, groups_true)

# Compare moretrees estimates to maximum likelihood -----------------------------------------
plot(beta_ml$est1, beta_moretrees$est1)
abline(a = 0, b = 1, col = "red")
cbind(beta_ml$est1, beta_moretrees$est1)
plot(beta_ml$est2, beta_moretrees$est2)
abline(a = 0, b = 1, col = "red")
cbind(beta_ml$est2, beta_moretrees$est2)

# Compare hyperparameter estimates to truth -------------------------------------------------
if (family == "gaussian") {
  cbind(mod1$hyperparams[2:5], c(omega, sigma2, tau, rho))
} else {
  cbind(mod1$hyperparams[2:4], c(omega, tau, rho))
}


# ELBO when hyperparams updated
plot_start <- 1
plot(plot_start:length(mod1$ELBO_track),
     mod1$ELBO_track[plot_start:length(mod1$ELBO_track)],
     type = "l")
plot(plot_start:length(mod1$rho_track),
     mod1$rho_track[plot_start:length(mod1$rho_track)],
     type = "l")
plot(plot_start:length(mod1$tau_track),
     mod1$tau_track[plot_start:length(mod1$tau_track)],
     type = "l")
if (family == "gaussian") {
  plot(plot_start:length(mod1$sigma2_track),
       mod1$sigma2_track[plot_start:length(mod1$sigma2_track)],
       type = "l")
}

