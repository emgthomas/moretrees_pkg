# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

# set.seed(98647)
devtools::load_all() # Sources all files in R/

# Pick one ---------------------------------------------------------------------------
family <- "gaussian"
# family <- "bernoulli"

# Input parameters -------------------------------------------------------------------
G <- 4 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
K <- sample(1:4, size = G, replace = T)
m <- 4
tau <- 3
rho <- 0.5
omega <- 2
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau))
# s_true <- rbinom(n = G, size = 1, prob = rho)
s_true <- sample(c(0, 0, 1, 1), size = )
beta <- sapply(1:G, function(g) gamma_true[[g]] * s_true[[g]])
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
sigma2 <- 2
n <- 500
# Generate some data -----------------------------------------------------------------
X <- sapply(K, FUN = function(k) Matrix::Matrix(rnorm(k * n, sd = 0.5), nrow = n))
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (g in 1:G) {
  lp <- lp + X[[g]] %*% beta[[g]]
}
lp <- as.numeric(lp)
if(family == "bernoulli") {
  expit <- 1 / (1 + exp(-lp))
  y <- sapply(expit, rbinom, n = 1, size = 1)
} else {
  y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
}

# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab(y, X, W, family = family,
                       update_hyper = T, update_hyper_freq = 50,
                              tol = 1E-8, max_iter = 5000)
                              # hyperparams_init = list(omega = omega,
                              #                         rho = rho,
                              #                         tau = tau,
                              #                         sigma2 = sigma2))
beta_est <- mod1$sparse_est
theta_est <- mod1$nonsparse_est
mod1 <- mod1$mod

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
ss_est <- sapply(beta_est, function(b) b$est) %>% unlist
plot(ss_est, unlist(beta))
abline(a = 0, b = 1, col = "red")

# Compare non-sparse effect estimates to truth ----------------------------------------------
plot(theta_est$est, theta)
abline(a = 0, b = 1, col = "red")

# Compare estimates to maximum likelihood ------------------------------------------------
X1 <- X[[1]]
for (g in 2:G) {
  X1 <- cbind(X1, X[[g]])
}
if (family == "bernoulli") {
  family2 <- "binomial"
} else {
  family2 <- family
}
mod2 <- glm(y ~ 0 + as.matrix(W) + as.matrix(X1), family = family2)
plot(mod2$coefficients, c(theta, unlist(beta)))
abline(a = 0, b = 1, col = "red")
plot(mod2$coefficients, c(theta_est$est,
                          ss_est))
abline(a = 0, b = 1, col = "red")

# Compare hyperparameter estimates to truth -------------------------------------------------
if (family == "gaussian") {
  cbind(mod1$hyperparams[2:5], c(omega, sigma2, tau, rho))
} else {
  cbind(mod1$hyperparams[2:4], c(omega, tau, rho))
}
