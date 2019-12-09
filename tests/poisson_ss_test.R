# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 3
G <- 10 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
# beta <- matrix(sample(c(-1, 1), replace = T, size = G * K,
#                       prob = c(0.25, 0.75)), nrow = G)
tau <- 2
# rho <- 0.5
alpha <- 1
tau_alpha <- 1
gamma_true <- matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
# s_true <- c(1,1,0)
rho <- mean(s_true)
beta <- gamma_true * s_true
n <- 500
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K * G * n, sd = 1), dim = c(G, n, K))
lp <- numeric(n) + 0
for (g in 1:G) {
  lp <- lp + matrix(X[g, , ], nrow = n) %*% matrix(beta[g, ], nrow = K)
}
lambda <- exp(alpha + as.numeric(lp))
y <- sapply(lambda, rpois, n = 1)

# Other data inputs ------------------------------------------------------------------
logfac <- function(x) {
  if(x == 0) return(0)
  sum(log(1:x))
}
sum_log_y_fac <- sum(sapply(y, logfac))
Xy <- matrix(nrow = G, ncol = K)
for (g in 1:G) {
  Xy[g, ] <- t(X[G, , ]) %*% y
}

# Initial values ---------------------------------------------------------------------
# prob <- 0.9 * s_true + 0.1 * (1 - s_true)
# u <- log(prob / (1 - prob))
# mu <- beta
# Sigma <- array(dim = c(G, K, K))
# for (g in 1:G){
#   Sigma[g, , ] <- diag(1, nrow = K)
# }
# mu_alpha <- alpha
# tau_t_alpha <- 0.1
# expA <- matrix(nrow = n, ncol = G)
# b <- numeric(n)
# for (g in 1:G) {
#   for (i in 1:n) {
#     b[i] <- t(X[g, i, ]) %*% Sigma[g, , ] %*% X[g, i, ]
#   }
#   expA[, g] <- exp(X[g, , ] %*% mu[g, ]  + b / 2)
# }

# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_poisson(y, X, tol = 1E-8, max_iter = 100,
                               update_hyper_freq = 1,
                              hyperparams_init = list(rho = rho,
                                                      tau = tau))

# Plot results -----------------------------------------------------------------------

# ELBO at every time step
ELBO_track <- mod1$ELBO_track
plot_start <- 5
plot_end <- length(ELBO_track)
# plot_end <- 80
plot(plot_start:plot_end,
     ELBO_track[plot_start:plot_end],
     type = "l")

# Compare parameter estimates to truth -----------------------------------------------

# Compare effect estimates to truth
prob_est <- 1 / (1 + exp(-mod1$vi_params$u))
mu_est <- mod1$vi_params$mu
mu_alpha_est <- mod1$vi_params$mu_alpha
Sigma_est <- mod1$vi_params$Sigma

c(mu_alpha_est, alpha)
tapply(prob_est, s_true, summary)
plot(prob_est, s_true)
plot(mu_est * prob_est, beta)
cbind(as.numeric(mu_est * prob_est), as.numeric(beta))
cbind(Sigma_est, prob_est)
