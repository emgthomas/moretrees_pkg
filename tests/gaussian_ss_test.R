# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

# set.seed(98647)
devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 4
G <- 20 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
m <- 5
tau <- 3
rho <- 0.5
omega <- 2
gamma_true <- Matrix::Matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
beta <- gamma_true * s_true
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
sigma2 <- 2
n <- 300
# Generate some data -----------------------------------------------------------------
X <- sapply(1:G, FUN = function(i) Matrix::Matrix(rnorm(K * n, sd = 0.5), nrow = n))
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (g in 1:G) {
  lp <- lp + X[[g]] %*% beta[g, ]
}
lp <- as.numeric(lp)
y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_normal(y, X, W, update_hyper = T, update_hyper_freq = 50,
                              tol = 1E-8, max_iter = 5000,
                              hyperparams_init = list(omega = omega,
                                                      rho = rho,
                                                      tau = tau,
                                                      sigma2 = sigma2))

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
plot(mod1$vi_params$mu * mod1$vi_params$prob, beta)
abline(a = 0, b = 1, col = "red")

# Compare non-sparse effect estimates to truth ----------------------------------------------
plot(mod1$vi_params$delta, theta)
abline(a = 0, b = 1, col = "red")

# Compare moretrees estimates to maximum likelihood -----------------------------------------
X1 <- X[[1]]
for (g in 2:G) {
  X1 <- cbind(X1, X[[g]])
}
# mod2 <- lm(y ~ 0 + as.matrix(W) + as.matrix(X1))
mod2 <- lm(y ~ 0 + as.matrix(X1))
plot(mod2$coefficients, c(theta, as.numeric(t(beta))))
abline(a = 0, b = 1, col = "red")
plot(mod2$coefficients, c(as.numeric(mod1$vi_params$delta), as.numeric(t(mod1$vi_params$mu))))
abline(a = 0, b = 1, col = "red")

# Compare hyperparameter estimates to truth -------------------------------------------------

cbind(mod1$hyperparams[2:5], c(omega, sigma2, tau, rho))

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
plot(plot_start:length(mod1$sigma2_track),
     mod1$sigma2_track[plot_start:length(mod1$sigma2_track)],
     type = "l")

