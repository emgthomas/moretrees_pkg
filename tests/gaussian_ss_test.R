# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 2
G <- 100 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
# beta <- matrix(sample(c(-1, 1), replace = T, size = G * K,
#                       prob = c(0.25, 0.75)), nrow = G)
tau <- 3
rho <- 0.3
gamma_true <- matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
beta <- gamma_true * s_true
sigma2 <- 2
n <- 1000
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K * G * n), dim = c(G, n, K))
lp <- numeric(n) + 0
for (g in 1:G) {
  lp <- lp + matrix(X[g, , ], nrow = n) %*% matrix(beta[g, ], nrow = K)
}
lp <- as.numeric(lp)
y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_normal(y, X, update_hyper = F, update_hyper_freq = 50,
                              tol = 1E-8,
                              hyperparams_init = list(rho = rho,
                                                      tau = tau,
                                                      sigma2 = sigma2))
# Compare hyperparams to truth --------------------------------------------------------

cbind(mod1$hyperparams[2:4], c(sigma2, tau, rho))

# Plot results -----------------------------------------------------------------------

# ELBO at every time step
ELBO_track <- mod1$ELBO_track2
plot_start <- 1
plot_end <- length(ELBO_track)
plot_end <- 30
plot(plot_start:plot_end,
     ELBO_track[plot_start:plot_end],
     type = "l")
min(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)])
which(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)] < 0)

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

# Compare effect estimates to truth
tapply(mod1$vi_params$prob, s_true, summary)
plot(mod1$vi_params$prob, s_true)
plot(mod1$vi_params$mu * mod1$vi_params$prob, beta)
