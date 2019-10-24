# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 2
G <- 200 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
beta <- matrix(sample(c(-1, 1), replace = T, size = G * K,
                      prob = c(0.25, 0.75)), nrow = G)
s_true <- sample(c(0, 1), replace = T, size = G)
beta <- beta * s_true
sigma2 <- 0.5
n <- 500
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K * G * n), dim = c(G, n, K))
lp <- numeric(n) + 0
for (g in 1:G) {
  lp <- lp + matrix(X[g, , ], nrow = n) %*% matrix(beta[g, ], nrow = K)
}
lp <- as.numeric(lp)
y <- lp + rnorm(n, sigma2)
# Run algorithm ----------------------------------------------------------------------
mod1 <- moretrees_normal(y, X, update_hyper = T)
# Plot results -----------------------------------------------------------------------
plot_start <- 50
plot(plot_start:length(mod1$ELBO),
     mod1$ELBO[plot_start:length(mod1$ELBO)],
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
tapply(mod1$params$prob, s_true, summary)
plot(mod1$params$prob, s_true)
plot(mod1$params$mu, beta)
