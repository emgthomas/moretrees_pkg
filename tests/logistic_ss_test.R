# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

# set.seed(98647)
devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
G <- 10 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
K <- sample(1:4, size = G, replace = T)
m <- 10
tau <- 3
rho <- 0.5
omega <- 2
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau))
s_true <- rbinom(n = G, size = 1, prob = rho)
print(s_true)
beta <- sapply(1:G, function(g) gamma_true[[g]] * s_true[[g]])
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
n <- 1000
# Generate some data -----------------------------------------------------------------
X <- sapply(K, FUN = function(k) Matrix::Matrix(rnorm(k * n, sd = 0.5), nrow = n))
W <- matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (g in 1:G) {
  lp <- lp + X[[g]] %*% beta[[g]]
}
lp <- as.numeric(lp)
expit <- 1 / (1 + exp(-lp))
y <- sapply(expit, rbinom, n = 1, size = 1)
y[y == 0] <- -1
# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_logistic(y, X, W, update_hyper = F, update_hyper_freq = 10,
                              tol = 1E-8, max_iter = 5000)
                              # hyperparams_init = list(omega = omega,
                              #                         rho = rho,
                              #                         tau = tau))

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
moretrees_est <- unlist(sapply(1:G, 
                               function(g) as.numeric(mod1$vi_params$mu[[g]] * mod1$vi_params$prob[g])))
plot(moretrees_est, unlist(beta))
abline(a = 0, b = 1, col = "red")

# Compare non-sparse effect estimates to truth ----------------------------------------------
plot(mod1$vi_params$delta, theta)
abline(a = 0, b = 1, col = "red")

# Compare moretrees estimates to ML estimates -----------------------------------------------
X1 <- X[[1]] * (mod1$vi_params$prob[1] > 0.5)
for (g in 2:G) {
  X1 <- cbind(X1, X[[g]] * (mod1$vi_params$prob[g] > 0.5) )
}
y2 <- y
y2[y == -1] <- 0
mod2 <- glm(y2 ~ 0 + as.matrix(W) + as.matrix(X1), family = binomial)
# mod2 <- glm(y2 ~ 0 + as.matrix(X1), family = binomial)
mod2$coefficients
plot(mod2$coefficients, c(theta, unlist(beta)))
abline(a = 0, b = 1, col = "red")
plot(mod2$coefficients, c(as.numeric(mod1$vi_params$delta),
                          moretrees_est))
abline(a = 0, b = 1, col = "red")

# Compare hyperparameter estimates to truth -------------------------------------------------

cbind(mod1$hyperparams[2:4], c(omega, tau, rho))

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
plot(plot_start:length(mod1$omega_track),
     mod1$omega_track[plot_start:length(mod1$omega_track)],
     type = "l")

