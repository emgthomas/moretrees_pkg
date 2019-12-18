# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 2
G <- 5 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
# beta <- matrix(sample(c(-1, 1), replace = T, size = G * K,
#                       prob = c(0.25, 0.75)), nrow = G)
tau <- 1
rho <- 0.5
alpha <- 1
tau_alpha <- 1
gamma_true <- matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
# s_true <- c(1,0)
beta <- gamma_true * s_true
n <- 200
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

# Initial values ---------------------------------------------------------------------
# prob <- 0.9 * s_true + 0.1 * (1 - s_true)
prob <- rep(0.5, G)
mu <- beta + matrix(rnorm(K * G, sd = 1), nrow = G, ncol = K)
Sigma <- matrix(0.1, nrow = G, ncol = K)
tau_t <- tau + rnorm(1, sd = 0.5)
mu_alpha <- alpha + rnorm(1, sd = 0.5)
tau_t_alpha <- 0.1
par <- c(log(prob / (1 - prob)), mu, log(Sigma), log(tau_t), mu_alpha, log(tau_t_alpha))

elbo_poisson(par = par, 
             X = X, y = y, n = n, K = K, G = G, sum_log_y_fac = sum_log_y_fac, # data
             tau = tau, rho = rho, tau_alpha = tau_alpha)

# Run algorithm ----------------------------------------------------------------------
mod1 <- optim(par, elbo_poisson, 
      X = X, y = y, n = n, K = K, G = G, sum_log_y_fac = sum_log_y_fac, # data
      tau = tau, rho = rho, tau_alpha = tau_alpha, # hyperparameters
      control = list(maxit = 20000, fnscale = -1) # fnscale = -1 ensures function 
      # is maximized, not minimized
      )
mod1$convergence

# Compare parameter estimates to truth -----------------------------------------------

par <- mod1$par
prob_est <- 1 / (1 + exp(-par[1:G]))
mu_est <- matrix(par[(G + 1):(G + G * K)], nrow = G)
Sigma_est <- exp(matrix(par[(G + G * K + 1):(G + G * K + G * K)], nrow = G))
tau_t_est <- exp(par[G + G * K + G * K + 1])
mu_alpha_est <- par[G + G * K + G * K + 2]
tau_t_alpha_est <- exp(par[G + G * K + G * K + 3])

c(mu_alpha_est, alpha)
tapply(prob_est, s_true, summary)
plot(prob_est, s_true)
abline(a = 0, b = 1, col = "red")
plot(mu_est * prob_est, beta)
abline(a = 0, b = 1, col = "red")
cbind(as.numeric(mu_est * prob_est), as.numeric(beta))
cbind(Sigma_est, prob_est)

# Compare parameter estimates to maximum likelihood ----------------------------------
X2 <- X[1, , ]
for (g in 2:G) {
  X2 <- cbind(X2, X[g, , ] * (prob_est[g] > 0.5))
}
mod2 <- glm(y ~ X2, family = poisson)
plot(c(alpha, as.numeric(t(beta))), mod2$coeff)
abline(a = 0, b = 1, col = "red")
cbind(c(mu_alpha_est, as.numeric(t(mu_est * prob_est))), mod2$coeff)
plot(c(mu_alpha_est, as.numeric(t(mu_est * prob_est))), mod2$coeff)
abline(a = 0, b = 1, col = "red")
