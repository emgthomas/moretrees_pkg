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
tau <- 0.1
rho <- 0.5
alpha <- 1
tau_alpha <- 1
gamma_true <- matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
# s_true <- rep(1, G)
# s_true[1] <- 0
# s_true <- c(1,1,0)
# rho <- mean(s_true)
rho <- 0.5
beta <- gamma_true * s_true
# beta[1,] <- rep(0, K)
n <- 100
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K * G * n, sd = 1), dim = c(G, n, K))
lp <- numeric(n) + 0
for (g in 1:G) {
  lp <- lp + matrix(X[g, , ], nrow = n) %*% matrix(beta[g, ], nrow = K)
}
lambda <- exp(alpha + as.numeric(lp))
y <- sapply(lambda, rpois, n = 1)
summary(y)

# # Other data inputs ------------------------------------------------------------------
# logfac <- function(x) {
#   if(x == 0) return(0)
#   sum(log(1:x))
# }
# sum_log_y_fac <- sum(sapply(y, logfac))
# Xy <- matrix(nrow = G, ncol = K)
# for (g in 1:G) {
#   Xy[g, ] <- t(X[G, , ]) %*% y
# }

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
mod1 <- spike_and_slab_poisson(y, X, tol = 1E-16, max_iter = 500,
                               update_hyper_freq = 1,
                              hyperparams_init = list(rho = rho,
                                                      tau = tau),
                              s_true = s_true)

# Try maximising ELBO directly -------------------------------------------------------
{
# Computing some constant values
sum_log_y_fac <- sum(sapply(y, logfac))
# Initial values
u <- rep(5, G)
u[s_true == 0] <- -5
# mu <- matrix(rnorm(G * K, sd = 1), nrow = G)
mu <- beta
Sigma <- array(dim = c(G, K, K))
for (g in 1:G){
  Sigma[g, , ] <- matrix(0.1, nrow = K, ncol =K) + diag(0.9, nrow = K)
}
log_Sigma <- log(Sigma)
# mu_alpha <- log(mean(y))
mu_alpha <- alpha
log_tau_t_alpha <- log(0.1)
# Put VI parameters in list
VIparams <- c(u = u, mu = mu, log_Sigma = log_Sigma,
                  mu_alpha = mu_alpha, log_tau_t_alpha = log_tau_t_alpha)

# maximise
ELBO_wrap <- function(VIparams, X, y, n, K, G, sum_log_y_fac, rho, tau) {
  nm <- names(VIparams)
  mu_alpha <- VIparams["mu_alpha"]
  tau_t_alpha <- exp(VIparams["log_tau_t_alpha"])
  u_nm <- sapply(1:G, function(g) paste0("u", g))
  u <- VIparams[u_nm]
  mu_nm <- sapply(1:(G * K), function(i) paste0("mu", i))
  mu <- matrix(VIparams[mu_nm], nrow = G)
  Sigma_nm <- sapply(1:(G * K * K), function(i) paste0("log_Sigma", i))
  Sigma <- array(exp(VIparams[Sigma_nm]), dim = c(G, K, K))
  elbo_poisson_maxim(X = X, y = y, n = n, K = K, G = G,
               sum_log_y_fac = sum_log_y_fac, 
               u = u, mu = mu, Sigma = Sigma,
               mu_alpha = mu_alpha, tau_t_alpha = tau_t_alpha,
               tau = tau, rho = rho)
}

ELBO_wrap(VIparams, X, y, n, K, G, sum_log_y_fac, rho, tau)

mod2 <- optim(par = VIparams, fn = ELBO_wrap, 
              X = X, y = y, n = n, K = K, G = G,
              sum_log_y_fac = sum_log_y_fac, 
              rho = rho, tau = tau,
              control = list(maxit = 10000))
mod2$convergence

res <- mod2$par
nm <- names(res)
mu_alpha_res <- res["mu_alpha"]
tau_t_alpha_res <- exp(res["log_tau_t_alpha"])
u_nm <- sapply(1:G, function(g) paste0("u", g))
u_res <- res[u_nm]
prob_res <- 1/(1+exp(-u_res))
mu_nm <- sapply(1:(G * K), function(i) paste0("mu", i))
mu_res <- matrix(res[mu_nm], nrow = G)
Sigma_nm <- sapply(1:(G * K * K), function(i) paste0("log_Sigma", i))
Sigma_res <- array(exp(res[Sigma_nm]), dim = c(G, K, K))
}

# Plot results -----------------------------------------------------------------------

# ELBO at every time step
ELBO_track <- mod1$ELBO_track
plot_start <- 1
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

lm(as.numeric(mu_est) ~ 0 + as.numeric(beta))
# error is a factor of around 2.8

X2 <- matrix(nrow = n, ncol = K * G)
for (g in 1:G) {
  for (k in 1:K) {
    X2[ , K * (g - 1) + k] <- X[g, , k]
  }
}
mod2 <- glm(y ~ X2, family = poisson)
cbind(mod2$coefficients, c(mu_alpha_est, as.numeric(t(mu_est))))
plot(mod2$coefficients, c(mu_alpha_est, as.numeric(t(mu_est))))

# plot(c(alpha, as.numeric(t(beta))), c(mu_alpha_est, as.numeric(t(mu_est))))
# plot(mod2$coefficients, c(alpha, as.numeric(t(beta))))


