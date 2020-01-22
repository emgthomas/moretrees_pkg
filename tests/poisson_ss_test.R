# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 2
G <- 20 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
# beta <- matrix(sample(c(-1, 1), replace = T, size = G * K,
#                       prob = c(0.25, 0.75)), nrow = G)
tau <- 0.06
rho <- 0.5
alpha <- 1
tau_alpha <- 1
gamma_true <- matrix(rnorm(G * K, mean = 0, sd = sqrt(tau)), nrow = G)
s_true <- rbinom(n = G, size = 1, prob = rho)
# s_true <- rep(1, G)
# s_true[1] <- 0
# s_true <- c(1,1,0)
# rho <- mean(s_true)
beta <- gamma_true * s_true
# beta[1,] <- rep(0, K)
n <- 1000
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K * G * n, sd = 1), dim = c(G, n, K))
lp <- numeric(n) + 0
for (g in 1:G) {
  lp <- lp + matrix(X[g, , ], nrow = n) %*% matrix(beta[g, ], nrow = K)
}
lambda <- exp(alpha + as.numeric(lp))
y <- sapply(lambda, rpois, n = 1)
summary(y)

# Initial values ---------------------------------------------------------------------
# u <- rep(5, G)
# u[s_true == 0] <- -5
# u <- rnorm(G, sd = 5)
u <- rep(5, G)
# mu <- matrix(rnorm(G * K, sd = 1), nrow = G)
# mu <- beta + matrix(rnorm(G * K, sd = 1), nrow = G)
# mu <- beta
X2 <- matrix(nrow = n, ncol = K * G)
for (g in 1:G) {
  for (k in 1:K) {
    X2[ , K * (g - 1) + k] <- X[g, , k]
  }
}
mod3 <- glm(y ~ X2, family = poisson)
mu <- matrix(mod3$coefficients[2:(K * G + 1)], nrow = G, byrow = T)
mu_alpha <- mod3$coefficients[1]
log_tau_t_alpha <- log(0.1)
tau_t_alpha <- exp(log_tau_t_alpha)
Sigma <- array(dim = c(G, K, K))
a <- 0.01
b <- 0.1
for (g in 1:G){
  Sigma[g, , ] <- matrix(a, nrow = K, ncol =K) + diag(b - a, nrow = K)
}
log_Sigma <- log(Sigma)
Sigma2 <- matrix(b, nrow = G, ncol = K)
log_Sigma2 <- log(Sigma2)
# mu_alpha <- log(mean(y))
# mu_alpha <- alpha + rnorm(1, sd = 0.1)

# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_poisson(y, X, tol = 1E-16, max_iter = 1000,
                               update_hyper_freq = 1,
                              hyperparams_init = list(rho = rho,
                                                      tau = tau,
                                                      tau_alpha = tau_alpha),
                              vi_params_init = list(u = u, mu = mu, Sigma = Sigma,
                                                    mu_alpha = mu_alpha, 
                                                    tau_t_alpha = tau_t_alpha))

# Try maximising ELBO directly -------------------------------------------------------
{
# Computing some constant values
sum_log_y_fac <- sum(sapply(y, logfac))
# Put VI parameters in list
VIparams <- c(u = u, mu = mu, log_Sigma = log_Sigma,
                  mu_alpha = mu_alpha, log_tau_t_alpha = log_tau_t_alpha)
VIparams2 <- c(u = u, mu = mu, log_Sigma = log_Sigma2,
              mu_alpha = mu_alpha, log_tau_t_alpha = log_tau_t_alpha)

# maximise
ELBO_wrap <- function(VIparams, X, y, n, K, G, sum_log_y_fac, rho, tau, tau_alpha) {
  nm <- names(VIparams)
  mu_alpha <- VIparams["mu_alpha"]
  tau_t_alpha <- exp(VIparams["log_tau_t_alpha"])
  u_nm <- sapply(1:G, function(g) paste0("u", g))
  u <- VIparams[u_nm]
  mu_nm <- sapply(1:(G * K), function(i) paste0("mu", i))
  mu <- matrix(VIparams[mu_nm], nrow = G)
  Sigma_nm <- sapply(1:(G * K * K), function(i) paste0("log_Sigma", i))
  Sigma2 <- array(exp(VIparams[Sigma_nm]), dim = c(G, K, K))
  # Sigma_nm <- sapply(1:(G * K), function(i) paste0("log_Sigma", i))
  # Sigma_mat <- matrix(exp(VIparams[Sigma_nm]), nrow = G, ncol = K)
  # Sigma <- array(0, dim = c(G, K, K))
  # for (g in 1:G) {
  #   for (k in 1:K) {
  #     Sigma[g, k, k] <- Sigma_mat[g, k]
  #   }
  # }
  elbo_poisson_maxim(X = X, y = y, n = n, K = K, G = G,
               sum_log_y_fac = sum_log_y_fac, 
               u = u, mu = mu, Sigma = Sigma,
               mu_alpha = mu_alpha, tau_t_alpha = tau_t_alpha,
               tau = tau, rho = rho, tau_alpha = tau_alpha)
}

ELBO_wrap2 <- function(VIparams, X, y, n, K, G, sum_log_y_fac, rho, tau, tau_alpha) {
  nm <- names(VIparams)
  mu_alpha <- VIparams["mu_alpha"]
  tau_t_alpha <- exp(VIparams["log_tau_t_alpha"])
  u_nm <- sapply(1:G, function(g) paste0("u", g))
  u <- VIparams[u_nm]
  prob <- 1 / (1 + exp(-u))
  mu_nm <- sapply(1:(G * K), function(i) paste0("mu", i))
  mu <- matrix(VIparams[mu_nm], nrow = G)
  Sigma_nm <- sapply(1:(G * K), function(i) paste0("log_Sigma", i))
  Sigma <- matrix(exp(VIparams[Sigma_nm]), nrow = G, ncol = K)
  elbo_poisson_maxim2(X = X, y = y, n = n, K = K, G = G,
                     sum_log_y_fac = sum_log_y_fac, 
                     prob = prob, mu = mu, Sigma = Sigma,
                     mu_alpha = mu_alpha, tau_t_alpha = tau_t_alpha,
                     tau = tau, rho = rho, tau_alpha = tau_alpha)
}

ELBO_wrap(VIparams2, X, y, n, K, G, sum_log_y_fac, rho, tau, tau_alpha)
ELBO_wrap2(VIparams2, X, y, n, K, G, sum_log_y_fac, rho, tau, tau_alpha)

mod2 <- optim(par = VIparams, fn = ELBO_wrap, 
              X = X, y = y, n = n, K = K, G = G,
              sum_log_y_fac = sum_log_y_fac, 
              rho = rho, tau = tau, tau_alpha = tau_alpha,
              control = list(maxit = 50000, fnscale = -1))
# fnscale = -1 ensures maximisation rather than minimisation
mod2$convergence

mod4 <- optim(par = VIparams2, fn = ELBO_wrap2, 
              X = X, y = y, n = n, K = K, G = G,
              sum_log_y_fac = sum_log_y_fac, 
              rho = rho, tau = tau, tau_alpha = tau_alpha,
              control = list(maxit = 50000, fnscale = -1))
# fnscale = -1 ensures maximisation rather than minimisation
mod4$convergence

res <- mod2$par
nm <- names(res)
mu_alpha_res <- res["mu_alpha"]
tau_t_alpha_res <- exp(res["log_tau_t_alpha"])
u_nm <- sapply(1:G, function(g) paste0("u", g))
u_res <- res[u_nm]
prob_res <- 1/(1+exp(-u_res))
mu_nm <- sapply(1:(G * K), function(i) paste0("mu", i))
mu_res <- matrix(res[mu_nm], nrow = G)
# Sigma_nm <- sapply(1:(G * K * K), function(i) paste0("log_Sigma", i))
# Sigma_res <- array(exp(res[Sigma_nm]), dim = c(G, K, K))

cbind(as.numeric(mu_res * prob_res), as.numeric(beta))
cbind(prob_res, s_true)
plot(as.numeric(mu_res) * prob_res, as.numeric(beta))
abline(a = 0, b = 1, col = "red")
}

# Plot results -----------------------------------------------------------------------

# ELBO at every time step
ELBO_track <- mod1$ELBO_track
plot_start <- 1
plot_end <- length(ELBO_track)
# plot_end <- 12
plot(plot_start:plot_end,
     ELBO_track[plot_start:plot_end],
     type = "l")
sum(ELBO_track[2:(length(ELBO_track))] - ELBO_track[1:(length(ELBO_track) - 1)] < 0)
which(ELBO_track[2:(length(ELBO_track))] - ELBO_track[1:(length(ELBO_track) - 1)] < 0)
ELBO_track[length(ELBO_track)]

# Compare parameter estimates to truth -----------------------------------------------

# compare estimates to starting values
plot(mod1$vi_params$u, u)

# Compare effect estimates to truth
prob_est <- 1 / (1 + exp(-mod1$vi_params$u))
mu_est <- mod1$vi_params$mu
mu_alpha_est <- mod1$vi_params$mu_alpha
Sigma_est <- mod1$vi_params$Sigma

c(mu_alpha_est, alpha)
tapply(prob_est, s_true, summary)
plot(prob_est, s_true)
plot(mu_est * prob_est, beta)
abline(a = 0, b = 1, col = "red")
cbind(as.numeric(mu_est * prob_est), as.numeric(beta))
cbind(Sigma_est, prob_est)

lm(as.numeric(mu_est) ~ 0 + as.numeric(beta))

X2 <- matrix(nrow = n, ncol = K * G)
for (g in 1:G) {
  for (k in 1:K) {
    X2[ , K * (g - 1) + k] <- X[g, , k] * (prob_est[g] > 0.5)
  }
}
mod3 <- glm(y ~ X2, family = poisson)
cbind(mod3$coefficients, c(mu_alpha_est, as.numeric(t(mu_est))))
plot(mod3$coefficients, c(mu_alpha_est, as.numeric(t(mu_est))))
abline(a = 0, b = 1, col = "red")
plot(mod3$coefficients, c(alpha, as.numeric(t(beta))))
abline(a = 0, b = 1, col = "red")
plot(mod3$coefficients, c(mu_alpha_res, as.numeric(t(mu_res * prob_res))))
abline(a = 0, b = 1, col = "red")

# plot(c(alpha, as.numeric(t(beta))), c(mu_alpha_est, as.numeric(t(mu_est))))
# plot(mod2$coefficients, c(alpha, as.numeric(t(beta))))


