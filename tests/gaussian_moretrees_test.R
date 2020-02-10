# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian outcomes  ------------------------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #
# set.seed(98647)
devtools::load_all() # Sources all files in R/

require(icd)
require(igraph)
require(stringr)
require(Matrix)

n <- 100
K_g <- 2

group <- "7.4"
tr <- ccs_tree("7")$tr
plot.igraph(tr, layout = layout_as_tree, root = group)
X <- matrix(rnorm(n * K_g), nrow = n)
outcomes <- sample(names(V(tr)[V(tr)$leaf]), size = n, replace = T)

# Create MORETreeS design matrix
dsgn <- moretrees_design_matrix(X, tr, outcomes)
Xstar <- dsgn$Xstar
A <- dsgn$A
rm(dsgn)

# Input parameters -------------------------------------------------------------------
G <- length(Xstar)
K <- rep(K_g, G)
m <- 0
tau <- 3
rho <- 0.2
omega <- 2
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
s_true <- rbinom(n = G, size = 1, prob = rho)
s_true[1] <- 1
xi <- sapply(1:G, function(g) Matrix::Matrix(gamma_true[[g]] * s_true[[g]],
             ncol = 1))
xi1 <- sapply(xi, function(xi) xi[1 , 1])
xi2 <- sapply(xi, function(xi) xi[2 , 1])
xi <- cbind(xi1, xi2)
beta1 <- A[leaves, ] %*% xi1
beta2 <- A[leaves, ] %*% xi2
beta <- cbind(beta1, beta2)
row.names(beta) <- leaves
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
sigma2 <- 2

# Generate some data -----------------------------------------------------------------
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (g in 1:G) {
  lp <- lp + Xstar[[g]] %*% xi[g, ]
}
lp <- as.numeric(lp)
y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))
# Run algorithm ----------------------------------------------------------------------
mod1 <- spike_and_slab_normal(y, Xstar, W, update_hyper = T, update_hyper_freq = 50,
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
plot_start <- 502
plot_end <- length(ELBO_track)
# plot_end <- 240
plot(plot_start:plot_end,
     ELBO_track[plot_start:plot_end],
     type = "l")

# Get betas
xi_est <- sapply(1:p, 
    function(v) mod1$vi_params$mu[[v]] * (mod1$vi_params$prob[v] >= 0.5))
xi_est1 <- sapply(xi_est, function(xi) xi[1 , 1])
xi_est2 <- sapply(xi_est, function(xi) xi[2 , 1])
beta_est1 <- A[leaves, ] %*% xi_est1
beta_est2 <- A[leaves, ] %*% xi_est2

# Compare sparse effect estimates to truth --------------------------------------------------
# compare estimated probabilities of variable inclusion to true variable inclusion indicators
tapply(mod1$vi_params$prob, s_true, summary) 
plot(mod1$vi_params$prob, s_true)
# compare estimated coefficients to true coefficients
plot(beta_est1, beta1)
abline(a = 0, b = 1, col = "red")
plot(beta_est2, beta2)
abline(a = 0, b = 1, col = "red")

# Compare non-sparse effect estimates to truth ----------------------------------------------
plot(mod1$vi_params$delta, theta)
abline(a = 0, b = 1, col = "red")

# Compare moretrees estimates to maximum likelihood -----------------------------------------
mod2 <- lm(y ~ 0 + as.numeric(unlist(x_splt[[1]])) + as.numeric(unlist(x_splt[[2]])))
beta_ml <- mod2$coefficients
plot(mod2$coefficients, c(theta, xi_unlist))
abline(a = 0, b = 1, col = "red")
plot(mod2$coefficients, c(as.numeric(mod1$vi_params$delta),
                          moretrees_est))
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

