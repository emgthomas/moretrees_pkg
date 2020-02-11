# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian outcomes  ------------------------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #
devtools::load_all() # Sources all files in R/

require(igraph)
require(Matrix)

# Input parameters -------------------------------------------------------------------
group <- "7"
tr <- ccs_tree(group)$tr
leaves <- names(V(tr)[V(tr)$leaf])
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- expm(Matrix::t(A))
A[A > 0 ] <- 1 
G <- length(V(tr))
n <- 3000
K_g <- 2 # number of variables
K <- rep(K_g, G)
m <- 0
tau <- 3
rho <- 0.1
omega <- 2
sigma2 <- 2

# Generate randomly grouped beta (groups follow tree)
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
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
groups_true <- as.integer(as.factor(as.numeric(beta1)))
table(groups_true)

# Generate some data -----------------------------------------------------------------
X <- matrix(rnorm(n * K_g), nrow = n)
outcomes <- sample(leaves, size = n, replace = T)

# Create non-sparse design matrix
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)

# Get linear predictor
lp <- W %*% theta
for (v in leaves) {
  which_v <- outcomes == v
  lp[which_v] <- lp[which_v] + X[which_v, ] %*% beta[v , ]
}
lp <- as.numeric(lp)

# Simulate outcomes
y <- lp + rnorm(n, mean = 0, sd = sqrt(sigma2))

# Run algorithm ----------------------------------------------------------------------
# Create MORETreeS design matrix
dsgn <- moretrees_design_matrix(X = X, tr = tr, y = y, outcomes = outcomes)
mod1 <- spike_and_slab_normal(dsgn$y, dsgn$Xstar, W, update_hyper = T, update_hyper_freq = 50,
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
xi_est <- sapply(1:G, 
                 function(v) as.numeric(mod1$vi_params$mu[[v]] * (mod1$vi_params$prob[v] >= 0.5)),
                 simplify = T)
xi_est1 <- xi_est[1, ]
xi_est2 <- xi_est[2, ]
plot(c(xi_est1, xi_est2), c(xi1, xi2))
abline(a = 0, b = 1, col = "red")
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

# # Compare non-sparse effect estimates to truth --------------------------------------------
# plot(mod1$vi_params$delta, theta)
# abline(a = 0, b = 1, col = "red")

# Compare estimated groups to truth ---------------------------------------------------------
groups_est <- as.integer(as.factor(as.numeric(beta_est1)))
table(groups_est, groups_true)

# Compare moretrees estimates to maximum likelihood -----------------------------------------
groups_list <- sapply(1:max(groups_est), function(i) leaves[groups_est == i])
beta_ml <- matrix(nrow = max(groups_est), ncol = K_g)
beta_est <- matrix(nrow = max(groups_est), ncol = K_g)
beta_true <- matrix(nrow = max(groups_est), ncol = K_g)
for (i in 1:max(groups_est)) {
  which_i <- outcomes %in% groups_list[[i]]
  y_i <- y[which_i]
  y_i[y_i == -1] <- 0
  mod_ml <- lm(y_i ~ 0 + X[which_i, ])
  beta_ml[i, ] <- mod_ml$coefficients
  beta_est[i, ] <- c(unique(beta_est1[groups_est == i]),
                     unique(beta_est2[groups_est == i]))
  beta_true[i, ] <- c(mean(beta1[groups_est == i]),
                      mean(beta2[groups_est == i]))
}
plot(beta_ml, beta_est)
abline(a = 0, b = 1, col = "red")
cbind(beta_ml[ , 1], beta_est[ ,1], beta_true[ , 1])
cbind(beta_ml[ , 2], beta_est[ ,2], beta_true[ , 2])

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

