# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian and Bernoulli outcomes  ----------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

rm(list = ls(), inherits = T)
devtools::load_all() # Sources all files in R/

# Parameter choices to be tested -----------------------------------------------------
params_list <- list(K_g = 1:2, 
                    m = 0:2,
                    nrestarts = c(1, 2))
params <- do.call(expand.grid, params_list)
i <- 12
params[i, ]

# Input parameters -------------------------------------------------------------------
nrestarts <- params$nrestarts[i]
n <- 1E3
group <- "7.4"
tr <- ccs_tree(group)$tr
leaves <- names(igraph::V(tr)[igraph::V(tr)$leaf])
# If desired, specify levels 
igraph::V(tr)$levels <- 1
igraph::V(tr)$levels[igraph::V(tr)$leaf] <- 2
A <- igraph::as_adjacency_matrix(tr, sparse = T)
A <- Matrix::expm(Matrix::t(A))
A[A > 0 ] <- 1
G <- length(igraph::V(tr))
p <- G
pL <- sum(igraph::V(tr)$leaf)
K_g <- params$K_g[i] # number of variables
K <- rep(K_g, G)
m <- params$m[i]
tau <- 3
hyper_fixed <- list(a_rho = c(0.9, 0.5), b_rho = c(3 , 2))
rho1 <- rbeta(1, shape1 = hyper_fixed$a_rho[1], shape2 = hyper_fixed$b_rho[1]) # rho for internal nodes
rho2 <- rbeta(1, shape1 = hyper_fixed$a_rho[2], shape2 = hyper_fixed$b_rho[2]) # rho for leaf nodes
omega <- 2
doParallel::registerDoParallel(cores = nrestarts)
log_dir <- "./tests/"

# Generate randomly grouped beta (groups follow tree)
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau), simplify = F)
s_true <- c(1, rbinom(n = p - pL - 1, size = 1, prob = rho1), 
            rbinom(n = pL, size = 1, prob = rho2))
xi <- mapply(function(gamma, s) matrix(gamma * s, nrow = 1),
             gamma = gamma_true, s = s_true,
             SIMPLIFY = T) %>% t
if (K_g == 1) xi <- t(xi)
beta <- A[leaves, ] %*% xi
zeta <- matrix(rnorm(m * G, mean = 0, sd = sqrt(omega)), nrow = G, ncol = m)
theta <- A[leaves, ] %*% zeta
groups_true <- as.integer(as.factor(as.numeric(beta[ , 1])))
table(groups_true)

# Generate some data -----------------------------------------------------------------
X <- matrix(rnorm(n * K_g), nrow = n, ncol = K_g)
outcomes <- c(leaves, sample(leaves, size = n - pL, replace = T))

# Create non-sparse design matrix
if (m > 0) {
  W <- matrix(rnorm(m * n, sd = 0.5), nrow = n, ncol = m)
} else {
  W <- NULL
}

# Get linear predictor
lp <- numeric(n)
for (v in leaves) {
  which_v <- outcomes == v
  lp[which_v] <- lp[which_v] + X[which_v, , drop = F] %*% Matrix::t(beta[v, , drop = F])
  if (m > 0) {
    lp[which_v] <- lp[which_v] + W[which_v, , drop = F] %*% Matrix::t(theta[v, , drop = F])
  }
}

# Simulate outcomes
p_success <- expit(lp)
y <- runif(n)
y <- as.integer(y <= p_success)

# Run algorithm ----------------------------------------------------------------------
require(gdata)
keep(X, W, y, outcomes, tr, nrestarts, hyper_fixed,
    s_true, groups_true, beta, theta, sure = T)
# Run model without W
mod_start <- moretrees(X = X, W = NULL, y = y, outcomes = outcomes,
                       tr = tr,
                       update_hyper_freq = 50,
                       hyper_fixed = hyper_fixed,
                       tol = 1E-8, 
                       tol_hyper = 1E-4,
                       max_iter = 3E4,
                       print_freq = 50,
                       nrestarts = nrestarts,
                       parallel = F,
                       get_ml = T,
                       log_dir = "./tests/")
mod_end <- moretrees(X = X, W = W, y = y, outcomes = outcomes,
                     initial_values = mod_start,
                     tr = tr,
                     update_hyper_freq = 50,
                     hyper_fixed = hyper_fixed,
                     tol = 1E-8, 
                     tol_hyper = 1E-4,
                     max_iter = 3E4,
                     print_freq = 50,
                     nrestarts = nrestarts,
                     get_ml = T,
                     log_dir = log_dir)
beta_est <- mod_end$beta_est
beta_moretrees <- mod_end$beta_moretrees
beta_ml <- mod_end$beta_ml
theta_est <- mod_end$theta_est
mod_restarts <- mod_end$mod_restarts
mod1 <- mod_end$mod

# Compare ELBOs for random restarts --------------------------------------------------
c(mod1$ELBO_track[length(mod1$ELBO_track)],
  sapply(mod_restarts, function(mod) mod$ELBO_track[length(mod$ELBO_track)]))
for (i in 1:nrestarts) {
  fl <- paste0(log_dir, "restart_", i, "_log.txt")
  if (file.exists(fl)) {
    file.remove(fl)
  }
}

# Plot results -----------------------------------------------------------------------

# Check if the ELBO decreases
# ELBO_track <- mod1$ELBO_track
ELBO_track <- c(mod_start$mod$ELBO_track, mod_end$mod$ELBO_track[2:length(mod_end$mod$ELBO_track)])
if(min(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)]) < 0) {
  print("ELBO decreases at these time points:")
  which(ELBO_track[2:length(ELBO_track)] - ELBO_track[1:(length(ELBO_track)-1)] < 0)
} else {
  print("ELBO always increases!")
}

# ELBO at every time step
plot_start <- 20000
plot_end <- length(ELBO_track)
# plot_end <- 2519
plot(ELBO_track[plot_start:plot_end],
     type = "l")

# Compare sparse effect estimates to truth --------------------------------------------------
# compare estimated probabilities of variable inclusion to true variable inclusion indicators
tapply(mod1$vi_params$prob, s_true, summary)
plot(mod1$vi_params$prob, s_true)
# compare estimated coefficients to true coefficients
k <- 1
clmn <- paste0("est", k)
plot(beta_est[ , clmn <- paste0("est", k)], beta[ , k])
abline(a = 0, b = 1, col = "red")
if (length(mod1$vi_params$delta[[1]]) > 0) {
  j <- 1
  clmn <- paste0("est", j)
  plot(theta_est[, clmn], theta[ , j])
  abline(a = 0, b = 1, col = "red") 
}

# Compare estimated groups to truth ---------------------------------------------------------
groups_est <- beta_est$group
table(groups_est, groups_true)

# Compare moretrees estimates to maximum likelihood -----------------------------------------
k <- 1
clmn <- paste0("est", k)
plot(beta_ml[, clmn], beta_moretrees[ , clmn])
abline(a = 0, b = 1, col = "red")
# cbind(beta_ml[, clmn], beta_moretrees[ , clmn])


