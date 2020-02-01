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

# Get data.frame showing mapping from ICD9 to multilevel CCS
ccs_icd9 <- data.frame(icd9 = unlist(icd9_map_multi_ccs[[1]]), stringsAsFactors = F)
for (i in 1:4) {
  ccs_list <- icd9_map_multi_ccs[[i]]
  ccs_df <- data.frame(icd9 = unlist(ccs_list), 
                       stringsAsFactors = F)
  ccs_df$ccs <- ccs_list %>% 
    names %>% # names of the list entries are the CCS codes
    sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
    unlist
  names(ccs_df)[2] <- paste0("ccs_l", i)
  ccs_icd9 <- merge(ccs_icd9, ccs_df, by = "icd9", all.x = T, all.y = F)
}

# Keep only diseases of the circulatory system
ccs_icd9 <- subset(ccs_icd9, ccs_l1 == "7")

# CCS codes only
ccs <- ccs_icd9
ccs$icd9 <- NULL
ccs <- ccs[!duplicated(ccs), ]

# Order CCS codes appropriately
ccs$ccs_l4 <- sapply(1:nrow(ccs),
      function(i) if(ccs$ccs_l4[i] == " ") paste0(ccs$ccs_l3[i],".0") else ccs$ccs_l4[i])
ccs_levels <- matrix(nrow = nrow(ccs), ncol = 4)
for (i in 1:nrow(ccs_levels)) {
  ccs_levels[i, ] <- as.integer(str_split(ccs$ccs_l4[i], "\\.")[[1]])
}
ccs_levels <- as.data.frame(ccs_levels)
names(ccs_levels) <- c("l1", "l2", "l3", "l4")
ccs_levels <- cbind(ccs_levels, ccs)
ccs_levels <- ccs_levels[order(ccs_levels$l1, ccs_levels$l2,
                               ccs_levels$l3, ccs_levels$l4) , ]
ccs <- ccs_levels[ , names(ccs)]

# Make tree
edges <- rbind(as.matrix(ccs[ , c(1, 2)]),
               as.matrix(ccs[ , c(2, 3)]),
               as.matrix(ccs[ , c(3, 4)]))
edges <- edges[edges[ , 2] != " ", ]
edges <- edges[!duplicated(edges), ]
tr <- igraph::graph_from_edgelist(e = edges, directed = T)
igraph::plot.igraph(tr, layout = layout_as_tree, root = 1)
leaves <- igraph::V(tr)[igraph::degree(tr, mode = "out") == 0]
igraph::V(tr)$leaf <- FALSE
igraph::V(tr)$leaf[V(tr) %in% leaves] <- TRUE
ccs_levels_sub <- subset(ccs_levels, l2 == 4)
v_sub <- unique(c(ccs_levels_sub$ccs_l2,
                  ccs_levels_sub$ccs_l3,
                  ccs_levels_sub$ccs_l4))
v_sub <- v_sub[!(v_sub == " ")]
subtr <- igraph::induced_subgraph(graph = tr,
                          v = v_sub)
plot.igraph(subtr, layout = layout_as_tree, root = "7.4")

# Extract relevant parameters for analysis
D <- igraph::as_adjacency_matrix(subtr, sparse = T)
A <- expm(Matrix::t(D))
A[A > 0 ] <- 1 
p <- length(igraph::V(subtr))
leaves <- V(subtr)$leaf
pL <- sum(leaves)
codes <- names(V(subtr))

# Construct design matrix
K_g <- 2
splt <- 0
n_leaf <- rep(100, pL)
n <- sum(n_leaf)
x <- sapply(n_leaf, rnorm, simplify = F) # exposure for each outcome
x_splt <- list()
x_splt[[1]] <- sapply(x, function(x, splt) {
  x[x > splt] <- 0
  x
}, splt = splt, simplify = F)
x_splt[[2]] <- sapply(x, function(x, splt) {
  x[x <= splt] <- 0
  x
}, splt = splt, simplify = F)
Xmat <- list()
Xstar <- list()
for (k in 1:K_g) {
  Xmat_k <- Matrix(0, nrow = n, ncol = p)
  Xmat_k[ , leaves] <- Matrix::bdiag(x_splt[[k]])
  colnames(Xmat_k) <- codes
  Xmat[[k]] <- Xmat_k
  Xstar[[k]] <- Xmat[[k]] %*% A
}

# Generate data
X <- rep(list(), p)
for (v in 1:p) {
  Xmat_v <- Matrix::Matrix(0, nrow = n, ncol = K_g)
  for (k in 1:K_g) {
    Xmat_v[ , k] <- Xstar[[k]][ , v]
  }
  X[[v]] <- Xmat_v
}

# Input parameters -------------------------------------------------------------------
G <- p
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
beta1 <- A[leaves, ] %*% xi1
beta2 <- A[leaves, ] %*% xi2
beta <- list(beta1, beta2)
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
sigma2 <- 2

# Generate some data -----------------------------------------------------------------
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (k in 1:K_g) {
  lp <- lp + Matrix::bdiag(x_splt[[k]]) %*% beta[[k]]
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
plot_start <- 100
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

