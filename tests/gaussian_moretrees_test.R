# --------------------------------------------------------------------------------- #
# -------- ssMOReTreeS for Gaussian outcomes  ------------------------------------- #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

require(icd)
require(igraph)

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
ccs_l4 <- sapply(1:nrow(ccs),
      function(i) if(ccs$ccs_l4[i] == " ") paste0(ccs$ccs_l3[i],".0") else ccs$ccs_l4[i])
ccs_levels <- matrix(nrow = nrow(ccs), ncol = 4)
for (i in 1:nrow(ccs_levels)) {
  ccs_levels[i, ] <- as.integer(str_split(ccs_l4[i], "\\.")[[1]])
}
ccs_levels <- as.data.frame(ccs_levels)
names(ccs_levels) <- c("l1", "l2", "l3", "l4")
ccs_levels <- cbind(ccs_levels, ccs)
ccs_levels <- ccs_levels[order(ccs_levels$l1, ccs_levels$l2,
                               ccs_levels$l3, ccs_levels$l4) , ]
ccs <- ccs_levels[ , names(ccs)]

# Make adjacency matrix
codes_all <- unique(unlist(ccs))
codes_all <- codes_all[!(codes_all == " ")]
p <- length(codes_all)
D <- Matrix(0, nrow = p, ncol = p)
rownames(D) <- codes_all
colnames(D) <- codes_all
for (l in 1:3) {
  lvl <- paste0("ccs_l",l)
  lvl_blw <- paste0("ccs_l",l + 1)
  codes_l <- unique(ccs[ , lvl])
  for (v in codes_l) {
    chldrn <- unique(ccs[ccs[ , lvl] == v, lvl_blw])
    chldrn <- chldrn[!(chldrn == " ")]
    if (length(chldrn) > 0) {
      D[rownames(D) == v, colnames(D) %in% chldrn] <- 1
    }
  }
}


# set.seed(98647)
devtools::load_all() # Sources all files in R/
# Input parameters -------------------------------------------------------------------
G <- 20 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
K <- sample(1:4, size = G, replace = T)
m <- 5
tau <- 3
rho <- 0.5
omega <- 2
gamma_true <- sapply(K, rnorm, mean = 0, sd = sqrt(tau))
s_true <- rbinom(n = G, size = 1, prob = rho)
beta <- sapply(1:G, function(g) gamma_true[[g]] * s_true[[g]])
theta <- rnorm(m, mean = 0, sd = sqrt(omega))
sigma2 <- 2
n <- 300
# Generate some data -----------------------------------------------------------------
X <- sapply(K, FUN = function(k) Matrix::Matrix(rnorm(k * n, sd = 0.5), nrow = n))
W <- Matrix::Matrix(rnorm(m * n, sd = 0.5), nrow = n)
lp <- W %*% theta
for (g in 1:G) {
  lp <- lp + X[[g]] %*% beta[[g]]
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

# Compare moretrees estimates to maximum likelihood -----------------------------------------
X1 <- X[[1]]
for (g in 2:G) {
  X1 <- cbind(X1, X[[g]])
}
mod2 <- lm(y ~ 0 + as.matrix(W) + as.matrix(X1))
# mod2 <- lm(y ~ 0 + as.matrix(X1))
plot(mod2$coefficients, c(theta, unlist(beta)))
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

