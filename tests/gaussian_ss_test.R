# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

devtools::load_all() # Sources all files in R/
# Input parameters ------------------------------------------------------------------
K <- 2
G <- 200 # note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
beta <- matrix(sample(c(-1,1),replace=T,size=G*K,prob=c(0.25,0.75)),nrow=G)
s_true <- sample(c(0,1),replace=T,size=G)
beta <- beta*s_true
sigma2 <- 0.5
n <- 500
# Generate some data -----------------------------------------------------------------
X <- array(rnorm(K*G*n),dim=c(G,n,K))
lp <- numeric(n) + 0
for(g in 1:G){
  lp <- lp + matrix(X[g,,],nrow=n) %*% matrix(beta[g,],nrow=K)
}
lp <- as.numeric(lp)
y <- lp + rnorm(n,sigma2)
# Run algorithm ----------------------------------------------------------------------
plot.start <- 40
plot(plot.start:length(ELBO.track), ELBO.track[plot.start:length(ELBO.track)], type = "l")
plot(plot.start:length(rho.track), rho.track[plot.start:length(rho.track)], type = "l")
plot(plot.start:length(tau.track), tau.track[plot.start:length(tau.track)], type = "l")
plot(plot.start:length(sigma2.track), sigma2.track[plot.start:length(sigma2.track)], type = "l")

tapply(params$prob, s_true, summary)
plot(params$prob, s_true)
plot(params$mu, beta)
