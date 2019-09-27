# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

source("./scripts/gaussian_ss.R")
require(plyr)

#### Parameters ###
K <- 1
G <- 200 ## note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
beta <- matrix(sample(c(-1,0,1),replace=T,size=G,prob=c(0.25,0.5,0.25)),nrow=G)
sigma2 <- 0.5
n <- 1000

#### Generate some data ####

# generate X
X <- array(rnorm(K*G*n),dim=c(G,n,K))

# get linear predictor
lp <- numeric(n) + 0
for(g in 1:G){
  lp <- lp + matrix(X[g,,],nrow=n) %*% matrix(beta[g,],nrow=K)
}
lp <- as.numeric(lp)

# generate y
y <- lp + rnorm(n,sigma2)

#### Pre-calcs ####
XtX <- aaply(X,1,function(X) crossprod(X,X),.drop=F)

#### Initial values ####
# Hyperparams (for now, assume known)
rho <- 0.5
tau <- var(as.numeric(beta[beta!=0]))

# Variational params
Sigma_inv <- aaply(.data=XtX/sigma2,.margins=1,.fun=function(A,tau,K) A + diag(1/tau,K),tau=tau,K=K,.drop=F)
Sigma <- aaply(.data=Sigma_inv,.margins=1,.fun=solve,.drop=F)
sigma2 <- 2
mu <- matrix(rnorm(G*K,sd=10),nrow=G)
prob <- runif(G)
tau_t <- tau

# Other inputs
ELBO.old <- -1E16
ELBO.new <- 1E16
tol <- 1E-4
update.hyper.freq <- 10

# put in list
params <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)

#### Run algorithm ####

ELBO.track <- c()
i <- 0
while(abs(ELBO.new-ELBO.old)>tol){
  ELBO.old <- ELBO.new
  i <- i+1
  update.hyper <- (i %% update.hyper.freq == 0)
  update.hyper.last <- (i %% update.hyper.freq == 1)
  # update.hyper <- FALSE
  # update.hyper.last <- FALSE
  params <- update_params_normal(ELBO.old=ELBO.old,X=X,XtX=XtX,y=y,
                                 prob=params$prob,mu=params$mu,Sigma=params$Sigma,Sigma_inv=params$Sigma_inv,
                                 tau_t=params$tau_t,sigma2=params$sigma2,rho=params$rho,tau=params$tau,
                                 update.hyper=update.hyper,update.hyper.last=update.hyper.last)
  ELBO.new <- params$ELBO
  ELBO.track <- c(ELBO.track,ELBO.new)
  print(i)
}

plot.start <- 3
plot(plot.start:length(ELBO.track),ELBO.track[plot.start:length(ELBO.track)],type="l")

# prob <- params$prob
# mu <- params$mu
# Sigma <- params$Sigma
# Sigma_inv <- params$Sigma_inv
# tau_t <- params$tau_t
# sigma2 <- params$sigma2
# rho <- params$rho
# tau <- params$tau

plot(params$prob,beta)
plot(params$mu,beta)
