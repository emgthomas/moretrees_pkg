# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

source("./R/gaussian_ss.R")
require(plyr)

#### Parameters ###
K <- 1
G <- 200 ## note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
beta <- matrix(sample(c(-1,0,1),replace=T,size=G,prob=c(0.25,0.5,0.25)),nrow=G)
sigma2 <- 0.5
n <- 5000

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
attributes(XtX)$dimnames <- NULL

#### Initial values ####
# Hyperparams (for now, assume known)
rho <- 0.5
tau <- var(as.numeric(beta[beta!=0]))

# Variational params
Sigma_inv <- aaply(.data=XtX,.margins=1,.fun=function(XtX,tau,K,sigma2) XtX/sigma2 + diag(1/tau,K),tau=tau,K=K,sigma2=sigma2,.drop=F)
attributes(Sigma_inv)$dimnames <- NULL
Sigma_inv.check <- array(0,dim=c(G,K,K))
for(g in 1:G){
  Sigma_inv.check[g,,] <- XtX[g,,]/sigma2 +  diag(1/tau,K)
}
all.equal(Sigma_inv,Sigma_inv.check)
Sigma <- aaply(.data=Sigma_inv,.margins=1,.fun=solve,.drop=F)
attributes(Sigma)$dimnames <- NULL
Sigma.check <- array(0,dim=c(G,K,K))
for(g in 1:G){
  Sigma.check[g,,] <- solve(matrix(Sigma_inv[g,,],nrow=K))
}
all.equal(Sigma,Sigma.check)

mu <- matrix(rnorm(G*K,sd=10),nrow=G)
prob <- runif(G)
tau_t <- tau

# Other inputs
ELBO.old <- -1E16
ELBO.new <- 1E16
tol <- 1E-4
update.hyper.freq <- 10

# put in list
# params.first <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)
# params <- params.first
params.second <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)
params <- params.second
# print(params)

#### Run algorithm ####

ELBO.track <- c()
rho.track <- c()
tau.track <- c()
i <- 0
while(abs(ELBO.new-ELBO.old)>tol){
  ELBO.old <- ELBO.new
  i <- i+1
  update.hyper <- (i %% update.hyper.freq == 0)
  update.hyper.last <- (i %% update.hyper.freq == 1)
  # update.hyper <- FALSE
  # update.hyper.last <- FALSE
  params <- update_params_normal(X=X,XtX=XtX,y=y,
                                 prob=params$prob,mu=params$mu,Sigma=params$Sigma,Sigma_inv=params$Sigma_inv,
                                 tau_t=params$tau_t,sigma2=params$sigma2,rho=params$rho,tau=params$tau,
                                 update.hyper=update.hyper,update.hyper.last=update.hyper.last)
  ELBO.new <- params$ELBO
  ELBO.track <- c(ELBO.track,ELBO.new)
  rho.track <- c(rho.track,params$rho)
  tau.track <- c(tau.track,params$tau)
  print(i)
}

plot.start <- 15
plot(plot.start:length(ELBO.track),ELBO.track[plot.start:length(ELBO.track)],type="l")
plot(plot.start:length(rho.track),rho.track[plot.start:length(rho.track)],type="l")
plot(plot.start:length(tau.track),tau.track[plot.start:length(tau.track)],type="l")

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
