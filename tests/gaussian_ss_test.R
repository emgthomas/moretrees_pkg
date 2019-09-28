# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# -------- Test code -------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #

source("./R/gaussian_ss.R")
require(plyr)

#### Parameters ###
K <- 3
G <- 2 ## note: for matrices/arrays indexed by g=1,...,G, g is always the first dimension
beta1 <- c(2,1,-1)
beta2 <- c(0,0,0)
beta <- rbind(beta1,beta2)
sigma2 <- 0.5
n <- 10000

#### Generate some data ####

# generate X
x1 <- matrix(data=0,nrow=n,ncol=K)
x1[,1] <- rep(1,n) # intercept
x1[,2] <- seq(-2,2,length.out=n)
x1[,3] <- x1[,2]^2
x2 <- matrix(data=0,nrow=n,ncol=K)
for(i in 1:K){
  x2[,i] <- rnorm(n,runif(1,-2,2)) + sample(-3:3,1)
}
X <- array(c(x1,x2),dim=c(n,K,G))
X <- aperm(X,c(3,1,2))

# get linear predictor
lp <- numeric(n) + 0
for(g in 1:G){
  lp <- lp + X[g,,] %*% beta[g,]
}
lp <- as.numeric(lp)

# generate y
y <- lp + rnorm(n,sigma2)

#### Pre-calcs ####
XtX <- aaply(X,1,function(X) crossprod(X,X))

#### Initial values ####
# Hyperparams (for now, assume known)
rho <- 0.5
tau <- var(as.numeric(beta[beta!=0]))

# Variational params
Sigma_inv <- aaply(.data=XtX/sigma2,.margins=1,.fun=function(A,tau,K) A + diag(1/tau,K),tau=tau,K=K)
Sigma <- aaply(.data=Sigma_inv,.margins=1,.fun=solve)
mu <- matrix(rnorm(G*K),nrow=G)
prob <- rep(0.5,G)
tau_t <- tau

# Other inputs
ELBO.old <- -1E16
ELBO.new <- 1E16
tol <- 1E-4
update.hyper.freq <- 10

# put in list
params <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)
# params.first <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)
# params <- params.first
params.second <- list(mu=mu,prob=prob,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,rho=rho,sigma2=sigma2,tau=tau)
params <- params.second

#### Run algorithm ####

ELBO.track <- c()
i <- 0
while(abs(ELBO.new-ELBO.old)>tol){
  ELBO.old <- ELBO.new
  i <- i+1
  # update.hyper <- (i %% update.hyper.freq == 0)
  # update.hyper.last <- (i %% update.hyper.freq == 1)
  update.hyper <- FALSE
  update.hyper.last <- FALSE
  params <- update_params_normal(X=X,XtX=XtX,y=y,
                                 prob=params$prob,mu=params$mu,Sigma=params$Sigma,Sigma_inv=params$Sigma_inv,
                                 tau_t=params$tau_t,sigma2=params$sigma2,rho=params$rho,tau=params$tau,
                                 update.hyper=update.hyper,update.hyper.last=update.hyper.last)
  ELBO.new <- params$ELBO
  ELBO.track <- c(ELBO.track,ELBO.new)
  print(i)
}

plot.start <- 1
plot(plot.start:length(ELBO.track),ELBO.track[plot.start:length(ELBO.track)],type="l")

prob <- params$prob
mu <- params$mu
Sigma <- params$Sigma
Sigma_inv <- params$Sigma_inv
tau_t <- params$tau_t
sigma2 <- params$sigma2
rho <- params$rho
tau <- params$tau


