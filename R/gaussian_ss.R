# --------------------------------------------------------------------------------- #
# -------- Spike & slab group variable selection for Gaussian outcome  ------------ #
# --------------------------------------------------------------------------------- #

######### ~.~ Package dependencies ~.~ #########
require(plyr)

######### ~.~ Sub-functions ~.~ #########
trace_prod <- function(A,B) sum(t(A) * B)

log1p.exp <- function(x){
  if(x > 20){
    return(x)
  } else {
    return(log1p(exp(x)))
  }
}

log1p.exp.vec <- function(x){
  y <- x
  which.small <- x <= 20
  y[which.small] <- log1p(exp(x[which.small]))
  return(y)
}

loglogit <- function(x){
  -log1p.exp.vec(-x)
}

######### ~.~ Compute ELBO ~.~ ###########
compute_ELBO_normal <- function(X,XtX,y,prob,mu,Sigma,tau_t,sigma2,rho,tau,update.hyper=F){
  ### -- Pre-calcs -- ###
  
  ## Sum of squared residuals
  lp <- numeric(n) + 0
  for(g in 1:G){
    lp <- lp + prob[g] * matrix(X[g,,],nrow=n) %*% matrix(mu[g,],nrow=K)
  }
  SSR <- sum((y-lp)^2)
  
  ## Expected sum of squared residuals
  prob_tr <- 0
  for(g in 1:G){
    prob_tr <- prob_tr + prob[g]*trace_prod(Sigma[g,,],XtX[g,,])
  }
  E_SSR <- prob_tr + SSR
  
  ## Expected sum of squared gammas
  E_SSgamma <- 0
  for(g in 1:G){
    E_SSgamma <- E_SSgamma + prob[g]*(sum(diag(matrix(Sigma[g,,],nrow=K))) + sum(mu[g,]^2))
  }
  E_SSgamma <- E_SSgamma + K*tau_t*sum(1-prob)
  
  ### -- Compute ELBO -- ###
  
  ## Line 1
  line1 <- -1/(2*sigma2)*E_SSR - (n/2)*log(2*pi*sigma2)
  
  ## Line 2
  line2 <- -E_SSgamma/(2*tau) - K*G*log(2*pi*tau)/2 + log(rho^sum(prob)) + log((1-rho)^sum(1-prob))
  
  ## Line 3
  line3 <- 0.5*(K*sum(prob) + G*K*log(2*pi) + sum(prob * apply(Sigma,1,function(A) log(det(A)))))
  
  ## Line 4
  line4 <- sum(1-prob)*K*0.5*(1 + log(2*pi*tau_t))
  
  ## Line 5
  line5 <- sum(prob[prob!=0]*log(prob[prob!=0])) + sum((1-prob[prob!=1])*log(1-prob[prob!=1]))
  
  ## ELBO
  ELBO <- line1 + line2 + line3 + line4 + line5
  
  ### -- Update hyperparams -- ###
  if(update.hyper){
    sigma2 <- E_SSR/n
    tau <- E_SSgamma/(K*G)
    rho <- mean(prob)
    
    return(list(ELBO=ELBO,sigma2=sigma2,tau=tau,rho=rho))
  } else {
    return(ELBO)
  }
}

######### ~.~ VI updates ~.~ ###########
update_params_normal <- function(X,XtX,y,prob,mu,Sigma,Sigma_inv,tau_t,sigma2,rho,tau,update.hyper=F,update.hyper.last=F){
  
  ### -- Sigma -- ###
  if(update.hyper.last){
    # Sigma only needs to updated if hyperparams were updated last step
    Sigma_inv <- aaply(.data=XtX/sigma2,.margins=1,.fun=function(A,tau,K) A + diag(1/tau,K),tau=tau,K=K,.drop=F)
    Sigma <- aaply(.data=Sigma_inv,.margins=1,.fun=solve,.drop=F)
  }
  
  ### -- mu -- ###
  pred_g <- matrix(0,nrow=n,ncol=G)
  for(g in 2:G){
    pred_g[,g] <- prob[g] * matrix(X[g,,],nrow=n) %*% matrix(mu[g,],nrow=K)
  }
  for(g in 1:G){
    mu[g,] <- (1/sigma2)*matrix(Sigma[g,,],nrow=K) %*% crossprod(matrix(X[g,,],nrow=n), y - rowSums(pred_g[,-g,drop=F]))
    if(g != G){
      pred_g[,g] <- prob[g] * matrix(X[g,,],nrow=n) %*% matrix(mu[g,],nrow=K)
    }
  }
  
  ### -- prob (pi) -- ###
  const <- log(rho/(1-rho)) - 0.5*K*log(tau_t)
  for(g in 1:G){
    u <- crossprod(matrix(mu[g,],nrow=K), matrix(Sigma_inv[g,,],nrow=K)) %*% matrix(mu[g,],nrow=K) + 0.5*log(det(matrix(Sigma[g,,],nrow=K))) + const
    prob[g] <- exp(loglogit(u[1,1]))
  }
  
  ### -- tau_t -- ###
  tau_t <- tau
  
  ### -- Compute ELBO -- ###
  ELBO <- compute_ELBO_normal(X=X,XtX=XtX,y=y,prob=prob,mu=mu,Sigma=Sigma,tau_t=tau_t,
                              sigma2=sigma2,rho=rho,tau=tau,update.hyper=update.hyper)
  if(update.hyper){
    sigma2 <- ELBO$sigma2
    rho <- ELBO$rho
    tau <- ELBO$tau
    ELBO <- ELBO$ELBO
  }
  
  ### -- Return -- ###
  return(list(ELBO=ELBO,prob=prob,mu=mu,Sigma=Sigma,Sigma_inv=Sigma_inv,tau_t=tau_t,
              sigma2=sigma2,rho=rho,tau=tau))
  
}

