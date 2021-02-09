## EM algorithm for a mixture of multivariate Gaussians

## we need this to calculate the density of a multivariate normal
library(mvtnorm)

## logsumexp(x) returns log(sum(exp(x))) but performs the computation in a more stable manner
logsumexp <- function(x) return(log(sum(exp(x - max(x)))) + max(x))

# compute the log-likelihood using expression ell(q,theta)
compute_ll <- function(xs,mus,Sigmas,lws,gammas) {
  ll <- 0
  n <- dim(xs)[1]
  K <- dim(mus)[1]
  for (i in 1:n) {
    for (k in 1:K) {
      if (gammas[i,k] > 0) {
        ll <- ll + gammas[i,k]*(lws[k]+dmvnorm(xs[i,],mus[k,],Sigmas[[k]],log=TRUE)-log(gammas[i,k]))
      }
    }
  }
  return(ll)
}

# compute the log-likelihood directly
compute_ll.direct <- function(xs,mus,Sigmas,lws) {
  ll <- 0
  n <- dim(xs)[1]
  K <- dim(mus)[1]
  for (i in 1:n) {
    s <- 0
    for (k in 1:K) {
      s <- s + exp(lws[k])*dmvnorm(xs[i,],mus[k,],Sigmas[[k]])
    }
    ll <- ll + log(s)
  }
  return(ll)
}

em_mix_gaussian <- function(xs,K,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  # start off with ws all equal
  lws <- rep(log(1/K),K)
  
  ## start off with a random selection of cluster means
  mus <- xs[sample(n,K),]
  
  ## let the initial covariances be simple estimates of the variance
  vs <- rep(0,p)
  for (j in 1:p) vs[j] <- var(xs[,j])
  Sigmas <- list()
  for (k in 1:K) Sigmas[[k]] <- diag(vs)
  
  # gammas will be set in the first iteration 
  gammas <- matrix(0,n,K)
  
  converged <- FALSE
  numit <- 0
  ll <- -Inf
  print("iteration : log-likelihood")
  while(!converged && numit < max.numit) {
    numit <- numit + 1
    mus.old <- mus
    ll.old <- ll
    
    ## E step - calculate gammas
    for (i in 1:n) {
      # the elements of lprs are log(w_k * p_k(x)) for each k in {1,...K}
      lprs <- rep(0,K)
      for (k in 1:K) {
        lprs[k] <- lws[k] + dmvnorm(xs[i,],mus[k,],Sigmas[[k]],log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,Sigmas,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,Sigmas,lws)
    
    # M step - update ws, mus and Sigmas
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      Sigmas[[k]] <- matrix(0,p,p)
      
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
        Sigmas[[k]] <- Sigmas[[k]] + gammas[i,k]/Ns[k]*xs[i,]%*%t(xs[i,])
      }
      Sigmas[[k]] <- Sigmas[[k]] - mus[k,]%*%t(mus[k,])
    }
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,Sigmas=Sigmas,gammas=gammas,ll=ll))
}
