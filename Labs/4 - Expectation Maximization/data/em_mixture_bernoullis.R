## EM algorithm for a mixture of Bernoullis

## logsumexp(x) returns log(sum(exp(x))) but performs the computation in a more stable manner
logsumexp <- function(x) return(log(sum(exp(x - max(x)))) + max(x))

prob <- function(x,mu,return.log=FALSE) {
  l <- sum(log(mu[x==1]))+sum(log(1-mu[x==0]))
  if (return.log) {
    return(l)
  } else {
    return(exp(l))
  }
}

compute_ll <- function(xs,mus,lws,gammas) {
  ll <- 0
  n <- dim(xs)[1]
  K <- dim(mus)[1]
  for (i in 1:n) {
    for (k in 1:K) {
      if (gammas[i,k] > 0) {
        ll <- ll + gammas[i,k]*(lws[k]+prob(xs[i,],mus[k,],return.log=TRUE)-log(gammas[i,k]))
      }
    }
  }
  return(ll)
}

compute_ll.direct <- function(xs,mus,lws) {
  ll <- 0
  n <- dim(xs)[1]
  K <- dim(mus)[1]
  for (i in 1:n) {
    s <- 0
    for (k in 1:K) {
      s <- s + exp(lws[k])*prob(xs[i,],mus[k,])
    }
    ll <- ll + log(s)
  }
  return(ll)
}

em_mix_bernoulli <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  # start off with ws all equal
  lws <- rep(log(1/K),K)
  
  if (is.null(start)) {
    mus <- .2 + .6*xs[sample(n,K),]
  } else {
    mus <- start
  }
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
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,lws)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # to avoid a numerical issue since each element of mus must be in [0,1]
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll))
}



em_mix_bernoulli_mu_4_6 <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  # start off with ws all equal
  lws <- rep(log(1/K),K)
  
  if (is.null(start)) {
    mus <- .4 + .2*xs[sample(n,K),]
  } else {
    mus <- start
  }
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
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,lws)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # to avoid a numerical issue since each element of mus must be in [0,1]
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll))
}









em_mix_bernoulli_W_prop <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  lws <- c(log(0.2835242),log(0.2166605), log(0.1635882), log(0.3362271))
  
  if (is.null(start)) {
    mus <- .2 + .6*xs[sample(n,K),]
  } else {
    mus <- start
  }
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
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,lws)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # to avoid a numerical issue since each element of mus must be in [0,1]
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll))
}
























em_mix_bernoulli_test <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  lws <- c(log(0.2835242),log(0.2166605), log(0.1635882), log(0.3362271))
  
  if (is.null(start)) {
    mus <- 1*xs[sample(n,K),]
  } else {
    mus <- start
  }
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
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,lws)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # to avoid a numerical issue since each element of mus must be in [0,1]
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll))
}















em_mix_bernoulli_W_mu_prop <- function(xs,K,start=NULL,max.numit=Inf) {
  p <- dim(xs)[2]
  n <- dim(xs)[1]
  
  # lws is log(ws)
  # we work with logs to keep the numbers stable
  lws <- c(log(0.2835242),log(0.2166605), log(0.1635882), log(0.3362271))
  
  if (is.null(start)) {
    mus <- mus_prop 
  } else {
    mus <- start
  }
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
        lprs[k] <- lws[k] + prob(xs[i,],mus[k,],return.log=TRUE)
      }
      # gammas[i,k] = w_k * p_k(x) / sum_j {w_j * p_j(x)}
      gammas[i,] <- exp(lprs - logsumexp(lprs))
    }
    
    ll <- compute_ll(xs,mus,lws,gammas)
    # we could also compute the log-likelihood directly below
    # ll <- compute_ll.direct(xs,mus,lws)
    
    # M step - update ws and mus
    Ns <- rep(0,K)
    for (k in 1:K) {
      Ns[k] <- sum(gammas[,k])
      lws[k] <- log(Ns[k])-log(n)
      
      mus[k,] <- rep(0,p)
      for (i in 1:n) {
        mus[k,] <- mus[k,]+gammas[i,k]/Ns[k]*xs[i,]
      }
    }
    # to avoid a numerical issue since each element of mus must be in [0,1]
    mus[which(mus > 1,arr.ind=TRUE)] <- 1 - 1e-15
    print(paste(numit,": ",ll))
    # we stop once the increase in the log-likelihood is "small enough"
    if (abs(ll-ll.old) < 1e-5) converged <- TRUE
  }
  return(list(lws=lws,mus=mus,gammas=gammas,ll=ll))
}

