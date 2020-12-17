#############################################
####    Function Definitions    #############
initial.state <- function(u, alpha)
{
  # u is the first observed uniform variate
  # alpha is the initial distribution
  m <- length(alpha)
  I <- matrix(NA, nrow = m, ncol = 2)
  I[1,1] <- 0
  for(i in 1:(m - 1))
  {
    I[i,2] <- I[i,1] + alpha[i]
    if(u < I[i,2])
    {
      Z_1 <- i
      return(Z_1)
    }
    I[i + 1,1] <- I[i,2]
  }
  return(m)
}
######################################

generate.mc <- function(x,alpha,PI) 
{
  # x is the vector of uniform observations
  # alpha is the initial distribution
  # PI is the transition probability
  Z <- c()
  m <- length(alpha)
  n <- length(x)
  Z[1] <- initial.state(x[1], alpha)
  
  for(i in 2:n)
    Z[i] <- initial.state(x[i], PI[Z[i-1],])
  return(Z)
}


###################Forward and backward recursion functions

forward_recursion <- function(numObs = N, numStates = K, initDist = c(1,rep(0,K-1)), a = NULL, b = NULL)
  {
  ct <- rep(1,numObs)
  fvar <- matrix(nrow = numObs,ncol = numStates)
  fvar_star <- matrix(nrow = numObs,ncol = numStates)
  fvar_tilde <- matrix(nrow = numObs,ncol = numStates)
  fvar[1,] <- initDist*b[1,]
  ct[1] <- 1/sum(fvar[1,])
  fvar_tilde[1,] <- fvar[1,]*ct[1]
  
  for(t in 2:numObs){
    for(j in 1:numStates){
      fvar[t,j] <- sum(fvar[t-1,]*a[,j])*b[t,j]
      fvar_star[t,j] <- sum(fvar_tilde[t-1,]*a[,j])*b[t,j]
      #print(t)
    }
    ct[t] <- 1/sum(fvar_star[t,])
    fvar_tilde[t,] <- fvar_star[t,]*ct[t]
  }
  output <- list(fvar_tilde, ct)
  return(output)
}


backward_recursion <- function(numObs = N, numStates = K, a = NULL, b = NULL, ct = NULL){
  bvar_star <- matrix(nrow = numObs,ncol = K)
  bvar_tilde <- matrix(nrow = numObs,ncol = K)
  bvar_tilde[numObs,] <- ct[numObs]
  for(t in (numObs-1):1){
    for(j in 1:numStates){
      bvar_star[t,j] <- sum(a[j,]*bvar_tilde[t+1,]*b[t+1,])
    }
    bvar_tilde[t,] <- bvar_star[t,]*ct[t]
  }
  return(bvar_tilde)
}

post_param = function(numStates = K, numMix = P, gamma.post = gamma_jp, delta.post = delta_jp, zeta = zeta_j, alpha = alpha_j, xi = xi_j)
{
  tmat.post <- matrix(nrow = numStates,ncol = numStates)
  lambda.post <- gamma.post/delta.post
  zeta.post <- matrix(nrow = numStates,ncol = numMix)
      for(j in 1:numStates){
      zeta.post[j,] <- zeta[j,]/sum(zeta[j,])
      tmat.post[j,] <- alpha[j,]/sum(alpha[j,])
      }
  pi.post <- xi/sum(pi)
  
  output <- list('InitDist' = pi.post, 'TransMat' = tmat.post, 'MixProb' = zeta.post, 'RainRate' = lambda.post)
  return(output)
}

simulate.hmm.uni <- function(numDays, numYears=1,numStates, numMix, initDist,tmat,zeta,lambda){
  v.sim <- matrix(nrow = numDays,ncol = numYears)
  duration <- numDays*numYears
  for(year in 1:numYears)
  {
    u <- runif(numDays)
    v.sim[,year] <- generate.mc(u,initDist,tmat)
  }
  v.sim.vector <- as.vector(v.sim)
  mix.unif <- matrix(runif(duration),nrow = duration)
  y.sim=rep(NA,duration)
  for(t in 1:duration){
    which.mixture <- initial.state(u=mix.unif[t],alpha = zeta[v.sim.vector[t],])
    y.sim[t] <- ifelse(which.mixture==1,0,rexp(1,lambda[v.sim.vector[t],which.mixture-1]))
  }
  return(y.sim)
}

simulate.hmm.mult <- function(numDays, numYears=1, numLoc=1, numStates, numMix, initDist,tmat,zeta,lambda){
  v.sim <- matrix(nrow = numDays,ncol = numYears)
  duration <- numDays*numYears
  for(year in 1:numYears)
  {
    u <- runif(numDays)
    v.sim[,year] <- generate.mc(u,initDist,tmat)
  }
  v.sim.vector <- as.vector(v.sim)
  mix.unif <- matrix(runif(duration),nrow = duration)
  y.sim=rep(NA,duration)
  for(t in 1:duration){
    which.mixture <- initial.state(u=mix.unif[t],alpha = zeta[v.sim.vector[t],])
    y.sim[t] <- ifelse(which.mixture==1,0,rexp(1,lambda[v.sim.vector[t],which.mixture-1]))
  }
  return(y.sim)
}

# DIC and ELBO calculations
uDIC <- function(numStates, numMix, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate){
   pd0 <- 0
   pd1 <- 0
   pd2 <- 0
  for(j in 1:numStates)
  {
     pd0 <- pd0 + sum(stateProb[,j,]*mixProb[,j,1,])*(log(zeta[j,1]) - log(sum(zeta[j,])) - 
                digamma(zeta[j,1]) + digamma(sum(zeta[j,])) )

    for(p in 2:numMix)
    {
      pd2 <- pd2 + sum(stateProb[,j,]*mixProb[,j,p,])*( log(gamma_shape[j,p-1]) - digamma(gamma_shape[j,p-1]) +
                   log(zeta[j,p]) - log(sum(zeta[j,])) - digamma(zeta[j,p]) + digamma(sum(zeta[j,])) )
     
    }
    for(k in 1:numStates)
    pd1 <- pd1 + rowSums(jtTransMat,dims = 2)[j,k]*( log(alpha[j,k]) - log(sum(alpha[j,])) - digamma(alpha[j,k]) + digamma(sum(alpha[j,])) )
   
  }
  pd3 <- sum(initProb*(log(xi) - log(sum(xi)) - digamma(xi) + digamma(sum(xi))))
  pd <- pd0+sum(pd1)+pd2+pd3  
  dic <- 4*(pd)  + 2*sum(log(ct))
  output = list('pd' = pd, 'dic' = dic)
  return(output)
}


uELBO = function(numStates, numMix, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate, obs, h){
  kl_C <- 0
  kl_A <- 0
  kl_theta <- 0
  kl_pi <- 0
  for(j in 1:numStates){
    kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,1,])*(digamma(zeta[j,1]) - digamma(sum(zeta[j,])) ) + 
      lgamma(sum(zeta[j,])) - sum(lgamma(zeta[j,])) #-
    #lgamma(sum(zeta_0[j,])) + sum(lgamma(zeta_0[j,]))
    for(m in 2:numMix){
      kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,m,])*(digamma(zeta[j,m]) - digamma(sum(zeta[j,])) ) 
      kl_theta <- kl_theta + sum(stateProb[,j,]*mixProb[,j,m,])*( digamma(gamma_shape[j,m-1]) - log(gamma_rate[j,m-1])) -
        sum(stateProb[,j,]*mixProb[,j,m,]*obs)*(gamma_shape[j,m-1]/gamma_rate[j,m-1] ) + h[j,m-1] #- h_j0[j,p-1]
    }
    kl_A <- kl_A + sum(rowSums(jtTransMat,dims = 2)[j,]*( digamma(alpha[j,]) - digamma(sum(alpha[j,])) )) +
      lgamma(sum(alpha[j,])) - sum(lgamma(alpha[j,]))   #-
    #lgamma(sum(alpha_0[j,])) + sum(lgamma(alpha_0[j,]))
  }
  kl_pi <- sum(apply(initProb, 2, sum)*(digamma(xi) - digamma(sum(xi)))) +
    lgamma(sum(xi)) - sum(lgamma(xi)) #- lgamma(sum(pi_0)) + sum(lgamma(pi_0)) 
  elbo <- - sum(log(ct)) - kl_theta - kl_A  - kl_C  - kl_pi
  return(elbo)
}

ELBO = function(numStates, numMix, numLoc, stateProb, mixProb, initProb, jtTransMat, ct, xi, alpha, zeta, gamma_shape, gamma_rate, obs, h){
  kl_C <- 0
  kl_A <- 0
  kl_theta <- 0
  kl_pi <- 0
  for(j in 1:numStates){
    for(l in numLoc){
      kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,1,,l])*(digamma(zeta[j,1,l]) - digamma(sum(zeta[j,,l])) ) + 
      lgamma(sum(zeta[j,,l])) - sum(lgamma(zeta[j,,l])) #-
      #lgamma(sum(zeta_0[j,])) + sum(lgamma(zeta_0[j,]))
      for(m in 2:numMix){
        kl_C <- kl_C + sum(stateProb[,j,]*mixProb[,j,m,,l])*(digamma(zeta[j,m,l]) - digamma(sum(zeta[j,,l])) ) 
        kl_theta <- kl_theta + sum(stateProb[,j,]*mixProb[,j,m,,l])*( digamma(gamma_shape[j,m-1,l]) - log(gamma_rate[j,m-1,l])) -
                                sum(stateProb[,j,]*mixProb[,j,m,,l]*obs[,,l])*(gamma_shape[j,m-1,l]/gamma_rate[j,m-1,l] ) + h[j,m-1,l] #- h_j0[j,p-1]
        }
    }
    kl_A <- kl_A + sum(rowSums(jtTransMat,dims = 2)[j,]*( digamma(alpha[j,]) - digamma(sum(alpha[j,])) )) +
            lgamma(sum(alpha[j,])) - sum(lgamma(alpha[j,]))   #-
            #lgamma(sum(alpha_0[j,])) + sum(lgamma(alpha_0[j,]))
  }
  kl_pi <-  sum(apply(initProb, 2, sum)*(digamma(xi) - digamma(sum(xi)))) +
            lgamma(sum(xi)) - sum(lgamma(xi)) #- lgamma(sum(pi_0)) + sum(lgamma(pi_0)) 
  elbo <- - sum(log(ct)) - kl_theta - kl_A  - kl_C  - kl_pi
  #print(paste(elbo,- sum(log(ct)),-kl_theta-kl_A-kl_C))
  return(elbo)
}
