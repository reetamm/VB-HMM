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

markov.generator <- function(x,alpha,PI) #look into alternative nomenclature
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
######################################
# This function provides a concordance matrix for an MMC
similarity.matrix <- function(states)
{
  dims <- dim(states)
  n <- dims[1]
  m <- dims[2]
  cormat <- diag(m)
  for(i in 1:(m - 1))
  {
    for(j in (i + 1):m)
    {
      cormat[i,j] <- sum(states[,i] == states[,j])/n
      cormat[j,i] <- cormat[i,j]
    }
  }
  return(cormat)
}

################################
# Function to compute Cramer's V for an MMC
cramer.corr <- function(states)
{
  dims <- dim(states)
  n <- dims[1]
  m <- dims[2]
  cormat <- diag(m)
  for(i in 1:(m-1))
  {
    for(j in (i+1):m)
    {
      cormat[i,j] <- cramerV(states[,i], states[,j])
      cormat[j,i] <- cormat[i,j]
    }
  }
  return(cormat)
}
#######################
kruskal.gamma <- function(states)
{
  require(DescTools)
  dims <- dim(states)
  n <- dims[1]
  m <- dims[2]
  cormat <- diag(m)
  for(i in 1:(m - 1))
  {
    for(j in (i + 1):m)
    {
      cormat[i,j] <- GoodmanKruskalGamma(states[,i], states[,j])
      cormat[j,i] <- cormat[i,j]
    }
  }
  return(cormat)
}


#############
# Simulate a multivariate Markov chain from a uniform sample when parameters are known
generate.mmc <- function(uniform.sample, alpha, PI, head = T)
{
  # alpha <- M x J and PI <- J x J x M
  #if(sum(dim(t(alpha))==dim(PI[1,,]))!=2)
    #print('Dimensions do not matchin alpha and PI')
 # if(ncol(uniform.sample)!= nrow(alpha))
    #print('Uniform sample of incorrect dimensions')
  n <- nrow(uniform.sample)
  M <- ncol(uniform.sample)
  states <- matrix(NA,nrow = n, ncol = M)
  for(i in 1:M)
  {
    states[,i] <- markov.generator(uniform.sample[,i], alpha[i,], PI[,,i])
    states[,i] <- (as.character(states[,i]))
  }
  if(head == T)
    print(head(states, 6))
  return(states)
}


#################################
createsigma <- function(rho, sd = c(1,1,1,1), M = 4)
{
  require(lme4)
    R <- sd*diag(M)
  for(i in 1:M)
    for(j in 1:M)
      if(i != j)
        R[i,j] <- rho 
  # Convert correlation matrix for var-cov matrix in Normal
  sigma <- sdcor2cov(R)
  output <- list(R = R, sigma = sigma)
  return(output)
}

##############################
getparams <- function(dataset)
{
  require(markovchain)
  M <- dim(dataset)[2]
  J <- length(table(dataset))
  alpha <- matrix(NA, ncol = J, nrow = M, byrow = T)
  PI <- array(NA, dim = c(J, J, M))
  
  for(i in 1:M)
  {
    temp.pi <- markovchainFit(dataset[,i])
    PI[,,i] <- temp.pi$estimate[1:J]
    alpha[i,] <- table(dataset[,i])/1656
  }
  
  output<-list(alpha=alpha,PI=PI)
  return(output)
}

#################################
creategamma <- function(numChains = 4, numStates = 3, shape = c(1,2,5), rate = c(1,1,1))
{
  gamma <- array(NA, dim = c(4,3,2))
  for(i in 1:4)
  {
    gamma[i,,1] <- shape
    gamma[i,,2] <- rate
  }
  return(gamma)
}

#####################################
# Alternate, less feature rich version of createsigma() above.
create.cor.matrix <- function(M,rho)
{
  R <- diag(M)
  for(i in 1:M)
    for(j in 1:M)
      if(i != j)
        R[i,j] <- rho 

  return(R)
}

makeMatPD <- function(corrMatrix) {
  eigen <- 0.00001
  cholStatus <- try(outChol <- chol(corrMatrix), silent=TRUE) 
  cholError  <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  
  # fix the correlation matrix if chol fails
  while(cholError) {
    
    # replace -ve eigen values with a small +ve number
    neweig  <- eigen(corrMatrix)
    neweig$values[neweig$values < 0] <- eigen
    
    # inv = transp for eig vectors
    corrMatrix <- neweig$vectors %*% diag(neweig$values) %*% t(neweig$vectors) 
    
    # try chol again
    cholStatus <- try(outChol <- chol(corrMatrix), silent=TRUE) 
    cholError  <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)
  }
  
  return (outChol)
}


extractmonth <- function(dataset, months)
{
  outdata <- dataset[1,]
  for(month in months)
  {
    loc <- which(month(dataset$date) == month)
    tmp <- dataset[loc,]
    outdata <- rbind(outdata,tmp)
  }
  outdata <- outdata[-1,]
  return(outdata)
}

RMSE <- function(m, o){
  sqrt(mean((m - o)^2))
}


###################Forward and backward recursion functions

forward_recursion <- function(numObs = N, numStates = K, initDist = c(1,rep(0,K-1)), a = NULL, b = NULL){
  ct <- rep(1,numObs)
  fvar = matrix(nrow = numObs,ncol = numStates)
  fvar_star = matrix(nrow = numObs,ncol = numStates)
  fvar_tilde = matrix(nrow = numObs,ncol = numStates)
  fvar[1,] <- initDist*b[1,]
  ct[1] <- 1/sum(fvar[1,])
  fvar_tilde[1,] <- fvar[1,]*ct[1]
  
  for(t in 2:numObs){
    for(j in 1:numStates){
      fvar[t,j] = sum(fvar[t-1,]*a[,j])*b[t,j]
      fvar_star[t,j] <- sum(fvar_tilde[t-1,]*a[,j])*b[t,j]
      #print(t)
    }
    ct[t] <- 1/sum(fvar_star[t,])
    fvar_tilde[t,] <- fvar_star[t,]*ct[t]
  }
  output = list(fvar_tilde, ct)
  return(output)
}

backward_recursion_2 <- function(numObs = N, numStates = K, a = NULL, b = NULL, ct = NULL){
  bvar_star = matrix(nrow = numObs,ncol = numStates)
  bvar_tilde = matrix(nrow = numObs,ncol = numStates)
  bvar_tilde[numObs,] <- 1/numStates
  for(t in (numObs-1):1){
    for(j in 1:numStates){
      bvar_star[t,j] <- sum(a[j,]*bvar_tilde[t+1,]*b[t+1,])
    }
    bvar_tilde[t,] <- bvar_star[t,]/sum(bvar_star[t,])
  }
  return(bvar_tilde)
}

backward_recursion <- function(numObs = N, numStates = K, a = NULL, b = NULL, ct = NULL){
  bvar_star = matrix(nrow = numObs,ncol = K)
  bvar_tilde = matrix(nrow = numObs,ncol = K)
  #bvar_tilde[N,] <- 1/K
  bvar_tilde[numObs,] <- ct[numObs]
  for(t in (numObs-1):1){
    for(j in 1:numStates){
      bvar_star[t,j] <- sum(a[j,]*bvar_tilde[t+1,]*b[t+1,])
    }
    #bvar_tilde[t,] <- bvar_star[t,]/sum(bvar_star[t,])
    bvar_tilde[t,] <- bvar_star[t,]*ct[t]
  }
  return(bvar_tilde)
}

post_param = function(numStates = K, numMix = P, gamma.post = gamma_jp, delta.post = delta_jp, zeta = zeta_j, alpha = alpha_j, pi_j = pi_j)
{
  tmat.post <- matrix(nrow = numStates,ncol = numStates)
  lambda.post <- gamma.post/delta.post
  zeta.post <- matrix(nrow = numStates,ncol = numMix)
      for(j in 1:numStates){
      zeta.post[j,] <- zeta[j,]/sum(zeta[j,])
      tmat.post[j,] <- alpha[j,]/sum(alpha[j,])
      }
  pi.post = pi_j/sum(pi_j)
  
  output = list('InitDist' = pi.post, 'TransMat' = tmat.post, 'MixProb' = zeta.post, 'RainRate' = lambda.post)
  return(output)
}

simulate.hmm.uni = function(duration, numChain=1,numStates, numMix, initDist,tmat,zeta,lambda){
  v <- c()
  sim.data <- c()
  v.sim <- c()
  u <- runif(duration)
  v.sim <- markov.generator(u,initDist,tmat)
  mix.unif <- runif(duration)
  y.sim <- c()
  for(t in 1:duration)
  {
    which.mixture <- initial.state(mix.unif[t],zeta[v.sim[t],])
    if(which.mixture==1)
      y.sim[t] <- 0 else
        y.sim[t] <- rexp(1,lambda[v.sim[t],which.mixture-1])
  }
  return(y.sim)
}
