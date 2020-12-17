library(ConnMatTools)
#set.seed(Sys.time())
source('./functionDefns.R')
#' Set read.data=T if you are running it for the first time and need to read the data in
#' The Chesapeake Data takes a while to read. Don't want to re-read it at every iteration
read.data <- T 

#' For checking ELBO convergence,
#' Seed is set to 0 then 20 seeds are sampled
#' [1] 836 679 129 930 509 471 299 270 978 187 307 597 277 874 950
#' [16] 494 330 775 841 591
#' Model is run for each of those seeds
#' ELBO for each run saved in ELBO.txt
#' Seed #8 (270) had highest ELBO
#' Random grid point under seed 0 is 1422

set.seed(0)
loc   <- sample(1:1927,1) #1422
seeds <- sample(1:1000,20) 

if(read.data==T){
  data  <- read.table('data/cbaywet')
  dim(data)
  y     <- data[,loc+1]
}

N         <- 1840  # Data length
M         <- 3          # Number of mixture components
K         <- 3          # Number of states
numDays   <- 92       # Number of days
numYears  <- N/numDays  # Number of years
y2        <- matrix(y, ncol = numYears, byrow = F) # Split y vector into a matrix where each year's data is a column
del_y0    <- ifelse(y2==0, 1, 0) # variable to store dry-wet days

set.seed(seeds[8]) 

#### Convergence criteria
maxiter   <- 300 # number of iterations to run the code for
dic       <- rep(0,maxiter)
dic_old   <- 50000 
dic[1]    <- 25000
elbo      <- rep(0,maxiter)
elbo_old  <- -50000 
elbo[1]   <- -25000
tol       <- 10^(-6)
iter      <- 1
improvement_dic   <- (dic_old-dic[1])/dic_old
improvement_elbo  <- (elbo_old-elbo[1])/elbo_old
numPrint  <- 10 # Print ELBO/DIC after a certain number of iterations

#### Initial estimates of hyperparameters
g0          <- c(runif(K,0,1), runif(K,1,20)) # Gamma parameters for the states
gamma_0     <- matrix(g0, nrow = K, ncol = M-1, byrow = F) # shape of exponential rate
delta_0     <- matrix(data = rep(2,K*(M-1)), nrow = K, ncol = M-1, byrow = T) # rate of exponential rate
zeta_0      <- matrix(data = c(4,4,4,4,4,4,4,4,4), nrow = K, ncol = M, byrow = T)/(M-1) # Dirichlet prior parameters for mixing probabilities
alpha_0     <- matrix(rep(10,K^2), nrow = K, byrow = T)/K # Dirichlet prior for transition matrix rows
xi_0        <- rep(1,K)/K # Dirichlet prior for initial probabilities
h_j0        <- gamma_0*log(delta_0) - lgamma(gamma_0) # constant terms in prior (log)

#### Initial estimates of latent variables
a_jk        <- matrix(NA, nrow = K, ncol = K) # posterior state kernel
b_tj        <- array(dim = c(numDays, K, numYears)) # posterior emission kernel
a_1j         <- matrix(nrow = numYears, ncol = K) # posterior initial probability kernel

ct     <- matrix(nrow = numDays, ncol = numYears)
fvar   <- array(dim = c(numDays, K, numYears))
bvar   <- array(dim = c(numDays, K, numYears))
b_star      <- array(dim = c(numDays, K, M))

q_tj        <- array(dim = c(numDays, K, numYears)) # posterior probability of stationary distribution
q_tjm       <- array(0, dim = c(numDays, K, M, numYears))  # posterior probability of mixture assignments (conditional)
q_jk <- array(dim=c(K, K, numDays-1, numYears)) # posterior joint transition probability matrix
q_1j         <- matrix(nrow = numYears, ncol = K) # posterior initial probability

#### Posterior hyperparameter variables
gamma_jm    <- gamma_0 # posterior shape of exponential rate
delta_jm    <- delta_0 # posterior rate of exponential rate
zeta_j      <- zeta_0 # posterior Dirichlet parameters for mixing probabilities
alpha_j     <- alpha_0 # posterior Dirichlet parameters for transition matrix rows
h_jm        <- matrix(nrow = K, ncol = M-1) # constant terms in posterior (log)
xi_j        <- xi_0 # posterior Dirichlet parameters for initial distribution

#### Variables to run and store VB results
gamma_iters  <- matrix(nrow = maxiter, ncol = K*(M-1))
delta_iters  <- matrix(nrow = maxiter, ncol = K*(M-1))
zeta_iters   <- matrix(nrow = maxiter, ncol = K*M)
alpha_iters  <- matrix(nrow = maxiter, ncol = K*K)
logh_iters   <- matrix(nrow = maxiter, ncol = K)

#### Variational optimization begins
while((abs(improvement_elbo) > tol | improvement_elbo < 0) & iter < maxiter){
  
  #### Update latent variables
  a_1j <- exp(digamma(xi_j) - digamma(sum(xi_j)))
  for(j in 1:K){
  a_jk[j,] <- exp(digamma(alpha_j[j,]) - digamma(sum(alpha_j[j,])))
    for(n in 1:numYears){
      tmp <- 0
      for(m in 2:M)
      tmp <- tmp + (1 - del_y0[, n])*exp(digamma(zeta_j[j,m]) - digamma(sum(zeta_j[j,])) + 
                    digamma(gamma_jm[j,m-1]) - log(delta_jm[j,m-1]) - y2[,n]*(gamma_jm[j,m-1]/delta_jm[j,m-1] ) )
      b_tj[,j,n] <- del_y0[,n]*exp( digamma(zeta_j[j,1]) - digamma(sum(zeta_j[j,])) ) + tmp 
      }
    }

  #### Forward and Backward recursions
  for(n in 1:numYears){
    fvar_out <- forward_recursion(numObs = numDays, numStates = K, initDist = a_1j, a = a_jk, b = b_tj[,,n])
    fvar[,,n] <- fvar_out[[1]]
    ct[,n] <- fvar_out[[2]]
  
    bvar[,,n] <- backward_recursion(numObs = numDays, numStates = K, a = a_jk, b = b_tj[,,n], ct=ct[,n])
    }

  #### Update the posterior state probabilities
  for(n in 1:numYears){
    q_tj[,,n] <- fvar[,,n]*bvar[,,n]/apply(bvar[,,n]*fvar[,,n], 1, sum)
    for(j in 1:K)
      for(l in 1:K){
        q_jk[j,l,,n] <- fvar[-numDays,j,n]*a_jk[j,l]*b_tj[-1,l,n]*bvar[-1,l,n]
      }
    #q_jk[,,,m] <- q_jk[,,,m]/rep(colSums(q_jk[,,,m],dims = 2),each=K^2) # Need to fix this expression
    for(t in 1:(numDays-1)){
      denom <- 0
      for(j in 1:K)
        for(k in 1:K){
          denom <- denom+q_jk[j,k,t,n]
        }
      q_jk[,,t,n] <- q_jk[,,t,n]/denom
    }
  }
  q_1j = t(q_tj[1,,])
  
  #### Update the mixing probabilties
  for(n in 1:numYears)
    for(j in 1:K){
      b_star[,j,1] <- del_y0[,n]
      for(m in 2:M)
        b_star[,j,m] <- (1-del_y0[,n])*exp(digamma(zeta_j[j,m]) - digamma(sum(zeta_j[j,])) + 
                        digamma(gamma_jm[j,m-1]) - log(delta_jm[j,m-1]) - y2[,n]*(gamma_jm[j,m-1]/delta_jm[j,m-1]) )
        sum_b <- rowSums(b_star[,j,], dims = 1)
        q_tjm[,j,-1,n] <- b_star[,j,-1]/sum_b
        q_tjm[,j,1,n] <- b_star[,j,1]
      }

  #### Update hyperparameters
  for(j in 1:K){
    xi_j[j] <- xi_0[j] + sum(q_1j[,j])
    zeta_j[j,1] <- zeta_0[j,1] + sum(q_tj[,j,]*q_tjm[,j,1,])
    for(m in 2:M){
      zeta_j[j,m] <- zeta_0[j,m] + sum(q_tj[,j,]*q_tjm[,j,m,])
      gamma_jm[j,m-1] <- gamma_0[j,m-1] + sum(q_tj[,j,]*q_tjm[,j,m,])
      delta_jm[j,m-1] <- delta_0[j,m-1] + sum(q_tj[,j,]*q_tjm[,j,m,]*y2)
      }
    alpha_j[j,] <- alpha_0[j,] + rowSums(q_jk,dims = 2)[j,]
    }
    h_jm <- gamma_jm*log(delta_jm) - lgamma(gamma_jm)
  
  #### Store paramter updates from this iteration
  gamma_iters[iter,]  <- gamma_jm
  delta_iters[iter,]  <- delta_jm
  zeta_iters[iter,]   <- as.vector(zeta_j)
  alpha_iters[iter,]  <- as.vector(alpha_j)
  
  #### DIC and ELBO calculations
  dic[iter+1]     <- uDIC(numStates = K, numMix = M, stateProb = q_tj, mixProb = q_tjm, initProb = q_1j, jtTransMat = q_jk,
                          ct = ct, xi = xi_j, alpha = alpha_j, zeta = zeta_j, gamma_shape = gamma_jm, gamma_rate = delta_jm)[[2]]
  dic_old         <- dic[iter]
  improvement_dic <- (dic_old-dic[iter+1])/dic_old

  elbo[iter+1]    <- uELBO(numStates = K, numMix = M, stateProb = q_tj, mixProb = q_tjm, initProb = q_1j, jtTransMat = q_jk,
                           ct = ct, xi = xi_j, alpha = alpha_j, zeta = zeta_j, gamma_shape = gamma_jm, gamma_rate = delta_jm,
                           obs = y2, h = h_jm)
  elbo_old        <- elbo[iter]
  improvement_elbo <- (elbo[iter+1]-elbo_old)/abs(elbo_old)


  if(iter%%numPrint==0)
  #print(paste(iter, dic[iter+1],improvement_dic))
  print(paste(iter, elbo[iter], improvement_elbo))
  iter <- iter+1
}

params      <- post_param(numStates = K, numMix = M, gamma.post = gamma_jm, delta.post = delta_jm, 
                          zeta = zeta_j, alpha = alpha_j, xi = xi_j)
params
pi.post     <- params$InitDist
tmat.post   <- params$TransMat
zeta.post   <- params$MixProb
lambda.post <- params$RainRate

y.sim <- simulate.hmm.uni(numDays = 92, numYears = 20, numStates = 3, numMix = 3, initDist = pi.post, tmat = tmat.post,
                         zeta = zeta.post, lambda = lambda.post)
sum(y.sim==0)
sum(y==0)
# write.table(y,'cbay1422',row.names = F,col.names = F)
# write.table(y.sim,'cbaysim',row.names = F,col.names = F)

#' Variables you can plot to visually check performance
#' plot(2:iter,elbo[2:iter],'l')
#' plot(2:iter,dic[2:iter],'l') # Deviance Information criterion approximation as per McGrory & Titterington (2009)
#' gamma_iters, delta_iters, zeta_iters, gamma_iters store each iteration's values of the hyperparamters
#' If you plot y and y.sim, plot for individual years since the data is disjoint
