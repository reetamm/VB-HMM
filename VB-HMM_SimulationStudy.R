library(ConnMatTools)
#set.seed(Sys.time())
source('./functionDefns.R')
numPrint=10
N <- 1800
M <- 3
K <- 3
source('./datasim.R')

pi.post = rep(0,K)
tmat.post = matrix(0,nrow = K,ncol = K)
zeta.post = matrix(0,nrow = K, ncol = M)
lambda.post = matrix(0,nrow = K, ncol = M-1)

#### Initial estimates of hyperparameters
gamma_0 <- matrix(data = c(.5,2,1.5,9,2,16),nrow = K,ncol = M-1,byrow = T) #shape of exponential rate
delta_0 <- matrix(data = c(2,2,2,2,2,2),nrow = K,ncol = M-1,byrow = T) #rate of exponential rate
zeta_0 <- matrix(data = c(6,8,6,6,7,7,8,6,6),nrow = K,ncol = M,byrow = T)/(M-1) #Dirichlet prior parameters for mixing probabilities
alpha_0 <- matrix(rep(10,K^2),nrow = K,byrow = T)/K #Dirichlet prior for transition matrix rows
xi_0 = rep(1,K)/K #Dirichlet prior for initial probabilities
h_j0 = gamma_0*log(delta_0) - lgamma(gamma_0) #constant terms in prior (log)
convergence = matrix(NA,300,20)
original = matrix(NA,N,1000)
simulated = matrix(NA,N,1000)
sim=1
for(sim in 1:1000){
 set.seed(sim)
y=simulate.hmm.uni(numDays=N,numYears = 1,numStates = K,numMix = M,initDist = pi.pre,tmat = tmat.pre,
                   zeta = zeta.pre,lambda = lambda.pre)
original[,sim] = y
numDays <- N
numYears <- N/numDays
y2 = matrix(y,ncol = numYears,byrow = F)
del_y0 <- ifelse(y2==0,1,0)

#### Variational Bayes iterations
maxiter <- 300 #number of iterations to run the code for
elbo <- rep(0,maxiter)
elbo_old <- -50000 
elbo[1] <- -25000
tol <- 10^(-6)
iter <- 1
improvement_elbo <- (elbo_old-elbo[1])/elbo_old


# initial estimates of latent variables
a_jk <- matrix(NA,nrow = K, ncol = K) #posterior state kernel
b_tj <- array(dim = c(numDays,K,numYears)) #posterior emission kernel
a_1j = matrix(nrow = numYears,ncol = K) #posterior initial probability kernel


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

  while((abs(improvement_elbo) > tol | improvement_elbo <0) & iter<maxiter)
    {
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

  #Forward and Backward recursions
  for(m in 1:numYears)
    {
  for(n in 1:numYears){
    fvar_out <- forward_recursion(numObs = numDays, numStates = K, initDist = a_1j, a = a_jk, b = b_tj[,,n])
    fvar[,,n] <- fvar_out[[1]]
    ct[,n] <- fvar_out[[2]]
  
    bvar[,,n] <- backward_recursion(numObs = numDays, numStates = K, a = a_jk, b = b_tj[,,n], ct=ct[,n])
  }
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
  
  ####Update the mixing probabilties
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

  ####Update hyperparameters
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

   #### DIC and ELBO calculations

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

pi.post = pi.post + params$InitDist
tmat.post = tmat.post + params$TransMat
zeta.post = zeta.post + params$MixProb
lambda.post = lambda.post + params$RainRate
if(sim%%10==0){
  print(paste('Results at Simulation',sim))
  print(pi.post/sim)
  print(tmat.post/sim)
  print(zeta.post/sim)
  print(lambda.post/sim)
}
 print(paste('Completed simulation',sim,'at iteration',iter))
 convergence[,sim] = elbo
 
 y.sim = simulate.hmm.uni(numDays=N,numYears = 1,numStates = K,numMix = M,initDist = params$InitDist,tmat = params$TransMat,
                  zeta = params$MixProb,lambda = params$RainRate)
 simulated[,sim] = y.sim
  }
# convergence2 = convergence[2:200,]
# convergence_df = data.frame(convergence2)
# library(ggplot2)
# library(tidyr)
# convergence_long = gather(convergence_df,seed, dic,X1:X20,factor_key = T)
# head(convergence_long)
# convergence_long$iter = rep(2:200,20)
# 
# ggplot(convergence_long,aes(y=dic,col=seed,x=iter)) + geom_line() + theme(legend.position = 'none') +
#   labs(x='Iteration count',y='DIC value')
# 
# write.table(original,'original',col.names = F,row.names = F)
# write.table(simulated,'simulated',col.names = F,row.names = F)
# del_y0.og = ifelse(original==0,0,1)
# del_y0.sim=ifelse(simulated[,1]==0,0,1)
# a=apply(del_y0.og,2,function(x)mean(x==0))
# b=apply(del_y0.sim,2,function(x)mean(x==0))
# plot(a,b)
# 
# a=apply(original,2,function(x)mean(x[x>0]))
# b=apply(simulated,2,function(x)mean(x[x>0]))
# plot(a,b)
# 
# 
# y.rle = rle(ifelse(y==0,0,1))
# y.sim.rle = rle(del_y0.sim)
# y.rle.0 = y.rle$lengths[y.rle$values==0]
# y.sim.rle.0 = y.sim.rle$lengths[y.sim.rle$values==0]
# summary(y.rle.0)
# summary(y.sim.rle.0)
# boxplot(y.rle.0,y.sim.rle.0)
