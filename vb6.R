library(ConnMatTools)
#set.seed(Sys.time())
source('./functionDefns.R')
source('./datasim.R')
N         <- length(y) # number of days
P         <- 3
K         <- 3
numDays   <- 1800
numYears  <- N/numDays
y2        <- matrix(y,ncol = numYears,byrow = F)
del_y0    <- ifelse(y2==0,1,0)

#### Initial estimates of hyperparameters
gamma_0     <- matrix(data = c(.5,2,1.5,9,2,16),nrow = K,ncol = P-1,byrow = T) # shape of exponential rate
delta_0     <- matrix(data = c(2,2,2,2,2,2),nrow = K,ncol = P-1,byrow = T) # rate of exponential rate
zeta_0      <- matrix(data = c(6,8,6,6,7,7,8,6,6),nrow = K,ncol = P,byrow = T)/(P-1) # Dirichlet prior parameters for mixing probabilities
alpha_0     <- matrix(rep(10,K^2),nrow = K,byrow = T)/K # Dirichlet prior for transition matrix rows
pi_0        <- rep(1,K)/K # Dirichlet prior for initial probabilities
h_j0        <- gamma_0*log(delta_0) - lgamma(gamma_0) # constant terms in prior (log)

#### Variational Bayes iterations
maxiter   <- 300 # number of iterations to run the code for
dic       <- rep(0,maxiter)
dic_old   <- 50000 
dic[1]    <- 25000
elbo      <- rep(0,maxiter)
elbo_old  <- -50000 
elbo[1]   <- -25000
tol       <- 10^(-7)
iter      <- 1
improvement_dic   <- (dic_old-dic[1])/dic_old
improvement_elbo  <- (elbo_old-elbo[1])/elbo_old

#### Initial estimates of latent variables
a_jl        <- matrix(NA,nrow = K, ncol = K) # posterior state kernel
b_tj        <- array(dim = c(numDays,K,numYears)) # posterior emission kernel
p_j         <- matrix(nrow = numYears,ncol = K) # posterior initial probability kernel

ct_wide     <- matrix(nrow = numDays,ncol = numYears)
fvar_wide   <- array(dim = c(numDays,K,numYears))
bvar_wide   <- array(dim = c(numDays,K,numYears))
b_star <- array(dim = c(numDays,K,P))

s_t         <- sample(1:K,N,replace = T) # random assignment of states
s_t         <- matrix(s_t,nrow = numDays)

q_tj        <- array(dim = c(numDays,K,numYears)) # posterior probability of stationary distribution
q_tjp       <- array(0,dim = c(numDays,K,P,numYears))  # posterior probability of mixture assignments (conditional)
jt_transition_mat <- array(dim=c(K,K,numDays-1,numYears)) # posterior joint transition probability matrix
q_p         <- matrix(nrow = numYears, ncol = K) # posterior initial probability

#### Initial assignment of observations to state and mixture component
for(m in 1:numYears)
for(j in 1:K){
  q_tj[,j,m] <- sum(s_t[,m]==j)/numDays
  q_tjp[,j,1,] <- del_y0
  for(p in 2:P)
    q_tjp[,j,p,m] <- (1-del_y0)/(P-1)
  for(t in 1:numDays){
  sum_q <- sum(q_tjp[t,j,,m])
  q_tjp[t,j,-1,m] <- q_tjp[t,j,-1,m]/sum_q
  }
}
q_p <- q_tj[1,,]

#### Posterior hyperparameter variables
gamma_jp    <- gamma_0 # posterior shape of exponential rate
delta_jp    <- delta_0 # posterior rate of exponential rate
zeta_j      <- zeta_0 # posterior Dirichlet parameters for mixing probabilities
alpha_j     <- alpha_0 # posterior Dirichlet parameters for transition matrix rows
h_jp        <- matrix(nrow = K,ncol = P-1) # constant terms in posterior (log)
pi_j        <- pi_0 # posterior Dirichlet paramters for initial distribution

#### Variables to run and store VB results
gamma_post  <- matrix(nrow = maxiter,ncol = K*(P-1))
delta_post  <- matrix(nrow = maxiter,ncol = K*(P-1))
zeta_post   <- matrix(nrow = maxiter,ncol = K*P)
alpha_post  <- matrix(nrow = maxiter,ncol = K*K)
logh_post   <- matrix(nrow = maxiter,ncol = K)

#### Variational optimization begins
while((abs(improvement_elbo) > tol | improvement_elbo <0) & iter<maxiter)
  #while(iter<maxiter)
    {
  
    #### Update latent variables
    p_j <- exp(digamma(pi_j) - digamma(sum(pi_j)))
    for(j in 1:K)
      {
      a_jl[j,] <- exp(digamma(alpha_j[j,]) - digamma(sum(alpha_j[j,])))
      for(m in 1:numYears)
        {
        tmp <- 0
        for(p in 2:P)
        tmp <- tmp + (1- del_y0[,m])*exp( digamma(zeta_j[j,p]) - digamma(sum(zeta_j[j,])) + digamma(gamma_jp[j,p-1]) - log(delta_jp[j,p-1]) - y2[,m]*(gamma_jp[j,p-1]/delta_jp[j,p-1] ) )
        b_tj[,j,m] <- del_y0[,m]*exp( digamma(zeta_j[j,1]) - digamma(sum(zeta_j[j,])) ) + tmp 
        }
      }

  #### Forward and Backward recursions
  for(m in 1:numYears)
    {
    #fvar_out = forward_recursion(numObs = numDays, numStates = K, initDist = q_tj[1,,m], a = a_jl, b = b_tj[,,m])
    fvar_out <- forward_recursion(numObs = numDays, numStates = K, initDist = p_j, a = a_jl, b = b_tj[,,m])
    fvar_wide[,,m] <- fvar_out[[1]]
    ct_wide[,m] <- fvar_out[[2]]
  
    bvar_wide[,,m] <- backward_recursion(numObs = numDays, numStates = K, a = a_jl, b = b_tj[,,m], ct=ct_wide[,m])
    }
  q_p = t(q_tj[1,,])
  
  #### Update the posterior state probabilities
  for(m in 1:numYears)
    {
    q_tj[,,m] <- fvar_wide[,,m]*bvar_wide[,,m]/apply(bvar_wide[,,m]*fvar_wide[,,m],1,sum)
    for(j in 1:K)
      for(l in 1:K){
        jt_transition_mat[j,l,,m] <- fvar_wide[-numDays,j,m]*a_jl[j,l]*b_tj[-1,l,m]*bvar_wide[-1,l,m]
      }
    #jt_transition_mat[,,,m] <- jt_transition_mat[,,,m]/rep(colSums(jt_transition_mat[,,,m],dims = 2),each=K^2)
    for(t in 1:(numDays-1)){
      denom <- 0
      for(j in 1:K)
        for(l in 1:K){
    denom <- denom+jt_transition_mat[j,l,t,m]
        }
      jt_transition_mat[,,t,m] <- jt_transition_mat[,,t,m]/denom
    }
  }
  
  #### Update the mixing probabilties
  for(m in 1:numYears)
    for(j in 1:K)
      {
      b_star[,j,1] <- del_y0[,m]
      for(p in 2:P)
        b_star[,j,p] <- (1-del_y0[,m])*exp(digamma(zeta_j[j,p]) - digamma(sum(zeta_j[j,])) + 
                            digamma(gamma_jp[j,p-1]) - log(delta_jp[j,p-1]) - y2[,m]*(gamma_jp[j,p-1]/delta_jp[j,p-1]) )
        sum_b <- rowSums(b_star[,j,],dims = 1)
        q_tjp[,j,-1,m] <- b_star[,j,-1]/sum_b
        q_tjp[,j,1,m] <- b_star[,j,1]
      }


  #### Update hyperparameters
  for(j in 1:K)
    {
    pi_j[j] <- pi_0[j] + mean(q_p[,j])
    zeta_j[j,1] <- zeta_0[j,1] + sum(q_tj[,j,]*q_tjp[,j,1,])
    for(p in 2:P)
      {
      zeta_j[j,p] <- zeta_0[j,p] + sum(q_tj[,j,]*q_tjp[,j,p,])
      gamma_jp[j,p-1] <- gamma_0[j,p-1] + sum(q_tj[,j,]*q_tjp[,j,p,])
      delta_jp[j,p-1] <- delta_0[j,p-1] + sum(q_tj[,j,]*q_tjp[,j,p,]*y2)
      }
    alpha_j[j,] <- alpha_0[j,] + rowSums(jt_transition_mat,dims = 2)[j,]
    }
    h_jp <- gamma_jp*log(delta_jp) - lgamma(gamma_jp)
  
    #### Store paramter updates from this iteration
    gamma_post[iter,] <- gamma_jp
    delta_post[iter,] <- delta_jp
    zeta_post[iter,] <- as.vector(zeta_j)
    alpha_post[iter,] <- as.vector(alpha_j)

    #### DIC and ELBO calculations
    dic[iter+1]     <- uDIC(numStates=K,numMix=P,stateProb=q_tj,mixProb=q_tjp,A_t=jt_transition_mat,
                       ct=ct_wide,pi_1=pi_j,alpha=alpha_j,zeta=zeta_j,shape=gamma_jp,rate=delta_jp,initProb=q_p)
    dic_old         <- dic[iter]
    improvement_dic <- (dic_old-dic[iter+1])/dic_old

    elbo[iter+1]    <- uELBO(numStates=K,numMix=P,stateProb=q_tj,mixProb=q_tjp,A_t=jt_transition_mat,
                         ct=ct_wide,pi_1=pi_j,alpha=alpha_j,zeta=zeta_j,shape=gamma_jp,rate=delta_jp,obs=y2,h=h_jp,initProb=q_p)
    elbo_old        <- elbo[iter]
    improvement_elbo <- (elbo[iter+1]-elbo_old)/abs(elbo_old)


    if(iter%%10==0)
    #print(paste(iter, dic[iter+1],improvement_dic))
    print(paste(iter,elbo[iter],improvement_elbo))
    iter <- iter+1
    }

params <- post_param(numStates = K,numMix = P,gamma.post = gamma_jp,delta.post = delta_jp, zeta = zeta_j, alpha = alpha_j, pi_j = pi_j)
params
pi.post <- params$InitDist
tmat.post <- params$TransMat
zeta.post <- params$MixProb
lambda.post <- params$RainRate

plot <- T
if(plot==T){
  for(p in 2:P)
    plot(1:maxiter,gamma_post[,p-1],'l')
  for(p in 2:P)
    plot(1:maxiter,delta_post[,p-1],'l')
  for(p in 1:P)
    plot(1:maxiter,zeta_post[,p],'l') 
}
# y.sim=simulate.hmm.uni(duration=N,numChain = 1,numStates = K,numMix = P,initDist = params$InitDist,tmat = params$TransMat,
#                         zeta = params$MixProb,lambda = params$RainRate)