library(ConnMatTools)
#set.seed(Sys.time())
source('./functionDefns.R')
#data = read.table('cbaywet')
#data = read.table('data1')
N <- 1800
P <- 3
K <- 3
pi.pre = c(.7,.2,.1)
tmat.pre = matrix(c(.45,.35,.20,.30,.40,.30,.30,.30,.40),ncol = 3,byrow = T)
zeta.pre = matrix(c(0.3,.5,.2,0.3,.3,.4,0.5,.2,.3),ncol = 3, byrow = T)
lambda.pre = matrix(c(.08,1,.6,5,1,8),ncol = 2,byrow = T)

pi.post = rep(0,K)
tmat.post = matrix(0,nrow = K,ncol = K)
zeta.post = matrix(0,nrow = K, ncol = P)
lambda.post = matrix(0,nrow = K, ncol = P-1)

#### Initial estimates of hyperparameters
# gamma.mean = c(.1,2,.5,4,2,10)
# gamma.sd = c(.1,1,1,1,1,1)
# gammaparams = gammaParamsConvert(mean = gamma.mean, sd = gamma.sd)
# gamma.shape = gammaparams$shape
# gamma.rate = 1/gammaparams$scale
# gammaparams
# gamma_0 <- matrix(data = gamma.shape,nrow = K,ncol = P-1,byrow = T) #shape of exponential rate
# delta_0 <- matrix(data = gamma.rate,nrow = K,ncol = P-1,byrow = T) #rate of exponential rate
gamma_0 <- matrix(data = c(.5,2,1.5,9,2,16),nrow = K,ncol = P-1,byrow = T) #shape of exponential rate
delta_0 <- matrix(data = c(2,2,2,2,2,2),nrow = K,ncol = P-1,byrow = T) #rate of exponential rate
zeta_0 <- matrix(data = c(6,8,6,6,7,7,8,6,6),nrow = K,ncol = P,byrow = T)/(P-1) #Dirichlet prior parameters for mixing probabilities
alpha_0 <- matrix(rep(10,K^2),nrow = K,byrow = T)/K #Dirichlet prior for transition matrix rows
pi_0 = rep(1,K)/K #Dirichlet prior for initial probabilities
h_j0 = gamma_0*log(delta_0) - lgamma(gamma_0) #constant terms in prior (log)
convergence = matrix(NA,300,20)
original = matrix(NA,N,1000)
simulated = matrix(NA,N,1000)

  for(sim in 1:1000){
 set.seed(sim)
#y <- data$V2#[1:180]
y=simulate.hmm.uni(duration=N,numChain = 1,numStates = K,numMix = P,initDist = pi.pre,tmat = tmat.pre,
                   zeta = zeta.pre,lambda = lambda.pre)
original[,sim] = y
# N <- length(y) #number of days
numDays <- N
numYears <- N/numDays
# sum(y==0)
y2 = matrix(y,ncol = numYears,byrow = F)

del_y0 <- ifelse(y2==0,1,0)
# K <- 3 #number of states
# P <- 3 #number of mixture components (first component would be dry)

#### Variational Bayes iterations
maxiter <- 300 #number of iterations to run the code for
dic <- rep(0,maxiter)
dic_old <- 50000 
dic[1] <- 25000
elbo <- rep(0,maxiter)
elbo_old <- -50000 
elbo[1] <- -25000
tol <- 10^(-6)
iter <- 1
# improvement <- (dic_old-dic[1])#/dic_old
improvement <- (elbo_old-elbo[1])/elbo_old


# initial estimates of latent variables

a_jl <- matrix(NA,nrow = K, ncol = K) #posterior state kernel
b_tj <- array(dim = c(numDays,K,numYears)) #posterior emission kernel
p_j = matrix(nrow = numYears,ncol = K) #posterior initial probability kernel


ct_wide = matrix(nrow = numDays,ncol = numYears)
fvar_wide = array(dim = c(numDays,K,numYears))
bvar_wide = array(dim = c(numDays,K,numYears))


s_t <- sample(1:K,N,replace = T) #random assignment of states
# s_t = v.sim
s_t = matrix(s_t,nrow = numDays)

q_tj <- array(dim = c(numDays,K,numYears)) #posterior probability of stationary distribution
q_tjp <- array(0,dim = c(numDays,K,P,numYears))  #posterior probability of mixture assignments (conditional)
jt_transition_mat <- array(dim=c(K,K,numDays-1,numYears)) #posterior joint transition probability matrix
q_p <- matrix(nrow = numYears, ncol = K) #posterior initial probability
# y.median = median(y)
# y.components = list()
# y.components[[1]] = which(y>y.median & y>0)
# y.components[[2]] = which(y<y.median & y>0)
for(m in 1:numYears)
for(j in 1:K){
  q_tj[,j,m] <- sum(s_t[,m]==j)/numDays
  q_tjp[,j,1,] <- del_y0
  for(p in 2:P)
    #q_tjp[y.components[[p-1]],j,p,m] <- 1
    q_tjp[,j,p,m] <- (1-del_y0)/(P-1)
  for(t in 1:numDays){
  sum_q = sum(q_tjp[t,j,,m])
  q_tjp[t,j,-1,m] = q_tjp[t,j,-1,m]/sum_q
  }
}
q_p <- q_tj[1,,]
b_star = array(dim = c(numDays,K,P))


gamma_jp <- gamma_0 #posterior shape of exponential rate
delta_jp <- delta_0 #posterior rate of exponential rate
zeta_j <- zeta_0 #posterior Dirichlet parameters for mixing probabilities
alpha_j <- alpha_0 #posterior Dirichlet parameters for transition matrix rows
h_jp <- matrix(nrow = K,ncol = P-1) #constant terms in posterior (log)
pi_j = pi_0 #posterior Dirichlet paramters for initial distribution

# # Variables to run and store VB results
# gamma_post <- matrix(nrow = maxiter,ncol = K*(P-1))
# delta_post <- matrix(nrow = maxiter,ncol = K*(P-1))
# zeta_post <- matrix(nrow = maxiter,ncol = K*P)
# alpha_post <- matrix(nrow = maxiter,ncol = K*K)
# logh_post <- matrix(nrow = maxiter,ncol = K)



  while((abs(improvement) > tol | improvement <0) & iter<maxiter)
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

  #Forward and Backward recursions
  for(m in 1:numYears)
    {
    #fvar_out = forward_recursion(numObs = numDays, numStates = K, initDist = q_tj[1,,m], a = a_jl, b = b_tj[,,m])
    fvar_out = forward_recursion(numObs = numDays, numStates = K, initDist = p_j, a = a_jl, b = b_tj[,,m])
    fvar_wide[,,m] = fvar_out[[1]]
    ct_wide[,m] = fvar_out[[2]]
  
    bvar_wide[,,m] = backward_recursion(numObs = numDays, numStates = K, a = a_jl, b = b_tj[,,m], ct=ct_wide[,m])
    }
  q_p = t(q_tj[1,,])
  
  # Update the posterior state probabilities

  for(m in 1:numYears)
    {
    q_tj[,,m] <- fvar_wide[,,m]*bvar_wide[,,m]/apply(bvar_wide[,,m]*fvar_wide[,,m],1,sum)
    for(j in 1:K)
      for(l in 1:K){
        jt_transition_mat[j,l,,m] <- fvar_wide[-numDays,j,m]*a_jl[j,l]*b_tj[-1,l,m]*bvar_wide[-1,l,m]
      }
    #jt_transition_mat[,,,m] <- jt_transition_mat[,,,m]/rep(colSums(jt_transition_mat[,,,m],dims = 2),each=K^2)
    for(t in 1:(numDays-1)){
      denom = 0
      for(j in 1:K)
        for(l in 1:K){
    denom = denom+jt_transition_mat[j,l,t,m]
        }
      jt_transition_mat[,,t,m] = jt_transition_mat[,,t,m]/denom
    }
  }
  
  ####Update the mixing probabilties
  for(m in 1:numYears)
    for(j in 1:K)
      {
      b_star[,j,1] <- del_y0[,m]
      for(p in 2:P)
        b_star[,j,p] <- (1-del_y0[,m])*exp(digamma(zeta_j[j,p]) - digamma(sum(zeta_j[j,])) + 
                            digamma(gamma_jp[j,p-1]) - log(delta_jp[j,p-1]) - y2[,m]*(gamma_jp[j,p-1]/delta_jp[j,p-1]) )
        #b_star[,j,p] <- (1-del_y0[,m])*(zeta_j[j,p]/sum(zeta_j[j,]))*(gamma_jp[j,p-1]/delta_jp[j,p-1])*exp(-y2[,m]*(gamma_jp[j,p-1]/delta_jp[j,p-1]) )
        sum_b <- rowSums(b_star[,j,],dims = 1)
        q_tjp[,j,-1,m] <- b_star[,j,-1]/sum_b
        q_tjp[,j,1,m] = b_star[,j,1]
      }


  ####Update hyperparameters
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
  
    # # Store paramter updates from this iteration
    # gamma_post[iter,] <- gamma_jp
    # delta_post[iter,] <- delta_jp
    # zeta_post[iter,] <- as.vector(zeta_j)
    # alpha_post[iter,] <- as.vector(alpha_j)
    # logh_post[iter,] = h
  
    # DIC and ELBO calculations
    # pd0 <- 0
    # pd1 <- 0
    # pd2 <- 0
    kl_C <- 0
    kl_A <- 0
    kl_theta <- 0
    kl_pi <- 0
    for(j in 1:K)
      {
      # pd0 <- pd0 + sum(q_tj[,j,]*q_tjp[,j,1,])*(log(zeta_j[j,1]) - log(sum(zeta_j[j,])) - 
      #             digamma(zeta_j[j,1]) + digamma(sum(zeta_j[j,])) )
      kl_C = kl_C + sum(q_tj[,j,]*q_tjp[,j,1,])*(digamma(zeta_j[j,1]) - digamma(sum(zeta_j[j,])) ) + 
                    lgamma(sum(zeta_j[j,])) - sum(lgamma(zeta_j[j,])) -
                    lgamma(sum(zeta_0[j,])) + sum(lgamma(zeta_0[j,]))
      for(p in 2:P)
        {
        # pd2 <- pd2 + sum(q_tj[,j,]*q_tjp[,j,p,])*( log(gamma_jp[j,p-1]) - digamma(gamma_jp[j,p-1]) +
        #             log(zeta_j[j,p]) - log(sum(zeta_j[j,])) - digamma(zeta_j[j,p]) + digamma(sum(zeta_j[j,])) )
        kl_C = kl_C + sum(q_tj[,j,]*q_tjp[,j,p,])*(digamma(zeta_j[j,p]) - digamma(sum(zeta_j[j,])) ) 
        kl_theta = kl_theta + sum(q_tj[,j,]*q_tjp[,j,p,])*( digamma(gamma_jp[j,p-1]) - log(delta_jp[j,p-1])) -
                              sum(q_tj[,j,]*q_tjp[,j,p,]*y2)*(gamma_jp[j,p-1]/delta_jp[j,p-1] ) + h_jp[j,p-1] - h_j0[j,p-1]
      }
      #for(l in 1:K)
      # pd1 <- pd1 + rowSums(jt_transition_mat,dims = 2)[j,l]*( log(alpha_j[j,l]) - log(sum(alpha_j[j,])) - digamma(alpha_j[j,l]) + digamma(sum(alpha_j[j,])) )
        kl_A <- kl_A + sum(rowSums(jt_transition_mat,dims = 2)[j,]*( digamma(alpha_j[j,]) - digamma(sum(alpha_j[j,])) )) +
                    lgamma(sum(alpha_j[j,])) - sum(lgamma(alpha_j[j,]))   -
                    lgamma(sum(alpha_0[j,])) + sum(lgamma(alpha_0[j,]))
      }
  
    # pd3 = sum(q_p*(log(pi_j) - log(sum(pi)) - digamma(pi_j) + digamma(sum(pi_j))))
    kl_pi = sum(apply(q_p, 2, sum)*(digamma(pi_j) - digamma(sum(pi_j)))) +
            lgamma(sum(pi_j)) - sum(lgamma(pi_j)) - lgamma(sum(pi_0)) + sum(lgamma(pi_0)) 
    # dic[iter+1] <- 4*(pd0+sum(pd1)+pd2+pd3)  + 2*sum(log(ct_wide))
    elbo[iter+1] = -sum(log(ct_wide)) - kl_theta - kl_A  - kl_C  - kl_pi
    elbo_old <- elbo[iter]
    improvement <- (elbo[iter+1]-elbo_old)/abs(elbo_old)

    # dic_old <- dic[iter]
    # improvement <- (dic_old-dic[iter+1])#/dic_old
    #if(iter%%10==0)
    #print(paste(iter, dic[iter+1],improvement))
    #print(paste(iter,elbo[iter],improvement))
    #print(paste(iter,sum(log(ct_wide)),kl_C,kl_A,kl_theta,kl_pi,improvement))  
    #print(paste(gamma_jp[2,1],delta_jp[2,1],improvement))
    iter <- iter+1
    }

plot <- F
if(plot==T){
  for(p in 2:P)
    plot(1:maxiter,gamma_post[,p-1],'l')
  for(p in 2:P)
    plot(1:maxiter,delta_post[,p-1],'l')
  for(p in 1:P)
    plot(1:maxiter,zeta_post[,p],'l') 
}

params = post_param(numStates = K,numMix = P,gamma.post = gamma_jp,delta.post = delta_jp, zeta = zeta_j, alpha = alpha_j, pi_j = pi_j)
params
# pi.pre
# tmat.pre
# zeta.pre
# lambda.pre
# pi.post = params$InitDist
# tmat.post = params$TransMat
# zeta.post = params$MixProb
# lambda.post = params$RainRate
pi.post = pi.post + params$InitDist
tmat.post = tmat.post + params$TransMat
zeta.post = zeta.post + params$MixProb
lambda.post = lambda.post + params$RainRate
# if(sim%%10==0){
#   print(paste('Results at Simulation',sim))
#   print(pi.post/sim)
#   print(tmat.post/sim)
#   print(zeta.post/sim)
#   print(lambda.post/sim)
# }
 print(paste('Completed simulation',sim,'at iteration',iter))
# convergence[,sim] = elbo
 y.sim=simulate.hmm.uni(duration=N,numChain = 1,numStates = K,numMix = P,initDist = params$InitDist,tmat = params$TransMat,
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
write.table(original,'original',col.names = F,row.names = F)
write.table(simulated,'simulated',col.names = F,row.names = F)
del_y0.og = ifelse(original==0,0,1)
del_y0.sim=ifelse(simulated[,1]==0,0,1)
a=apply(del_y0.og,2,function(x)mean(x==0))
b=apply(del_y0.sim,2,function(x)mean(x==0))
plot(a,b)

a=apply(original,2,function(x)mean(x[x>0]))
b=apply(simulated,2,function(x)mean(x[x>0]))
plot(a,b)


y.rle = rle(ifelse(y==0,0,1))
y.sim.rle = rle(del_y0.sim)
y.rle.0 = y.rle$lengths[y.rle$values==0]
y.sim.rle.0 = y.sim.rle$lengths[y.sim.rle$values==0]
summary(y.rle.0)
summary(y.sim.rle.0)
boxplot(y.rle.0,y.sim.rle.0)
