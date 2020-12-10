pi.pre      <- c(.7,.2,.1)
tmat.pre    <- matrix(c(.45,.35,.20,.30,.40,.30,.30,.30,.40),ncol = 3,byrow = T)
zeta.pre    <- matrix(c(0.3,.5,.2,0.3,.3,.4,0.5,.2,.3),ncol = 3, byrow = T)
lambda.pre  <- matrix(c(.08,1,.6,5,1,8),ncol = 2,byrow = T)

y           <- simulate.hmm.uni(numDays,numYears,numStates = K,numMix = P,initDist = pi.pre,tmat = tmat.pre,
                              zeta = zeta.pre,lambda = lambda.pre)

