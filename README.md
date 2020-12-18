# Variational Bayes Estimation of Hidden Markov Model for Daily Precipitation at a Single Location
Stochastic precipication generators which use HMMs for single site and multi-site rainfall are fairly common (Hughes and Guttorp, 1994; Robertson et al., 2006). However, they all use the EM-like Baum-Welch algorithm. Being a maximum likelihood method, it doesn't deal with model complexity well. I have been looking at adapting variational optimization to parameter estimation in HMMs for precipitation.  
Precipitation in our model is specified as a semi-continuous distribution with a point mass at zero and a mixture of two exponential distributions for postitive precipitation.  
Instead of the Baum-Welch algorithm, we use structured variational inference for parameter estimated. More details on the methodology can be found in this document.  
Things you can change in the code:
1. Input data, y
2. numYears and numDays, which correspond to the number of independent chains (temporally independent, not spatially) and length of individual chain. Data length N = numYears x numDays. For the simulation study, numYears = 1 and N = 1800. For Chesapeake, numYears = 20, and numDays = 92 (Jul1-Sep30).  
3. Number of states, K
4. Number of exponential mixture components, M. Note that m=0,1,..,M, where component 0 corresponds to dry days. In the code, m is indexed from 0, and so M includes the dry component as well.  
5. Maximum number of iterations to run, maxiter.  
6. The tolerance is used as a convergence criteria, and is based on the relative change in the Evidence Lower Bound (ELBO) between iterations. Tolerance in set to 10^-6.  
Priors:
1. The exponential distribution rate (lambda) has a Gamma prior. The prior'r shape for each component (gamma_jm) is ordered for identifiability. The first shape is <1 in each state (high rainfall) and the second shape is >1 (low rainfall). The rate parameter is kept constant at 2.  
2. Prior for the initial distribution (pi_1) is a symmetric Dirichlet with a concentration of 1. 
3. Prior for each row of the transition matrix (a_j) is a symmetric Dirichlet with a concentration of 10.
4. Prior for each row of mixture components (c_j) is Dirichlet. For the simulation study is asymmetric with a concentration of 10, for Chesapeake it is symmetric with a concentration of 12.  
