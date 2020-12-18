# Variational Bayes Estimation of a Hidden Markov Model for Daily Precipitation at a Single Location
Stochastic precipication generators which use HMMs for single site and multi-site rainfall are fairly common ([Hughes and Guttorp (1994)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/93WR02983); [Robertson et al. (2006)](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1256/qj.05.75)). However, they all use the EM-like Baum-Welch algorithm. Being a maximum likelihood method, it doesn't deal with model complexity well. I have been looking at adapting variational optimization to parameter estimation in HMMs for precipitation. Variational inference for HMMs are covered in [Ji et al. (2006)](https://ieeexplore.ieee.org/document/1597110) and [McGrory and Titterington (2009)](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1467-842X.2009.00543.x). A tutorial for the Forward-Backward algorithm for HMMs, which is also used in the variational framework, can be found in [Rabiner (1989)](https://web.ece.ucsb.edu/Faculty/Rabiner/ece259/Reprints/tutorial%20on%20hmm%20and%20applications.pdf).
Precipitation in our model is specified as a semi-continuous distribution with a point mass at zero and a mixture of two exponential distributions for postitive precipitation.  
Instead of the Baum-Welch algorithm, we use structured variational inference for parameter estimated. More details on the methodology can be found in [this document](https://github.com/reetamm/VB-HMM/blob/main/vb-hmm.pdf).  

## List of variables in the code and their interpertation
Variable | Description
---------|------------
`y` | Input data.
`numYears`, `numDays` | the number of independent chains (temporally, not spatially) and length of individual chain, respectively .
`N` | Data length, equal to `numYears` x `numDays`. [1]
`K` | Number of states (indices j,k = 1,..,K).
`M` | Number of exponential mixture components (Indices m=0,1,..,M, where component 0 corresponds to dry days). [2]
`maxiter` | Maximum number of iterations to run.
`tolerance` | convergence criteria based on the relative change in the Evidence Lower Bound (ELBO) between iterations. Tolerance in set to `10^-6`.  
`gamma_jm` | Shape parameter of Gamma prior for exponential rainfall. [3]
`delta_jm` | Rate parameter of Gamma prior for exponential rainfall; kept constant at 2.
`xi_j` | Parameters of a symmetric Dirichlet prior for the initial distribution (`pi_1`); concentration of 1.
`alpha_j` | Parameters of a symmetric Dirichlet prior for each row of the transition matrix (`a_j`); concentration of 10.
`zeta_j` | Parameters of a Dirichlet prior for each row of mixture components (`c_j`).[4] 

[1] For the simulation study, `numYears` = 1 and `N` = 1800. For Chesapeake, `numYears` = 20, and `numDays` = 92 (Jul 1-Sep 30).  
[2] In the code, `m` is indexed from 0, and so `M` includes the dry component as well.  
[3] Ordered for identifiability. The first shape is <1 in each state (high rainfall) and the second shape is >1 (low rainfall)
[4] For the simulation study is asymmetric with a concentration of 10, for Chesapeake it is symmetric with a concentration of 12.  

## Work in progress:  
1. Viterbi algorithm for this model  
2. Deviance Information Criterion (DIC) calculations as a model selection criteria  
3. Extension to multi-site precipitation  
4. Converting this into an R package
