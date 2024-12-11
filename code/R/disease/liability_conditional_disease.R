library(mvtnorm)
# Sample from a Gaussian distribution conditional on disease status
# Input: 
# D - vector of disease status (binary)
# K - vector of prevalences
# Sigma - matrix of genetic correlations for the scores. (We assume Var(Y_t)=1, E[Y_t]=0)
# Sigma.eps - matrix of correlations for the non-score part (genetic + environment). 
# iters - num. of importance sampling iterations
# sampling.method - rejection (default) - rejection sampling (can be very slow for large T), 
#                   importance: move the Gaussians to have as means Z_K, the thresholf for the liability
# Output: 
# x - matrix of samples from the conditional density f_X| D where X ~ N(0, Sigma), 
#     and D|X given by the liability threshold model.
# w - vector of weights for the samples, that is w_i is the weights given to x_i 
#     These samples can represent an the empirical conditional distribution f-hat-cond(x) = \sum_i w_i 1_{x=x_i}
liability_conditional_disease <- function(D, K, Sigma, Sigma.eps, iters, sampling.method="rejection") {
  T <- length(D)  # number of traits 
  # Try naive sampling (rejection sampling)
  ctr <- 0
  x <- matrix(rep(0, iters*T), nrow=iters)
  w <- rep(1, iters) # Default: constant weights
  if(sampling.method == "importance")
  {
    mu <- D*qnorm(K/2) + (1-D)*qnorm(max((1+K)/2, 1-K/2))  # Better to get the conditional expectation !!! 
#    Sigma.eps.samp <- 
  }
   else
  {
    mu <- rep(0, T)
#    Sigma.eps.samp 
  }  
  
  while(ctr < iters)
  {
    scores <- rmvnorm(iters, mean=mu, sigma=Sigma)
    eps <- rmvnorm(iters, sigma=Sigma.eps)
    D.mat <- scores + eps < rep(K, iters, 1) 
    good.ind <- which(apply(D.mat == matrix(rep(D, iters), nrow=iters, byrow = TRUE), 1, all))
    if(length(good.ind)>0)
    {
      good.ind <- good.ind[1:min(length(good.ind),iters-ctr)]
      x[(ctr+1):(ctr+length(good.ind)),] <- scores[good.ind,]
      if(sampling.method == "importance") # compute the likelihood ratio of the two Gaussians
      {
        w[(ctr+1):(ctr+length(good.ind))] <- 
          dmvnorm(scores[good.ind,], mean=rep(0, T), sigma=Sigma) / 
          dmvnorm(scores[good.ind,], mean=mu, sigma=Sigma)   # QU: Maybe we need also to take epsilon in the likelihood? 
#        w[(ctr+1):(ctr+length(good.ind))] <- 
#          dmvnorm(scores[good.ind,], mean=rep(0, T), sigma=Sigma+Sigma.eps) / 
#          dmvnorm(scores[good.ind,], mean=mu, sigma=Sigma+Sigma.eps)   # QU: Maybe we need also to take epsilon in the likelihood? This gives more stable results but wrong!!!
      }
      ctr <- ctr + length(good.ind)
      print(c("Ctr: ", ctr))
    }
  }
  
  return(list(x=x,w=w/sum(w)))  # w always normalized to sum to one
}
  
  
# Usage example:  (easy case)
Sigma <- matrix(c(0.8, 0.5, 0.5, 0.7), nrow=2)
Sigma.eps <- matrix(c(0.2, 0.1, 0.1, 0.3), nrow=2)
K <- c(0.01, 0.1)
iters <- 1000
D <- c(1,1) # sick in both diseases
epdf <- liability_conditional_disease(D, K, Sigma, Sigma.eps, iters)

print(colMeans(epdf$x))  # unweighted version 

# Hard case (reject too many)
T <- 2
Sigma <- matrix(rep(0.5/T, T*T), nrow=T)
diag(Sigma) <- 0.5
Sigma.eps <- matrix(rep(0.25/T, T*T), nrow=T)
diag(Sigma.eps) <- 0.5
K <- rep(0.1, T) # c(0.01, 0.1)
iters <- 10000
D <- c(rep(1, T/2), rep(0, T/2)) # c(1,1) # sick in both diseases
epdf <- liability_conditional_disease(D, K, Sigma, Sigma.eps, iters)

print(colMeans(epdf$x))  # un-weighted version 
print(colSums(epdf$x * matrix(rep(epdf$w, T), nrow=iters)))  # weighted version (should be equal here)

# Very unstable !!! need a different, better importance distribution 
epdf.importance <- liability_conditional_disease(D, K, Sigma, Sigma.eps, iters, sampling.method = "importance")
#print(colMeans(epdf.importance$x))  # un-weighted version 
print(colSums(epdf.importance$x * matrix(rep(epdf.importance$w, T), nrow=iters)))  # weighted version (should be equal here)
