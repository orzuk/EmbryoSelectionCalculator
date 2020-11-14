# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)


chr.lengths <- c(0.0821,0.0799,0.0654,0.0628,0.0599,0.0564,0.0526,0.0479,0.0457,0.0441,
                                     0.0446,0.0440,0.0377,0.0353,0.0336,0.0298,0.0275,0.0265,0.0193,0.0213,0.0154,0.0168,0.0515)

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)



# Simulate a tensor of polygenic scores risk 
simulate_PS_chrom_disease_risk <- function(M, C, T, Sigma.T, sigma.blocks, prev) # vectorize later    Z, K, E, n_sims = 1000) {
{
  X = array(0, dim=c(M, C, T))
  
  for(i in 1:M)
    for(j in 1:C)
    {
      X[i,j,] <- rmvnorm(1, mean = rep(0, T), sigma = Sigma.T) * sigma.blocks[i]
    }
  
  # Next simulate disease D
#  E <- rmvnorm(1, mean = rep(0, T), sigma = Sigma.T) * h2 # add envirounmental noise 
#  D = X + G_X + E 
  return(X)    
#  return(list(X=X, D=D))
}



# Reduce a tensor of size M*C*T and a matrix of size M*C to a risk vector of length T
compute_X_c_vec <- function(X, c.vec)
{
  T <- dim(X)[3]
  M = dim(X)[1]
  X.c <- rep(0, T)
  for( i in 1:M)
    X.c = X.c + X[i,c.vec[i],]
  
  return(X.c)
}

# Reduce a tensor of size M*C*T and a matrix of size M*C to a risk vector of length T
compute_X_C_mat <- function(X, C.mat)
{
  T <- dim(X)[3]
  M = dim(X)[1]
  X.c <- rep(0, T)
  for( i in 1:M)
    X.c = X.c + colSums(sweep(X[i,,], MARGIN=1, C.mat[i,], `*`))  
  
  return(X.c)
}



# The gradient for the loss 
grad_loss_PS <- function(X, C, loss.type, loss.params)
{
  if(loss.type == "balancing")
  {
    
  }
  
  if(loss.type == 'disease')
  {
    
  }
}  

hessian_loss_PS <- function()
{
  
}



# a function for optimizing the selection of chromosomes:
# Run the optimization to find the optimal C:
optimize_C <- function(X, C_init, loss_C, loss.params)
{
  
}


# Probability that a random Gaussian vector i.i.d. of dimension k is pareto out of n vectors 
pareto_P <- function(n, k)
{
  r <- 0
  for(i in c(0:n))
  {
    r <- r + choose(n,i) * (-1)^i / (1+i)^k
  }
  return(r)
}
