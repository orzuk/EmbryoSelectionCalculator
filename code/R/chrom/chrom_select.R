# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)

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
      X[i,j,] <- rmvnorm(1, mu = rep(0, T), sigma = Sigma.T) * sigma.blocks[i]
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




# The loss 
loss_PS <- function(X.c, loss.type, loss.params)
{
#  X.c <- compute_X_C_mat(X, C)
  if(loss.type == "quant")
    loss <- sum(X.c * loss.params$theta)
  
  if(loss.type == "balancing")
  {
    loss <- sum(X.c^2 * loss.params$theta)    
  }
  
  if(loss.type == 'disease')  # weighted disease probability 
  {
    z.K <- qnorm(loss.params$K)*sqrt(loss.params$h.ps)
    loss <- sum(dnorm( (X.c - z.K)/sqrt(loss.params$h.ps)  ) * loss.params$theta) 
  }
  
  return(loss)
}  

hessian_loss_PS <- function()
{
  
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


# Check if vector x is not dominated by any vector in X.mat
is_pareto_optimal <- function(x, X.mat)
{
  epsilon = 0.0000000000000000000000001
  if(isempty(X.mat))
    return(TRUE)
  return( max(colMins(t(replicate(dim(X.mat)[1], x))+epsilon - X.mat, value=TRUE)) >= 0 )
}

# Extract only pareto-optimal vectors in a matrix
get_pareto_optimal_vecs <- function(X)
{
  epsilon = 0.0000000000000000000000001
  n = dim(X)[1]
  print(n)
  is.pareto = rep(0,n)
  for(i in 1:n)
    is.pareto[i] = max(colMins(t(replicate(n, X[i,]))+epsilon - X, value=TRUE)) >= 0
#  which(is.pareto)
  return(X[which(as.logical(is.pareto)),])
}
  
# A branch and bound algorithm 
optimize_C_branch_and_bound <- function(X, loss.C, loss.params)
{
  M = dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
  print("Start B&B")
  cur.X <- X[1,,]
  cur.X <- get_pareto_optimal_vecs(cur.X) # Save only Pareto-optimal vectors 
  L <- dim(cur.X)[1]
  
  for( i in c(2:M))
  {  
    print(paste0("Start B&B i=", i))

    new.X <- c()
    for(j in c(1:L))  # loop over all vectors in the current stack      
      for  (c in c(1:C))  # loop over possible vectors to add 
      {
#        print("Start if")
        if(is_pareto_optimal(X[i,c,], new.X))
          new.X  <- rbind(new.X, X[i,c,])
      }
    cur.X <- new.X
    L <- dim(new.X)[1]  
    print(paste("Stack Size:", L))
  }
    
   
  # Finally find the cost-minimizer out of  the Pareto-optimal vectors
  L <- dim(cur.X)[1]
  loss.vec <- rep(0, L)
  for(i in 1:L)
    loss.vec[i] <- loss_PS(cur.X[i,], loss.C, loss.params)
  print(loss.vec)
  i.min <- which.min(loss.vec) # loss_PS(cur.X, loss.C, loss.params))
    
  return(cur.X[i.min,])
}


