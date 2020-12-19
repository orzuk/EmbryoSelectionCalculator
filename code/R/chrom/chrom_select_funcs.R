# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)
library(tensor)
library(olpsR) # for projection onto the simplex remotes::install_github("ngloe/olpsR")

chr.lengths <- c(0.0821,0.0799,0.0654,0.0628,0.0599,0.0564,0.0526,0.0479,0.0457,0.0441,
                                     0.0446,0.0440,0.0377,0.0353,0.0336,0.0298,0.0275,0.0265,0.0193,0.0213,0.0154,0.0168,0.0515)

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)



# Simulate a 3rd-order tensor of polygenic scores risk 
simulate_PS_chrom_disease_risk <- function(M, C, T, Sigma.T, Sigma.K, sigma.blocks, prev) # vectorize later    Z, K, E, n_sims = 1000) {
{
  X = array(0, dim=c(M, C, T))
  sim.vec <- 0
  if(sim.vec)
  {
    print("Sim vec")  
    for(i in 1:M)
      for(j in 1:C)
        X[i,j,] <- rmvnorm(1, mu = rep(0, T), sigma = Sigma.T) * sigma.blocks[i]
  } else
  {   # New: don't assume independence, use Kinship coefficient
    print("Sim mat")  
    for(i in 1:M)
      X[i,,] <- rmatnorm(1, Sigma.T, Sigma.K, M=matrix(0, nrow=T, ncol=C))  * sigma.blocks[i]  # V, M = matrix(0, nrow = nrow(U), ncol = nrow(V)))
  }   

# Next simulate disease D (not needed for now)
#  E <- rmvnorm(1, mu = rep(0, T), sigma = Sigma.T) * h2 # add envirounmental noise 
#  D = X + G_X + E 
  return(X)    
#  return(list(X=X, D=D))
}


# One-hot encoding as a matrix
c_vec_to_onehot <- function(c.vec, T)
{
  C.mat <- matrix(0, dim(c.vec)[1], T)
  for(i in 1:dim(c.vec)[1])
    C.mat[i, c.vec[i]] <- 1
  return(C.mat)
}


# One-hot encoding as a matrix
c_onehot_to_vec <- function(C.mat)
{
  return(max.col(C.mat))
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
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    loss <- sum(X.c^2 * loss.params$theta)    
  }
  
  if(loss.type == 'disease')  # weighted disease probability 
  {
    z.K <- qnorm(loss.params$K)
    loss <- sum(pnorm( (z.K-X.c)/sqrt(1-loss.params$h.ps)  ) * loss.params$theta) 
  }
  
  return(loss)
}  





# The gradient for the loss.
# Return a matrix of size: M*C
grad_loss_PS <- function(X, C, loss.type, loss.params)
{
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    return (2 * tensor_vector_prod(X, loss.params$theta * rep(1, M) %*% tensor_matrix_prod(X, C, 2)) )
  }
  
  if(loss.type == 'disease')
  {
    z.K <- qnorm(loss.params$K)
    Sigma.eps.inv <- 1/(1-sqrt(loss.params$h.ps))
    X.c <- compute_X_C_mat(X, C)
    return( tensor_vector_prod(X, loss.params$theta * Sigma.eps.inv * dnorm( (z.K-X.c)*Sigma.eps.inv )) ) 
  }
}  


# The hessian matrix of the loss. 
# Return a 4-order tensor of size: (M*C)*(M*C)
hessian_loss_PS <- function(X, C.mat, loss.type, loss.params)
{
  M = dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
  H <- array(0,c(M,C,M,C))  # 4-th order tensor  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    for(i in c(1:M))
      for(j in c(1:C))
      {
        H[i,j,,] <- 2*loss.params$theta[1] * X[i,j,1] * X[,,1]
        for(k in c(2:T))
          H[i,j,,] <- H[i,j,,] + 2*loss.params$theta[k] * X[i,j,k] * X[,,k]
      }
  }
  if(loss.type == 'disease')
  {
    z.K <- qnorm(loss.params$K)
    Sigma.eps.inv <- 1/(1-sqrt(loss.params$h.ps))
    X.c <- compute_X_C_mat(X, C.mat)
    alpha <- loss.params$theta * Sigma.eps.inv^3 * (z.K - X.c) * dnorm ((z.K - X.c)*Sigma.eps.inv)
    for(i in c(1:M))
      for(j in c(1:C))
      {
        H[i,j,,] <- alpha[1] * X[i,j,1] * X[,,1]
        for(k in c(2:T))
          H[i,j,,] <- H[i,j,,] + alpha[k] * X[i,j,k] * X[,,k]
      }
  }
  return(H)    
}


# Mutiply a tensor by matrix 
# X - a 3rd-prder tensor
# M - a matrix 
tensor_matrix_prod <- function(X, A, ax=3)
{
  if(ax == 2)
  {
    
    R = sweep(X[,1,], MARGIN=1, A[,1], `*`)
    for(i in 2:dim(X)[ax])
      R <- R + sweep(X[,i,], MARGIN=1, A[,i], `*`) # X[,i,] * A[,i]
  }
  if(ax == 3)
  {
    R = X[,,1] * A
    for(i in 2:dim(X)[ax])
      R <- R + X[,,i] * A
  }
  return(R)
}


# Mutiply a tensor by vector
# X - a 3rd-order tensor
# v - a vector 
tensor_vector_prod <- function(X, v, ax=3)
{
  if(ax == 3)
  {
    R = X[,,1] * v[1]
    for(i in 2:length(v))
      R <- R + X[,,i] * v[i]
  }
  return(R)
}


# Project each row of a matrix onto the simplex 
project_stochastic_matrix <- function(C.mat)
{
  C.proj <- C.mat
  for(i in 1:dim(C.mat)[1])
    C.proj[i,] <- projsplx_2(C.mat[i,])
  return (C.proj)
}
  



# Probability that a random Gaussian vector i.i.d. of dimension k is pareto out of n vectors 
pareto_P <- function(n, k)
{
  r <- 0
  for(i in c(1:n))
  {
    r <- r + choose(n-1,i-1) * (-1)^(i-1) / i^k
  }
  return(r)
}


# Probability that a random Gaussian vector i.i.d. of dimension k is pareto out of n vectors 
pareto_P2 <- function(n, k)
{
  if( k == 1)
    return(1/n)
  p <- 0
  for(i in c(1:n))
    p <- p + pareto_P2(i, k-1)
  return(p/n)    
}


# Compute recursivelt for all k and n to save time 
pareto_P_mat <- function(max.n, max.k)
{
  P.mat <- matrix(0, nrow=max.n, ncol=max.k)
  P.mat[,1] <- 1 / (1:max.n)  # k=1
  P.mat[1,] <- 1 # n=1
#  for(n in 1:max.n)  # fill for k=1
#    P.mat[n,1] = 1/n
  for(k in 2:max.k)
    for(n in 2:max.n)
      P.mat[n,k] <- ((n-1)*P.mat[n-1,k] + P.mat[n,k-1]) / n
  return(P.mat)
}


# Asymptotic approximation
pareto_P_approx <- function(n, k, order=2)
{
  if(order==1)
    return (log(n)^(k-1) / (n*factorial(k-1)))
  r <- 0
  for(i in 1:k)
    r <- r  + (log(n)^(i-1) * (-digamma(1))^(k-i) / (n*factorial(i-1)))
  return(r)
}

# Asymptotic approximation - take 2nd order 
#pareto_P_approx2 <- function(n, k)
#{
#  r <- 0
#  for(i in 1:k)
#    r <- r  + (log(n)^(i-1) * (-digamma(1))^(k-i) / (n*factorial(i-1)))
#  return(r)
#}


# Compute pareto optimal probability under indepndence with simmmulations 
pareto_P_sim <- function(n, k, iters=1000)
{
  n.pareto <- 0
  for(i in 1:iters)
    n.pareto <- n.pareto + length(get_pareto_optimal_vecs(matrix(runif(n*k), nrow=n, ncol=k))$pareto.inds) # Simulate vectors 
  return(n.pareto / (iters*n))
}
  



# Compare different methods for calculating p_k(n)
compare_pareto_P <- function(n.vec, k, C=2, iters = 1000)
{
  num.n <- length(n.vec)
  p.k <- rep(0, num.n)
  p.k.asymptotic <- rep(0, num.n, 1)
  p.k.asymptotic2 <- rep(0, num.n, 2)
  p.k.sim <- rep(0, num.n)
  p.k <- pareto_P_mat(max(n.vec), k)[n.vec, k]
  p.k.blocks <- rep(0, num.n)
  for(i in 1:num.n)
  {
    n <- n.vec[i]
#    p.k[i] <- pareto_P2(n, k)
    p.k.asymptotic[i] <- pareto_P_approx(n, k, 1)
    p.k.asymptotic2[i] <- pareto_P_approx(n, k, 2)
    if(iters > 0)
    {
#      p.k.sim[i] <- pareto_P_sim(n, k, iters)
      M <- round(log(n)/log(C))
      print(paste0("M=", M, " C=", C))
      p.k.blocks[i] <- pareto_P_block(C, M, k, iters)
    }
    # add also simulation  
  }
  return( list(p.k=p.k, p.k.sim=p.k.sim, p.k.blocks=p.k.blocks, 
               p.k.asymptotic=p.k.asymptotic, p.k.asymptotic2=p.k.asymptotic2) )
}
  
  
  
# Enumerate all multisubsets. Should work only for small n,k
#pareto_P3 <- function(n, k)
#{
#  return ( sum(1/colprods(combn(1:n,k)+0.0))/n ) # This doesn't include equalities !! 
#}



# Check if vector x is not dominated by any vector in X.mat
is_pareto_optimal <- function(x, X.mat)
{
  epsilon = 0.0000000000000000000000001
  if(isempty(X.mat))
    return(TRUE)
  if(is.null(dim(X.mat)))
    n.row <- 1
  else
    n.row <- dim(X.mat)[1]
  return( min(rowMaxs(t(replicate(n.row, x))+epsilon - X.mat, value=TRUE)) >= 0 )
}

# Extract only pareto-optimal vectors in a matrix
get_pareto_optimal_vecs <- function(X.mat)
{
  n = dim(X.mat)[1]
  is.pareto = rep(0,n)
  for(i in 1:n)
    is.pareto[i] = is_pareto_optimal(X.mat[i,], X.mat) #   max(colMins(t(replicate(n, X[i,]))+epsilon - X, value=TRUE)) >= 0
  pareto.inds <- which(as.logical(is.pareto))
  return(list(pareto.X=X.mat[pareto.inds,], pareto.inds=pareto.inds))
}




