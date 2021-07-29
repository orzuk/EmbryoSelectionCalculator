# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)
library(tensor)
library(olpsR) # for projection onto the simplex remotes::install_github("ngloe/olpsR")
library(BH)
library(ecr)
library(Rfast)

Rcpp::sourceCpp("cpp/chrom_funcs.cpp")  # fast functions  


chr.lengths <- c(0.0821,0.0799,0.0654,0.0628,0.0599,0.0564,0.0526,0.0479,0.0457,0.0441,
                                     0.0446,0.0440,0.0377,0.0353,0.0336,0.0298,0.0275,0.0265,0.0193,0.0213,0.0154,0.0168,0.0515)

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)



# Simulate a 3rd-order tensor of polygenic scores risk 
# Inputs: 
# M - number of embryos
# Number of total chromosomes (blocks)
# T - number of traits
# Sigma.T - covariance of traits
# Sigma.K - covariance of individuals
# sigma.blocks - ??
# prev - vector of disease prevalences 
# 
# Output: 
# X - a 3rd order tensor with scores. X[i,j,k] is score of embryo i in block j for trait k
simulate_PS_chrom_disease_risk <- function(M, C, T, Sigma.T, Sigma.K, sigma.blocks, prev) # vectorize later    Z, K, E, n_sims = 1000) {
{
  X = array(0, dim=c(M, C, T))
  sim.vec <- 0
  if(sim.vec)
  {
#    print("Sim vec")  
    for(i in 1:M)
      for(j in 1:C)
        X[i,j,] <- rmvnorm(1, mu = rep(0, T), sigma = Sigma.T) * sigma.blocks[i]
  } else
  {   # New: don't assume independence, use Kinship coefficient
#    print("Sim mat")  
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


# Reduce a tensor of size M*C*T and a vector of length M to a risk vector of length T by summation 
# Input: 
# X - 3rd order tensor
# c.vec - indices
# Output: 
# X.c - a vector of length T that is sum of M vectors , one per each slice 
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
  # What if C.mat is of dimension 1?
  M = dim(X)[1]
  C = dim(X)[2]
  T <- dim(X)[3]
  X.c <- rep(0, T)
  if(is.null(dim(C.mat)[1]) || (C==1))  # C=1 , C.mat is a vector 
  {
    for(i in 1:M)
    {
#      print(paste0("i=", i))
#      print(X[i,,] * C.mat[i])
      X.c = X.c + X[i,,] * C.mat[i]  # colSums(sweep(X[i,,], MARGIN=1, C.mat[i], `*`))  
    }
  }
  else
    for(i in 1:M)
      X.c = X.c + colSums(sweep(X[i,,], MARGIN=1, C.mat[i,], `*`))  
  
  return(X.c)
}



# The loss 
# Input: 
# X.c - a vector of length C
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# loss - a scalar 
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



# The loss for vectors 
# Input: 
# X.c.mat - a matrix with multiple vectors of length C
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# loss.vec - a vector of losses for each row of X.c.mat  
loss_PS_mat <- function(X.c.mat, loss.type, loss.params)
{
  #  X.c <- compute_X_C_mat(X, C)
  if(loss.type == "quant")
    loss.vec <- X.c.mat %*% loss.params$theta # scalar product 
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    loss.vec <- (X.c.mat*X.c.mat) %*% loss.params$theta    
  }
  
  if(loss.type == 'disease')  # weighted disease probability 
  {
    z.K <- qnorm(loss.params$K)
    loss.vec <- t(pnorm( (z.K-t(X.c.mat))/sqrt(1-loss.params$h.ps)  )) %*% loss.params$theta 
  }
  return(loss.vec)
}  





# The gradient for the loss.
# Return a matrix of size: M*C
# Can also add a regularizer 
grad_loss_PS <- function(X, C, loss.type, loss.params)
{
  M = dim(X)[1]

  if(!("eta" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$eta <- 0 
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    return (2 * tensor_vector_prod(X, loss.params$theta * rep(1, M) %*% tensor_matrix_prod(X, C, 2))  - loss.params$eta)
  }
  
  if(loss.type == 'disease')
  {
    z.K <- qnorm(loss.params$K)
    Sigma.eps.inv <- 1/(1-sqrt(loss.params$h.ps))
    X.c <- compute_X_C_mat(X, C)
    return( tensor_vector_prod(X, loss.params$theta * Sigma.eps.inv * dnorm( (z.K-X.c)*Sigma.eps.inv )) - loss.params$eta) # added regularizer  
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
# X - a 3rd-prder tensor of dims: M*C*T
# A - a matrix. Dimension depends on axis to multiply by: 
# ax=3: A is a matrix of size M*C
# ax=2: A is a matrix of size M*T
# ax=1: A is a matrix of size C*T
tensor_matrix_prod <- function(X, A, ax=3)
{
  if(ax == 2)
  {
    if(dim(X)[ax] == 1)
      R = sweep(X[,1,], MARGIN=1, A, `*`) # only one value
    else
    {
      R = sweep(X[,1,], MARGIN=1, A[,1], `*`)
      for(i in 2:dim(X)[ax])
        R <- R + sweep(X[,i,], MARGIN=1, A[,i], `*`) # X[,i,] * A[,i]
    }
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
  log.binom <- lchoose(n-1, 0:(n-1)) - k * log(1:n)
  max.exp <- max(log.binom)
#  log.binom <- log.binom - max.exp
  return (   sum(exp(log.binom - max.exp) * (-1)^c(0:(n-1))) * exp(max.exp) )
  
  # naive implementation  
#  r <- 0
#  for(i in c(1:n))
#  {
#    r <- r + choose(n-1,i-1) * (-1)^(i-1) / i^k
#  }
#  return(r)
}


# Probability that a random Gaussian vector i.i.d. of dimension k is pareto out of n vectors 
# Recursive, slower. Not used.
pareto_P2 <- function(n, k)
{
  if( k == 1)
    return(1/n)
  p <- 0
  for(i in c(1:n))
    p <- p + pareto_P2(i, k-1)
  return(p/n)    
}


# Probability that a random Gaussian vector i.i.d. of dimension k is pareto out of n vectors. 
# Compute recursively for all k and n to up to max.n, max.k respectively to save time 
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


# Asymptotic approximation for p_n,k (first or other orders )
pareto_P_approx <- function(n, k, order=2)
{
  if(order==1)
    return (log(n)^(k-1) / (n*factorial(k-1)))
  r <- 0
  for(i in 1:k)
    r <- r  + (log(n)^(i-1) * (-digamma(1))^(k-i) / (n*factorial(i-1)))
  return(r)
}



# Compute pareto optimal probability under indepndence with simulations 
pareto_P_sim <- function(n, k, iters=1000)
{
  n.pareto <- 0
  for(i in 1:iters)
    n.pareto <- n.pareto + length(get_pareto_optimal_vecs(matrix(runif(n*k), nrow=n, ncol=k))$pareto.inds) # Simulate vectors 
  return(n.pareto / (iters*n))
}
  



# Compare different methods for calculating p_k(n)
# Input: 
# n.vec - values of n
# k - value of k
# C - ???
# iters - how many simulations to run 
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
  cpp.flag <- FALSE
  if(cpp.flag)
  {
    par <- get_pareto_optimal_vecs_rcpp(X.mat)
    par$pareto.inds <- par$pareto.inds + 1 # set one-based indices for R
    return(par)
  } else
  {
    pareto.inds <- which.nondominated(-t(X.mat)) # new: use ecr package !! 
#    n = dim(X.mat)[1] # internal implementation 
#    if(is.null(n)) # here X.mat is a vector - only one vector 
#      return(list(pareto.X=X.mat, pareto.inds=1))
#    is.pareto = rep(0,n)
#    for(i in 1:n)
#      is.pareto[i] = is_pareto_optimal_rcpp(X.mat[i,], X.mat) # new: use rcpp   #   max(colMins(t(replicate(n, X[i,]))+epsilon - X, value=TRUE)) >= 0
#    pareto.inds <- which(as.logical(is.pareto))
    return(list(pareto.X=X.mat[pareto.inds,], pareto.inds=pareto.inds))
  }
}

# insert a new vector x if it is not dominated by any vector in X.mat. 
# Exclude all vectors dominated by x 
update_pareto_optimal_vecs <- function(X.mat, x)
{
  epsilon = 0.0000000000000000000000001
  if(isempty(X.mat))
    return(TRUE)
  if(is.null(dim(X.mat)))
    n.row <- 1
  else
    n.row <- dim(X.mat)[1]
  
  # next check if x dominates or is being dominated 
  
    return( min(rowMaxs(t(replicate(n.row, x))+epsilon - X.mat, value=TRUE)) >= 0 )
}



# Unite two lists of pareto-optinal vectors. Keep only pareto optimals in the joint list - 
union_pareto_optimal_vecs <- function(X.mat1, X.mat2, cpp.flag = TRUE)
{
  T <- dim(X.mat1)
  if(is.null(T))
    T <- 1
  else
    T <- T[2]
  if(cpp.flag)
  {
    union.X <- union_pareto_optimal_vecs_rcpp(matrix(X.mat1, ncol=T), matrix(X.mat2, ncol=T))
    union.X$pareto.inds1 = union.X$pareto.inds1+1  # change to one-based indices for R 
    union.X$pareto.inds2 = union.X$pareto.inds2+1   
  }
  else  
  {
    X.mat1 <- matrix(X.mat1, ncol=T)
    X.mat2 <- matrix(X.mat2, ncol=T)
    L1 <- dim(X.mat1)[1]
    L2 <- dim(X.mat2)[1]
    # New: unite the two and use the ecr package
    union.X <- c()
    pareto.inds <- which.nondominated(-t(rbind(X.mat1, X.mat2)))
    union.X$pareto.inds1 <- pareto.inds[pareto.inds <= L1]
    union.X$pareto.inds2 <- pareto.inds[pareto.inds > L1] - L1
    union.X$pareto.X <- rbind(X.mat1[union.X$pareto.inds1,], X.mat2[union.X$pareto.inds2,])
#    n1 = dim(X.mat1)[1]
#    n2 = dim(X.mat2)[1]
#    if(is.null(n1)) # here X.mat is a vector - only one vector 
#      return(list(pareto.X=X.mat2, pareto.inds=2:(n2+1)))
  
    # here still need to implement
  }
  return(union.X)
  
}
  
  


# Compute a vector of coordinate-wise lipschitz constants
compute_lipschitz_const <- function(loss.type, loss.params)
{
  T <- length(loss.params$theta)
  lip <- rep(0, T) # vector of coordinate wise lipschitz constants
  #  X.c <- compute_X_C_mat(X, C)
  if(loss.type == "quant")
    lip <- loss.params$theta
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))
  {
    lip <- loss.params$theta
  }
  
  if(loss.type == 'disease')  # weighted disease probability 
  {
    lip <- loss.params$theta/sqrt(2*pi)  # constant doesn't depend on prevalence 
  }
  return (lip)
}

# compute the maximal contribution of each vector to the loss
get_tensor_lipshitz_params <- function(X, loss.type, loss.params)  
{
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  lip <- compute_lipschitz_const(loss.type, loss.params)
    
  X.loss.mat <- matrix(rep(0, M*C), nrow=M, ncol=C)
  lip.pos.mat <- X.loss.mat
  lip.neg.mat <- X.loss.mat
  for(i in 1:M)
    for(j in 1:C)
    {
      X.loss.mat[i,j] <- loss_PS(X[i,j,], loss.type, loss.params)
      lip.pos.mat[i,j] <- pmax(X[i,j,], 0) %*% lip
      lip.neg.mat[i,j] <- -pmin(X[i,j,], 0) %*% lip
    }
  return(list(X.loss.mat=X.loss.mat, lip.pos.mat=lip.pos.mat, lip.neg.mat=lip.neg.mat))
}

# Check if loss function is monotonic
is_monotone_loss <- function(loss.type)
{
  return ( loss.type %in% c("disease", "quant") )
}
  

# Compute pareto probability for binary vectors (consider ties)
pareto_P_binary <- function(n, k, strong = FALSE)
{
  if(strong)
    return( sum( dbinom(0:k, k, 0.5) * (1 - 0.5^c(0:k))^(n-1)) )
  else
    return( sum( dbinom(0:k, k, 0.5) * (1 - 0.5^c(0:k) * (1 - 0.5^seq(k, 0, -1)))^(n-1)) )
}


# Compute pareto probability for binary vectors for a matrix of n and k values (consider ties)
pareto_P_binary_mat <- function(max.n, max.k, strong = FALSE)
{
  P.mat <- matrix(0, nrow=max.n, ncol=max.k)
  # for k == 1
  if(strong)
    P.mat[, 1] <- 0.5 ^ c(1:max.n) 
  else
    P.mat[, 1] <- 1 - 0.5 * (1 - 0.5 ^ c(0:(max.n-1))) 
  
  for(k in c(2:max.k))
  {
    binom.prob.vec <- dbinom(0:k, k, 0.5)
    if(strong)
        cond.prob.mat <- exp( log( (1 - 0.5^c(0:k))) %*% t(c(0:(max.n-1))) )
    else
        cond.prob.mat <- exp( log(1 - 0.5^c(0:k) * (1 - 0.5^seq(k, 0, -1))) %*% t(c(0:(max.n-1))) )
    
    P.mat[, k] <- colSums(binom.prob.vec * cond.prob.mat)
  }
  P.mat[1,] <- 1 # fix log-error for strong Pareto
  return(P.mat)
}
  
  
  
  
  