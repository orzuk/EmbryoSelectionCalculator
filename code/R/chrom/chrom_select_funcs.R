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

source("pareto_funcs.R")
Rcpp::sourceCpp("cpp/chrom_funcs.cpp")  # fast functions  


chr.lengths <- c(0.0821,0.0799,0.0654,0.0628,0.0599,0.0564,0.0526,0.0479,0.0457,0.0441,
                                     0.0446,0.0440,0.0377,0.0353,0.0336,0.0298,0.0275,0.0265,0.0193,0.0213,0.0154,0.0168,0.0515)

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)


###############################################################
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
###############################################################
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



###############################################################
# The loss 
# Input: 
# X.c - a vector of length C
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# loss - a scalar 
###############################################################
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



###############################################################
# The loss for vectors 
# Input: 
# X.c.mat - a matrix with multiple vectors of length C
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# loss.vec - a vector of losses for each row of X.c.mat  
###############################################################
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



###############################################################
# The loss for vectors 
# Input: 
# X - 3rd order tensor
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# upperbound - an upper-bound on the optimal loss  
# lowerbound - an lower-bound on the optimal loss  
###############################################################
bound_loss_PS_mat <- function(X, loss.type, loss.params)
{
  M <- dim(X)[1]
  T <- dim(X)[3]
  max.X <- rep(0, T)
  min.X <- rep(0, T)
  for(b in 1:M) # loop on blocks
  {
    max.X <- max.X + colMaxs(X[b,,], value = TRUE) # Get maximum at each coordinate 
    min.X <- min.X + colMins(X[b,,], value = TRUE) # Get minimum at each coordinate 
  }
  upper = loss_PS(min.X, loss.type, loss.params)
  lower = loss_PS(max.X, loss.type, loss.params)
  return(list(upperbound=upper, lowerbound=lower, max.X=max.X, min.X=min.X))
}  



###############################################################
# The gradient for the loss.
# Return a matrix of size: M*C
# Can also add a regularizer 
###############################################################
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

###############################################################
# Lipschitz constant for a loss
# Input: 
# X.c.mat - a matrix with multiple vectors of length C
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# alpha.vec - a vector of lipschitz constants for each input coordinate of the loss  
#
###############################################################

lipschitz_loss_PS <- function(loss.type, loss.params)
{
  if(loss.type == "quant")
    alpha.vec <- loss.params$theta
  
  if((loss.type == "stabilizing") || (loss.type == "balancing"))  # Global constant is infinity. May get local constants with bound on X
    alpha.vec <- rep(Inf, len(loss.params$theta))    

  if(loss.type == 'disease')  # weighted disease probability 
    alpha.vec <- loss.params$theta / sqrt(2*pi*(1-loss.params$h.ps))

  return(alpha.vec)  
}

# Multiply a tensor by matrix 
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


#############################################
# Multiply a tensor by vector
# Input: 
# X - a 3rd-order tensor
# v - a vector 
# ax - on which axis to run the summation of multiplication 
# Output: 
# A matrix R such that R_{i,j} = \sum_k X_{ijk} v_k
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
  


# Compute a vector of coordinate-wise lipschitz constants for a loss function
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


#############################################
# compute the maximal contribution of each vector to the loss
# Input: 
# X - a 3rd-prder tensor of dims: M*C*T
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# X.loss.mat - matrix of losses for each scores vector (for given block and copy)
# lip.pos.mat - matrix of A+ values
# lip.neg.mat - matrix of A- values
#############################################
get_tensor_lipschitz_params <- function(X, loss.type, loss.params)  
{
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  lip <- compute_lipschitz_const(loss.type, loss.params)  # vector of lipschitz constants
    
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
  


  