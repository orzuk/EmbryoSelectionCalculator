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
  
  if(loss.type == "balancing")
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
  if(loss.type == "balancing")
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
  if(loss.type == "balancing")
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
# X - a 3rd-prder tensor
# v - a vector 
tensor_vector_prod <- function(X, A, ax=3)
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
  

# a function for optimizing the selection of chromosomes:
# Run the optimization to find the optimal C:
optimize_C_relax <- function(X, C.init, loss.C, loss.params)
{
  M = dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
  # set defaults 
  if(isempty(C.init))
  {
    print("SET INIT")
    C.init <- matrix(1/C, M, C)
    print(dim(C.init))
  }
  if(!("mu.init" %in% names(loss.params)))  
    loss.params$mu.init <- 0.01 
  if(!("decay" %in% names(loss.params)))  
    loss.params$decay <- "inverse" 
  if(!("beta" %in% names(loss.params)))  
    loss.params$beta <- 0.9 
  if(!("epsilon" %in% names(loss.params)))  
    loss.params$epsilon <- 0.000001 
  if(!("max.iters" %in% names(loss.params)))  
    loss.params$max.iters <- 10000 
  
  loss.vec <- rep(0, loss.params$max.iters)
  loss.vec[1] <- loss_PS(compute_X_C_mat(X, C.init), loss.C, loss.params)
  delta.loss <- 999999
  mu.t <- loss.params$mu.init
  C.cur <- C.init
  t <- 1
  print("start while")
  while((abs(delta.loss) > loss.params$epsilon) & (t<=loss.params$max.iters)) # Projected gradient descent
  {
    if(t%%50 == 0)
      print(paste("while t=", t))
    C.next <- C.cur + mu.t * grad_loss_PS(X, C.cur, loss.C, loss.params)
#    print(paste("while project t=", t))
#    print(C.cur)
#    print("Grad:")
#    print(mu.t * grad_loss_PS(X, C.cur, loss.C, loss.params))
#    print("Next:")
#    print(C.next)
#    dim(C.next)
    C.cur <- project_stochastic_matrix(C.next)
    if(loss.params$decay == "exp") # exponential decay
      mu.t <- mu.t * loss.params$beta  # update step size 
    if(loss.params$decay == "inverse") # 1/t decay
      mu.t <- loss.params$mu.init / (1 + loss.params$beta*t)

    t <- t+1
    loss.vec[t] <- loss_PS(compute_X_C_mat(X, C.cur), loss.C, loss.params)
    delta.loss <- loss.vec[t] - loss.vec[t-1]    
        # otherwise constant learning rate 
  }
    
    
  c.vec <- max.col(C.cur) # Convert to zero-one 
  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.C, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_C_mat(X, C.cur)
      
  return(list(opt.X=opt.X, opt.loss=opt.loss, c.opt=c.vec, loss.vec=loss.vec[1:t]))
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
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]

  # Need to save also c-vec for each branch
    
  print("Start B&B")
  cur.X <- X[1,,]
  cur.c <- 1:C
  cur.X <- get_pareto_optimal_vecs(cur.X) # Save only Pareto-optimal vectors 
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  print(cur.X)
  print(L)
  print(paste("L=", L, " before ..."))
  
  L.vec <- rep(0, M)
  L.vec[1] = L
  print(paste("L=", L, " start loop"))
  for( i in c(2:M))
  {  
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    print(paste0("Start B&B i=", i))
    new.X <- c()
    new.c <- c()
    for(j in c(1:L))  # loop over all vectors in the current stack      
      for  (c in c(1:C))  # loop over possible vectors to add 
      {
        if(L == 1)
          v = cur.X+X[i,c,]
        else
          v = cur.X[j,]+X[i,c,]
#        print("Start if")
        if(is_pareto_optimal(v, new.X))
        {
          new.X <- rbind(new.X, v)
          new.c <- rbind(new.c, c(cur.c[j,], c) )
        }
      }
    cur.X <- new.X
    cur.c <- new.c
    L <- dim(new.X)[1]  
    L.vec[i] = L
    print(paste("Stack Size:", L))
  }
    
  # Finally find the cost-minimizer out ofthe Pareto-optimal vectors
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
  {
    L = 1
    loss.vec = loss_PS(cur.X, loss.C, loss.params)
  }  else
  {
    loss.vec <- rep(0, L)
    for(i in 1:L)
      loss.vec[i] <- loss_PS(cur.X[i,], loss.C, loss.params)
  }
  print(loss.vec)
  i.min <- which.min(loss.vec) # loss_PS(cur.X, loss.C, loss.params))
    
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X))
}

# A closed-form solution for the case of balancing selection 
optimize_C_balancing_exact <- function(X, loss.C, loss.params)
{
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
#  A <- grad_loss_PS(X, C, "balancing", loss.params)

  A <- matrix(0, nrow=M*C, ncol=M*C)
  for(k in c(1:T))
    A <- A + as.vector(X[,,k]) %*% t(as.vector(X[,,k]))
    
  E <- matrix(0, nrow=M, ncol=M*C)
  for(i in c(1:M))
    E[i,((i-1)*C+1):(i*C)] <- 1
  
  b <- c(rep(0, M*C), rep(1, M)) # free vector for linear system   
  
  Big.A <- rbind(cbind(A, t(E)), cbind(E, matrix(0, nrow=M, ncol=M)))
  
  v <- solve(Big.A, b) # Solve system
  
  c.vec <- max.col(matrix(v[1:(M*C)], nrow=M, ncol=C)) # Convert to matrix and take max of each row 
  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.C, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec)
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, c.opt=c.vec))
  
}
  