# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)
library(tensor)
library(olpsR) # for projection onto the simplex remotes::install_github("ngloe/olpsR")

source('chrom_select_funcs.R')

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
    loss.params$mu.init <- 1 
  if(!("decay" %in% names(loss.params)))  
    loss.params$decay <- "inverse" 
  if(!("beta" %in% names(loss.params)))  
  {
    if(loss.params$decay == "inverse")
      loss.params$beta <- 0.01 
    else 
      loss.params$beta <- 0.99 
  }
  if(!("epsilon" %in% names(loss.params)))  
    loss.params$epsilon <- 0.000001 
  if(!("max.iters" %in% names(loss.params)))  
    loss.params$max.iters <- 10000 
  
  loss.vec <- rep(0, loss.params$max.iters)
  loss.vec[1] <- loss_PS(compute_X_C_mat(X, C.init), loss.C, loss.params)
  delta.loss <- 999999
  mu.t <- loss.params$mu.init
  mu.t.vec <- rep(0, loss.params$max.iters)
  mu.t.vec[1] <- mu.t
  C.cur <- C.init
  t <- 1
  print("start while")
  while((abs(delta.loss) > loss.params$epsilon) & (t<=loss.params$max.iters)) # Projected gradient descent
  {
    if(t%%50 == 0)
      print(paste("while t=", t))
    C.next <- C.cur - mu.t * grad_loss_PS(X, C.cur, loss.C, loss.params) # move in minus gradient direction (minimization) 
    C.cur <- project_stochastic_matrix(C.next)
    if(loss.params$decay == "exp") # exponential decay
      mu.t <- mu.t * loss.params$beta  # update step size 
    if(loss.params$decay == "inverse") # 1/t decay
      mu.t <- loss.params$mu.init / (1 + loss.params$beta*t)
    
    t <- t+1
    mu.t.vec[t] <- mu.t
    loss.vec[t] <- loss_PS(compute_X_C_mat(X, C.cur), loss.C, loss.params)
    delta.loss <- loss.vec[t] - loss.vec[t-1]    
    # otherwise constant learning rate 
  }
  
  
  c.vec <- max.col(C.cur) # Convert to zero-one 
  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.C, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec) #  opt.X <- compute_X_C_mat(X, C.cur)
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, loss.vec=loss.vec[1:t], mu.t.vec=mu.t.vec[1:t], 
              C.mat=C.cur))
}



# Compute for block vectors (not i.i.d.)
# here: n = C^M
pareto_P_block <- function(C, M, k, iters=1000)
{
  print("start pareto P block")
  n <- C^M
  n.pareto <- 0
  Sigma.T <- eye(k) # No correlations between siblings 
  Sigma.K <- 0.5*diag(C) + matrix(0.5, nrow=C, ncol=C)   # kinship-correlations matrix 
  sigma.blocks <- ones(M, 1)
  loss.params <- c()
  loss.params$theta <- ones(k, 1)
  for(i in 1:iters)
  {
    print("Sim X")
    X = simulate_PS_chrom_disease_risk(M, C, k, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
    print("Solve B&B")
    sol.bb <- optimize_C_branch_and_bound(X, "quant", loss.params) # run B&B. Loss at the end doesn't matter. 
    n.pareto <- n.pareto + length(sol.bb$loss.vec)
  }
  return(n.pareto / (iters*n))
}



# Choose in a greedy manner the best X  
optimize_C_quant <- function(X, loss.C, loss.params)
{
  M <- dim(X)[1]
  T <- dim(X)[3]
  c.vec <- rep(0, M)
  opt.x <- rep(0, T)
  for( i in c(1:M))
  {
    v <- sum(X[i,,] * loss.params$theta)  
    c.vec[i] <- which.max(rowSums(X[i,,] * loss.params$theta)  )
    opt.x <- opt.x + X[i,c.vec[i],]
  }
  return(list(opt.X = opt.x, opt.c = c.vec, opt.loss = sum(opt.x*loss.params$theta)))
}



# A branch and bound algorithm 
optimize_C_branch_and_bound <- function(X, loss.C, loss.params)
{
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
  # Need to save also c-vec for each branch
  
  #  print("Start B&B")
  cur.X <- X[1,,]
  #  cur.c <- 1:C
  par.X <- get_pareto_optimal_vecs(cur.X) # Save only Pareto-optimal vectors . Needs fixing 
  cur.c <- t(t(par.X$pareto.inds))
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  L.vec <- rep(0, M)
  L.vec[1] = L
  #  print(paste("L=", L, " start loop"))
  for( i in c(2:M))
  {  
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    new.X <- c()
    new.c <- c()
    for(j in c(1:L))  # loop over all vectors in the current stack      
      for  (c in c(1:C))  # loop over possible vectors to add 
      {
        if(is.null(dim(cur.X)))
          v = cur.X+X[i,c,]
        else
          v = cur.X[j,]+X[i,c,]
        if(is_pareto_optimal(v, new.X))
        {
          new.X <- rbind(new.X, v)
          if(is.null(dim(cur.c)))
            new.c <- rbind(new.c, c(cur.c[j], c) )
          else
            new.c <- rbind(new.c, c(cur.c[j,], c) )
        }
      }
    cur.X <- new.X
    cur.c <- new.c
    L.vec[i] = dim(new.X)[1]  
    if(i == M)
      print(paste0("B&B C=", C, " i=", i, " out of ", M, " Stack Size:", dim(new.X)[1]))
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
  i.min <- which.min(loss.vec) # find vector minimizing loss 
  
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X))
}


###############################################################
# A closed-form solution for the case of stabilizing selection 
#
# Input: 
# X - tensor of polygenic scores 
# loss.C - string signifying loss type
# loss.params - parameters of the loss function
#
# Output: 
# opt.X - optimal X 
# opt.loss - 
# .c.opt - optimal value of the loss 
###############################################################
optimize_C_stabilizing_exact <- function(X, loss.C, loss.params)
{
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  #  A <- grad_loss_PS(X, C, "stabilizing", loss.params)
  
  A <- matrix(0, nrow=M*C, ncol=M*C)
  for(k in c(1:T))
    #    A <- A + 2 * loss.params$theta[k] *  as.vector((X[,,k])) %*% t(as.vector((X[,,k])))
    A <- A + 2 * loss.params$theta[k] *  as.vector(t(X[,,k])) %*% t(as.vector(t(X[,,k])))
  E <- matrix(0, nrow=M, ncol=M*C)
  for(i in c(1:M))
    E[i,((i-1)*C+1):(i*C)] <- 1
  b <- c(rep(0, M*C), rep(1, M)) # free vector for linear system   
  Big.A <- rbind(cbind(A, t(E)), cbind(E, matrix(0, nrow=M, ncol=M)))
  
  #  return(list(  Big.A=Big.A, b=b)) # temp debug
  
  if(T+2*M >= (C+1)*M)
    v <- solve(Big.A, b) # Solve system
  else  
    v <- pinv(Big.A) %*% b  # infinite solutions. Use pseudo-inverse
  C.mat <- matrix(v[1:(M*C)], nrow=M, ncol=C, byrow = TRUE)
  loss.mat <- loss_PS(compute_X_C_mat(X, C.mat), loss.C, loss.params)
  c.vec <- max.col(C.mat) # Convert to matrix and take max of each row 
  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.C, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec)
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, C.mat=C.mat, loss.mat=loss.mat, 
              Big.A=Big.A, b=b))
}


# A wrapper function for all optimizations



# Compute average gain using simulations 
compute_gain_sim <- function(params, loss.C, loss.params)
{
  for (t in 1:params$iters)
  {
    X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
    print("Solve B&B")
    
    if(loss.C == "quant")
      sol <- optimize_C_branch_and_bound(X, "quant", loss.params) # run B&B. Loss at the end doesn't matter. 
    if(loss.C == "disease")
      sol <- optimize_C_branch_and_bound(X, "disease", loss.params) # run B&B. Loss at the end doesn't matter. 
    if(loss.C == "stabilizing")
      sol <- optimize_C_stabilizing_exact(X, "stabilizing", loss.params) # run B&B. Loss at the end doesn't matter. 
    
    
    
  }
  
  return(gain)
}


