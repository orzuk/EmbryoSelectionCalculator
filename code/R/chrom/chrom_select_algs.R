# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)
library(tensor)
library(olpsR) # for projection onto the simplex remotes::install_github("ngloe/olpsR")
library(quadprog)  # for stabilizing selection loss 
library(CVXR)
library(Rcsdp)  # For semidefinite programming


Rcpp::sourceCpp("cpp/chrom_funcs.cpp")  # fast functions  
source('chrom_select_funcs.R')



###############################################################
# Optimize the selection of chromosomes:
# Run the optimization to find the optimal C:
# perform gradient-descent for the relaxed problem. Finally, round results
# Input: 
# X - a 3D tensor of block genetic scores
# C.init - initial condition (M*C matrix of probabilities)
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# 
# Output: 
# list with optimal choices 
###############################################################
optimize_C_relax <- function(X, C.init, loss.type, loss.params)
{
  print("Start optimize relax")
  if(!("eta" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$eta <- 0 
  
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  
  # set defaults 
  if(isempty(C.init))
    C.init <- matrix(1/C, M, C)
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
  if(!("epsilon" %in% names(loss.params)))  # default tolerance 
    loss.params$epsilon <- 0.000001 
  if(!("max.iters" %in% names(loss.params)))  
    loss.params$max.iters <- 10000 
  
  loss.vec <- rep(0, loss.params$max.iters)
  loss.vec[1] <- loss_PS(compute_X_C_mat(X, C.init), loss.type, loss.params) - loss.params$eta * sum(C.init**2)
  delta.loss <- 999999
  mu.t <- loss.params$mu.init
  mu.t.vec <- rep(0, loss.params$max.iters)
  mu.t.vec[1] <- mu.t
  C.cur <- C.init
  t <- 1
#  print("start while")
  while((abs(delta.loss) > loss.params$epsilon) & (t<=loss.params$max.iters)) # Run projected gradient descent
  {
    if(t%%100 == 0)
      print(paste("gradient t=", t, " loss=", loss.vec[t]))
    cur.grad <- grad_loss_PS(X, C.cur, loss.type, loss.params)
    cur.grad <- cur.grad / max(1, max(abs(cur.grad)))
    C.next <- C.cur - mu.t * cur.grad # grad_loss_PS(X, C.cur, loss.type, loss.params) # Gradient descent: move in minus gradient direction (minimization) 
    C.cur <- project_stochastic_matrix(C.next) # Project onto the simplex 
    if(loss.params$decay == "exp") # exponential decay
      mu.t <- mu.t * loss.params$beta  # update step size 
    if(loss.params$decay == "inverse") # 1/t decay
      mu.t <- loss.params$mu.init / (1 + loss.params$beta*t)
    # otherwise constant learning rate 
    
    t <- t+1
    mu.t.vec[t] <- mu.t
    
    loss.vec[t] <- loss_PS(compute_X_C_mat(X, C.cur), loss.type, loss.params) - loss.params$eta * sum(C.cur**2)
    delta.loss <- loss.vec[t] - loss.vec[t-1]    
  }
  
#  c.vec <- max.col(C.cur) # Rounding stage: convert to zero-one 
#  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.type, loss.params)  # cost of the rounded solution 
#  opt.X <- compute_X_c_vec(X, c.vec) #  opt.X <- compute_X_C_mat(X, C.cur)

  ret <- real_to_integer_solution(X, C.cur, loss.type, loss.params)
  ret$loss.vec <- loss.vec[1:t]
  ret$mu.t.vec <- mu.t.vec[1:t]
  ret$C.mat <- C.cur
  return(ret)
      
#  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, loss.vec=loss.vec[1:t], mu.t.vec=mu.t.vec[1:t], 
#              C.mat=C.cur))
}


###############################################################
# Helper function: convert real solution to integer solution
# here: n = C^M
# Input: 
# X - a 3D tensor of block genetic scores
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# Output: 
# opt.X - Resulting optimal X_c
# opt.loss - loss of optimal solution
# opt.c - binary solution in the simplex product
###############################################################
real_to_integer_solution <- function(X, C.mat, loss.type, loss.params)
{
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  if(!("n.candidates" %in% names(loss.params)))  
    loss.params$n.candidates <- 1000 # generate random candidates

  # First project onto the simplex: 
  C.mat.in <- project_stochastic_matrix(C.mat)
  c.vecs <- matrix(0, nrow = M, ncol = loss.params$n.candidates)
  for(i in 1:M)  # Sample according t weights 
    c.vecs[i,] <- sample(1:C, loss.params$n.candidates, replace = TRUE  , prob=C.mat.in[i,])

  i.min <- which.min(loss_PS_mat(compute_X_c_vecs(X, t(c.vecs)), loss.type, loss.params))
  c.vec <- c.vecs[,i.min]
#  c.vec <- max.col(C.mat) # Rounding stage: convert to zero-one 
  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.type, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec) #  opt.X <- compute_X_C_mat(X, C.cur)
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec))
}  



###############################################################
# Helper function: convert SDP real solution to integer solution
# here: n = C^M
# Input: 
# X - a 3D tensor of block genetic scores
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# Output: 
# opt.X - Resulting optimal X_c
# opt.loss - loss of optimal solution
# opt.c - binary solution in the simplex product
###############################################################
SDP_to_integer_solution <- function(X, C.mat.pos.def, loss.type, loss.params, method = "svd")
{
  M <- dim(X)[1]
  C <- dim(X)[2]
  T <- dim(X)[3]
  print(dim(X))
  
  rand.iters <- 200
  ctr <- 0
  if(method == "svd")
    rand.iters <- 2
  c.vecs <- matrix(0, nrow = rand.iters, ncol = M)
  if(method == "randomization")
  {
    z = mvrnorm(n = rand.iters, mu = rep(0, dim(C.mat.pos.def)[1]), Sigma = C.mat.pos.def)
    for(i in 1:M)
      c.vecs[,i] <- apply(z[,((i-1)*C+1):(i*C) ], 1, FUN = which.max)
    ctr <- rand.iters-2
  }
  # Always use svd!! 
  SDR.svd <- svd(C.mat.pos.def, 1, 1) # take the best rank-1 approximation
  c.vecs[ctr+1,] <- apply(matrix(head(SDR.svd$u, -1), nrow=M, byrow=TRUE), 1, FUN = which.max) # can be also min!!! 
  c.vecs[ctr+2,] <- apply(matrix(head(SDR.svd$u, -1), nrow=M, byrow=TRUE), 1, FUN = which.min) # can be also min!!! 
  
#  {
#    rand.iters <- 2
#    SDR.svd <- svd(C.mat.pos.def, 1, 1) # take the best rank-1 approximation
#    c.vecs <- matrix(0, nrow = M, ncol = 2)
#    c.vecs[,1] <- apply(matrix(head(SDR.svd$u, -1), nrow=M, byrow=TRUE), 1, FUN = which.max) # can be also min!!! 
#    c.vecs[,2] <- apply(matrix(head(SDR.svd$u, -1), nrow=M, byrow=TRUE), 1, FUN = which.min) # can be also min!!! 
#  }  
  print(c.vecs)
  i.min <- which.min(loss_PS_mat(compute_X_c_vecs(X, (c.vecs)), loss.type, loss.params))
  print("i.min:")
  print(i.min)
  c.vec <- c.vecs[i.min,]  # get best random vector
  print("c.vec:")
  print(c.vec)
  C.mat <- matrix(0, nrow=M, ncol=C)
  C.mat[cbind(1:M, c.vec)] <- 1

  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.type, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec) #  opt.X <- compute_X_C_mat(X, C.cur)
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, C.mat=C.mat))
}  



###############################################################
# Compute the probability of being Pareto-optimal for block vectors (not i.i.d.)
# here: n = C^M
# Input: 
# C,M,T - dimensions 
# iters - number of simulations to run 
# Output: 
# Empirical Pareto-front probability 
###############################################################
pareto_P_block <- function(C, M, T, iters=1000)
{
  print("start pareto P block")
  n <- C^M
  n.pareto <- 0
  Sigma.T <- eye(k) # No correlations between siblings 
  Sigma.K <- 0.5*diag(C) + matrix(0.5, nrow=C, ncol=C)   # kinship-correlations matrix 
  sigma.blocks <- ones(M, 1)
  loss.params <- c()
  loss.params$theta <- ones(T, 1)
  for(i in 1:iters)
  {
    X = simulate_PS_chrom_disease_risk(M, C, T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, T))
    sol.bb <- optimize_C_branch_and_bound(X, "quant", loss.params) # run B&B. Loss at the end doesn't matter. 
    n.pareto <- n.pareto + length(sol.bb$loss.vec)
  }
  return(n.pareto / (iters*n))
}



###############################################################
# Choose in a greedy manner the best X. Mininizing the additive quantitative loss   
# Input: 
# X - a 3D tensor of block genetic scores
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# Output: 
# list with optimal choices:
# - opt.X 
# - opt.c 
# - opt.loss 
###############################################################
optimize_C_quant <- function(X, loss.type, loss.params)
{
  M <- dim(X)[1]
  T <- dim(X)[3]
  c.vec <- rep(0, M)
  opt.x <- rep(0, T)
  for(i in c(1:M))
  {
    v <- sum(X[i,,] * loss.params$theta)  
    c.vec[i] <- which.min(rowSums(X[i,,] * loss.params$theta)  )
    opt.x <- opt.x + X[i,c.vec[i],]
  }
  return(list(opt.X = opt.x, opt.c = c.vec, opt.loss = sum(opt.x*loss.params$theta)))
}



###############################################################
# A branch and bound algorithm for finding the X combination with minimal loss for monotone loss 
# Run the optimization to find the optimal C:
# Input: 
# X - a 3D tensor of block genetic scores
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# 
# Output: 
# list with optimal choices:
# - opt.X 
# - opt.c 
# - opt.loss 
# - loss.vec 
# - L.vec 
# - pareto.opt.X 
# - pareto.opt.c
###############################################################
optimize_C_branch_and_bound <- function(X, loss.type, loss.params)
{
  if(!("cpp" %in% names(loss.params)))  
    loss.params$cpp <- FALSE  # default: run in R  

  #  print("Start optimize B&B in R") 
  if(length(dim(X)) == 2) # M == 1) # Nothing to optimize here!
  {
    M <- 1; C <- dim(X)[1]; T <- dim(X)[2]
#    print(paste0("M, C, T:", M, " ", C, " ", T))
    loss.vec <- loss_PS_mat(X, loss.type, loss.params) # TO FILL! !! 
    i.min <- which.min(loss.vec)
    par.X <- get_pareto_optimal_vecs(X) # Save only Pareto-optimal vectors. Needs fixing 
#    print("Pareto X:")
#    print(par.X$pareto.X)
#    print("L.vec, Pareto C:")
#    print(par.X$pareto.inds)
    loss.vec <- loss.vec[par.X$pareto.inds]
    if(loss.params$cpp) # new: run in cpp
      par.X$pareto.inds <- par.X$pareto.inds - 1
    return(list(opt.X = X[i.min,], opt.c = i.min, opt.loss = min(loss.vec), 
                loss.vec = loss.vec, L.vec = C, 
                pareto.opt.X= par.X$pareto.X, pareto.opt.c = t(t(par.X$pareto.inds))))
  }
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
#  print(paste0("M, C, T:", M, " ", C, " ", T))
  
  if(loss.params$cpp) # new: run in cpp
  {
#    print("Start optimize B&B CPP") 
    return(optimize_C_branch_and_bound_rcpp(X, loss.type, loss.params))
  }
  if(!("lipschitz" %in% names(loss.params)))
    loss.params$lipschitz <- FALSE  # Default: no lipschitz constant (vector)
#  print("Start optimize B&B R") 
  
  
  par.X <- get_pareto_optimal_vecs(X[1,,]) # Save only Pareto-optimal vectors. Needs fixing 
  cur.c <- t(t(par.X$pareto.inds))  # Format matrix 
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  if(loss.params$lipschitz)
  {
    A.plus <- rowMaxs(tensor_vector_prod(pmax(X,0), loss.params$lipschitz.alpha, 3), value = TRUE)
    A.minus <- rowMins(tensor_vector_prod(pmin(X,0), loss.params$lipschitz.alpha, 3), value = TRUE)  
    A.plus.cum <- rev(cumsum(rev(A.plus)))
    A.minus.cum <- rev(cumsum(rev(A.minus)))
  }
    
  L.vec <- rep(0, M)
  L.vec[1] = L
  #  print(paste("L=", L, " start loop"))
  for(i in 2:M) # loop on block
  {
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    #    if(i == M)
    print(paste0("B&B i=", i, " L=", L))
    
    # We know that the first vectors are Pareto-optimal
    if(L>1)
    {
      new.X <- sweep(cur.X, 2, X[i,1,], "+")
    } else
      new.X <- matrix(cur.X + X[i,1,], nrow=1) # check that it doesn't flip x
    new.c <- cbind(cur.c, rep(1, L) )
    # new version: create sums and take union
    for(c in 2:C)
    {
      if(L>1)
        temp.X <- sweep(cur.X, 2, X[i,c,], "+")
      else
        temp.X <- cur.X  + X[i,c,]
      
      #      ttt <- Sys.time()
      #      temp.X <- get_pareto_optimal_vecs(temp.X) # check which is time-consuming. No need to check! they're already pareto-optimal!
      #      get.time <- difftime(Sys.time() , ttt, units="secs") 
      #      print(paste0("get time: (sec.): ", get.time))
      ttt <- Sys.time()
      union.X <- union_pareto_optimal_vecs(new.X, temp.X) # $pareto.X)
      union.time <- difftime(Sys.time() , ttt, units="secs") 
      print(paste0("union time: (sec.): ", union.time))
      print(paste0("Num pareto1: ", length(union.X$pareto.inds1), " Num pareto2: ", length(union.X$pareto.inds2)))
      #      ttt <- Sys.time()
    #      union.X <- union_pareto_optimal_vecs(new.X, temp.X, 2) # $pareto.X)
      #      union.time <- difftime(Sys.time() , ttt, units="secs") 
      #      print(paste0("union time2: (sec.): ", union.time))
      #      print(paste0("Num pareto1: ", length(union.X$pareto.inds1), " Num pareto2: ", length(union.X$pareto.inds2)))
      #      ttt <- Sys.time()
      #      union.X <- union_pareto_optimal_vecs(new.X, temp.X, FALSE) # new: run in C!
      #      union.time.ecr <- difftime(Sys.time() , ttt, units="secs") 
      #      print(paste0("union time ecr: (sec.): ", union.time.ecr))
      #      print(paste0("Num pareto1: ", length(union.X$pareto.inds1), " Num pareto2: ", length(union.X$pareto.inds2)))
      
      new.X <- union.X$pareto.X
      add.c <- cbind(cur.c, rep(c, L))
      new.c <- rbind(new.c[union.X$pareto.inds1,], add.c[union.X$pareto.inds2,]) # need to modify here indices 
    }
    
    if(is.null(dim(new.X)[1])) L.vec[i] = 1 else L.vec[i] = dim(new.X)[1]  
    # New: get-rid of some vectors that can't develop into optimal ones: 
    if(loss.params$lipschitz)
    {
#      new.v <- sweep( B[[b+1]]$pareto.opt.X, 2, B[[b]]$pareto.opt.X[j,], "+")
      new.loss.vec <- loss_PS_mat_rcpp(new.X, loss.type, loss.params)
      print(paste("MAX LOSS:", max(new.loss.vec), "MIN LOSS:", min(new.loss.vec), "Diff: ", max(new.loss.vec)-min(new.loss.vec)))
      print(paste("A.plus:", A.plus.cum[i], "A.minus:", A.minus.cum[i]))
      good.inds <- which(new.loss.vec <= min(new.loss.vec) + A.plus.cum[i] - A.minus.cum[i])  
      new.c <- new.c[good.inds, ]
      new.X <- new.X[good.inds, ]
    }
        
    cur.X <- new.X
    cur.c <- new.c
    
    #    if(i == M)
    #      print(paste0("B&B C=", C, " i=", i, " out of ", M, " Stack Size:", dim(new.X)[1]))
  } # end loop on blocks 
  
  # Finally find the cost-minimizer out of the Pareto-optimal vectors
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
  {
    L = 1
    loss.vec = loss_PS(cur.X, loss.type, loss.params)
  }  else
  {
    #    loss.vec <- rep(0, L)
    #    for(i in 1:L)
    #      loss.vec[i] <- loss_PS(cur.X[i,], loss.type, loss.params)
    loss.vec <- loss_PS_mat(cur.X, loss.type, loss.params) # use cpp if cpp flag is on 
    #    print("Max error:")
    #    print(max(abs(loss.vec-loss.vec2)))
  }
  i.min <- which.min(loss.vec) # find vector minimizing loss 
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), 
              loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X, pareto.opt.c = cur.c))
}



###############################################################
# A branch and bound algorithm for finding the X combination with minimal loss for Lipschitz loss functions
# Run the optimization to find the optimal C:
# Input: 
# X - a 3D tensor of block genetic scores
# loss.type - type of loss function
# loss.params - parameters determining the loss function
# 
# Output: 
# list with optimal choices:
# - opt.X 
# - opt.c 
# - opt.loss 
# - loss.vec 
# - L.vec 
# - pareto.opt.X 
# - pareto.opt.c
###############################################################
optimize_C_branch_and_bound_divide_and_conquer <- function(X, loss.type, loss.params)
{
  # Add timing: 
  start.time <- Sys.time()
#  print("Start B&B Div. Conq. in R !!! ")
  M <- dim(X)[1];   C <- dim(X)[2];  T <- dim(X)[3]
  if(!("n.blocks" %in% names(loss.params)))  
    loss.params$n.blocks <- 2  # default: divide to two blocks 
  loss.params$n.blocks <- min(loss.params$n.blocks, M)  # can't divide to more blocks  ! 
  
#  lip.alpha <- compute_lipschitz_const(loss.type, loss.params)
#  lip.tensor <- get_tensor_lipschitz_params(X, loss.type, loss.params)  # get loss and bounds for individual vectors (not implemented?)
  L.vec <- rep(0, M)
  
  # Divide to blocks: 
  M.vec <- rep(floor(M/loss.params$n.blocks), loss.params$n.blocks)
  if( mod(M, loss.params$n.blocks)>0 )
    M.vec[1:mod(M, loss.params$n.blocks)] <- M.vec[1:mod(M, loss.params$n.blocks)] + 1 # correct to sum to M
  M.vec.cum <- c(0, cumsum(M.vec))
  B <- vector("list", loss.params$n.blocks) 
  n.pareto <- rep(0,  loss.params$n.blocks)
  opt.X.upperbound <- rep(0, T) # A 'good' vector with low score that is an upper-bound for L*
#  opt.c.upperbound <- c() # not needed
  for(b in 1:loss.params$n.blocks) # get Pareto-optimal vectors for each block 
  {
    B[[b]] <- optimize_C_branch_and_bound(X[(M.vec.cum[b]+1):M.vec.cum[b+1],,], loss.type, loss.params)  # compute Pareto optimal vectors for block
    opt.X.upperbound <- opt.X.upperbound + B[[b]]$opt.X  #    opt.c.upperbound <- c(opt.c.upperbound,  B[[b]]$opt.c)  # not needed !!! 
    B[[b]]$max.X <- colMaxs(B[[b]]$pareto.opt.X, value = TRUE) # Get maximum at each coordinate 
    n.pareto[b] <-  length(B[[b]]$loss.vec)  # number of vectors in each block
    
#    print("B-pareto-c-start:")
#    print(B[[b]]$pareto.opt.c)
    
  }
#  print("Sub-blocks # pareto-optimal vecs:")
#  print(n.pareto)

  start.bounds <- bound_monotone_loss_pareto_blocks_PS_mat(B, loss.type, loss.params)
#  print(paste0("Bounds: [", start.bounds$lowerbound, ", ", start.bounds$upperbound, "], eps=", start.bounds$upperbound-start.bounds$lowerbound))    
    
#  print("B-pareto-c-bli-filter:")
#  print(B[[1]]$pareto.opt.c)
  B <- filter_solutions(B, loss.type, loss.params)
  
  n.pareto <- B$n.pareto
  B <- B$sol
#  print("B-pareto-c-im-filter:")
#  print(B[[1]]$pareto.opt.c)
  
  L.vec[(M.vec.cum[1]+1):M.vec.cum[2]] = B[[1]]$L.vec  # update start 
#  print("L.vec now:")
#  print(L.vec)
  L.upperbound <- loss_PS(opt.X.upperbound, loss.type, loss.params) + 0.00000000001  # > Loss_*
  new.L.upperbound = L.upperbound
  bb.time <- difftime(Sys.time() , start.time, units="secs") 
  
  print(paste0("cpp=", loss.params$cpp, " b&b time (sec.): ", bb.time))
  
  ##############################################
  # Next loop from one side and merge blocks:
  ##############################################
  if(loss.params$n.blocks == 1)
  {
    run.blocks <- c()
    new.X <- B[[b]]$pareto.opt.X
    new.c <- B[[b]]$pareto.opt.c
  } else
    run.blocks <- 1:(loss.params$n.blocks-1)  # already discard last block
  
  for(b in run.blocks) # Now start to merge blocks 
  {
    block.start.time <- Sys.time()
    new.opt.X.upperbound = B[[b]]$opt.X
    max.X <- rep(0, T)
    for(b2 in c((b+1):loss.params$n.blocks))
    {
      max.X <- max.X + B[[b2]]$max.X # This reduction can be applied on both sides 
      new.opt.X.upperbound <- new.opt.X.upperbound + B[[b2]]$opt.X 
    }
    new.max.X = max.X - B[[b+1]]$max.X
    new.L.upperbound <- min(new.L.upperbound, loss_PS(new.opt.X.upperbound, loss.type, loss.params) + 0.00000000001)  # > Loss_*
    
#    print(paste0("Global upperbound=", L.upperbound, ", New upperbound=", new.L.upperbound))
    
    #    B[[b]]$L.lowerbound.vec <- rep(0, n.pareto[b])   # lower-bound on the final score L* of every selection of current X 
#    for(i in 1:n.pareto[b])
#      B[[b]]$L.lowerbound.vec[i] = loss_PS(B[[b]]$pareto.opt.X[i,] + max.X, loss.type, loss.params)  # here loss_PS should be vectorized 
    B[[b]]$L.lowerbound.vec = loss_PS_mat(B[[b]]$pareto.opt.X + t(replicate(n.pareto[b], max.X)), loss.type, loss.params)  # Vectorized version 
    cur.good.inds <- which(B[[b]]$L.lowerbound.vec <= new.L.upperbound)  # L.upperbound (old)
    
#    print(paste0("num. good inds: ", length(cur.good.inds), " out of: ", length(B[[b]]$L.lowerbound.vec)))
    new.c <- new.X <- c()
    ctr <- 0
#    print(paste0("Merge loop on good inds sub-block=", b, ", #good.inds=", length(cur.good.inds), " out of ", n.pareto[b]))
#    print("B-pareto-c")
#    print(B[[b]]$pareto.opt.c)
    if(is.vector(B[[b]]$pareto.opt.X) && (M.vec[b]>1))
      B[[b]]$pareto.opt.X <- as.matrix(t(B[[b]]$pareto.opt.X))
    if(is.vector(B[[b]]$pareto.opt.c) && (M.vec[b]>1))
      B[[b]]$pareto.opt.c <- as.matrix(t(B[[b]]$pareto.opt.c))
#    print("START LOOP")
#    print("B-pareto-c-now:")
#    print(B[[b]]$pareto.opt.c)
    
    # PROBLEM: VECTOR CAN BE FOR TWO REASONS:
    # 1. JUST ONE OPTIMAL VECTOR
    # 2. M=1
    # MUST DISTINGUISH THE TWO!!!
    for(j in cur.good.inds)  # loop over all vectors in the current stack (heavy loop!)
    {
#      print(paste0("j is:", j))
      ctr <- ctr + 1
      if(ctr%%1000 == 0)
        print(paste0("Run ind ", ctr, " of ", length(cur.good.inds)))
#      print(paste0("Comptue new v n.pareto=", n.pareto[b+1], " dim.b=", dim(B[[b]]$pareto.opt.X), 
#            " dim.b1=", dim(B[[b+1]]$pareto.opt.X), ", ", length(B[[b+1]]$loss.vec)))
#      print(paste("j=", j))
#      print("Bj")
#      print("B-pareto-x")
#      print(B[[b]]$pareto.opt.X)
      
      if(n.pareto[b]==1)
      {
        new.v <- matrix(rep(B[[b]]$pareto.opt.X, n.pareto[b+1]), nrow=n.pareto[b+1], byrow=TRUE) + B[[b+1]]$pareto.opt.X # B[[b]]$pareto.opt.X[j,] + B[[b+1]]$pareto.opt.X  # take all vectors together
      } else
      {
#        print(B[[b]]$pareto.opt.X[j,])
        new.v <- matrix(rep(B[[b]]$pareto.opt.X[j,], n.pareto[b+1]), nrow=n.pareto[b+1], byrow=TRUE) + B[[b+1]]$pareto.opt.X # B[[b]]$pareto.opt.X[j,] + B[[b+1]]$pareto.opt.X  # take all vectors together
      }      
      new.v <- list(pareto.X = new.v, pareto.inds = 1:n.pareto[b+1]) #      new.v <- get_pareto_optimal_vecs(new.v)

      new.L.lowerbound.vec = loss_PS_mat_rcpp(new.v$pareto.X + t(replicate(n.pareto[b+1], new.max.X)), loss.type, loss.params)    # Filter right away!!! 
      new.good.inds <- which(new.L.lowerbound.vec <= new.L.upperbound)  # L.upperbound (old)
      
#      print(paste0("Passed ", length(new.good.inds), " out of: ", n.pareto[b+1]))
      if(length(new.good.inds)>0)  # add vectors 
      {
        new.v$pareto.X = new.v$pareto.X[new.good.inds,]  # filter new ones !!! 
        new.v$pareto.inds = new.v$pareto.inds[new.good.inds]
            
#        print(paste0("MERGE NEW Passed ", length(new.good.inds), " out of: ", n.pareto[b+1]))
        # Next merge the two 
        new.X <- rbind(new.X, new.v$pareto.X)
#        print(paste0("MERGE C NEW Passed ", length(new.good.inds), " out of: ", n.pareto[b+1]))
#        print("INDS: ")
#        print(new.v$pareto.inds)
#        print("PARETO C:")
#        print(B[[b]]$pareto.opt.c)
#        print(B[[b+1]]$pareto.opt.c)
#        print("new.c:")
#        print(new.c)
        if(is.vector(B[[b]]$pareto.opt.c))  # Always keep a matrix when more than one pareto-optimal vector && (M.vec[b+1]>1))
          B[[b]]$pareto.opt.c <- t(matrix(B[[b]]$pareto.opt.c))
        if(is.vector(B[[b+1]]$pareto.opt.c))  # Always keep a matrix when more than one pareto-optimal vector && (M.vec[b+1]>1))
        {
          new.c <- rbind(new.c, cbind(matrix(rep(B[[b]]$pareto.opt.c[j,], length(new.v$pareto.inds)), nrow=length(new.v$pareto.inds), byrow=TRUE), 
                                      matrix(B[[b+1]]$pareto.opt.c, nrow= length(new.v$pareto.inds), byrow=FALSE)) ) # only one ! 
        } else
        {
          new.c <- rbind(new.c, cbind(matrix(rep(B[[b]]$pareto.opt.c[j,], length(new.v$pareto.inds)), nrow=length(new.v$pareto.inds), byrow=TRUE), 
                                  matrix(B[[b+1]]$pareto.opt.c[new.v$pareto.inds,], nrow= length(new.v$pareto.inds), byrow=FALSE)) )
        }
#        print(paste0("FINISHED MERGE C NEW Passed ", length(new.good.inds), " out of: ", n.pareto[b+1]))
      }    
      #      if( length(new.v$pareto.inds)<=1 )
#      {
#        print("1: ")
#        print(cbind(matrix(rep(B[[b]]$pareto.opt.c[j,], length(new.v$pareto.inds)), nrow=length(new.v$pareto.inds), byrow=TRUE)    ) )
#        print("2: ")
#        print( B[[b+1]]$pareto.opt.c[new.v$pareto.inds,] )
#      }
#      new.X <- get_pareto_optimal_vecs(new.X) # could be costly to run again all new.X against themselves 
#      if( length(new.X$pareto.inds)<=1 )
#        print(paste0("Num pareto union: ", length(new.X$pareto.inds)))
#      new.c <- new.c[new.X$pareto.inds,]
#      new.X <- new.X$pareto.X
    }  # loop on cur good inds 
  # New: Try uniting only at the end!!! 
#    print(paste0("Finished loop on good. inds. Total #vectors before Pareto: =", dim(new.c)[1]))
    new.X <- get_pareto_optimal_vecs(new.X) # could be costly to run again all new.X against themselves 
#    if( length(new.X$pareto.inds)<=1 )
#      print(paste0("Num pareto union: ", length(new.X$pareto.inds)))
    new.c <- new.c[new.X$pareto.inds,]
    new.X <- new.X$pareto.X
    
    if(is.vector(new.c))
    {
#      print("Only one!")
#      print(B[[b+1]]$pareto.opt.X)
      n.pareto[b+1] <- 1 # only one vector 
    } else
      n.pareto[b+1] <- dim(new.X)[1] # update number of vectors in next layer
    

#    print(paste0("Merge sub-block=", b, ", L=", n.pareto[b+1]))
    L.vec[(M.vec.cum[b+1]+1):(M.vec.cum[b+2]-1)] = L.vec[M.vec.cum[b+1]]  # update start 
#    print("SET L.vec now:")
#    print(L.vec)
    L.vec[M.vec.cum[b+2]] = n.pareto[b+1]  # update end 
#    print("SET2 L.vec now:")
#    print(L.vec)
    
    # update next layer: 
    B[[b+1]]$pareto.opt.X <- new.X
    B[[b+1]]$pareto.opt.c <- new.c
    if(n.pareto[b+1]==1)
    {
      B[[b+1]]$opt.X <- B[[b+1]]$max.X <- B[[b+1]]$pareto.opt.X
    } else
    {
      B[[b+1]]$max.X <- colMaxs(B[[b+1]]$pareto.opt.X, value = TRUE) # Update: Get maximum at each coordinate 
#      print(paste0("NOW Update OPT: Merge sub-block=", b, ", L=", n.pareto[b+1]))
      B[[b+1]]$opt.X <- B[[b+1]]$pareto.opt.X[which.min(loss_PS_mat(B[[b+1]]$pareto.opt.X, loss.type, loss.params)),] # Compute scores and update also the optimum 
#      print(paste0("Updated OPT: Merge sub-block=", b, ", L=", n.pareto[b+1]))
    }

    #    print(paste0("Num. Next layer: ", dim(new.X)[1]))
    #        cur.X <- new.X
    #    cur.c <- new.c
    #    L.vec[i] = dim(new.X)[1]  

#    print(paste0("Finished Merge sub-block=", b, ", L=", n.pareto[b+1]))
    
    now.bounds <- bound_monotone_loss_pareto_blocks_PS_mat(B[(b+1):loss.params$n.blocks], loss.type, loss.params)
#    print("BBB")
#    print(paste0("New-Bounds: [", now.bounds$lowerbound, ", ", now.bounds$upperbound, "], eps=", now.bounds$upperbound-now.bounds$lowerbound))    
    
        
  } # end loop on blocks 
  
#  print("merge!")
  merge.time <- difftime(Sys.time() , start.time, units="secs") - bb.time
#  print(paste0("merge time (sec.):", merge.time))
  
  
  # Finally find the cost-minimizer out of the Pareto-optimal vectors
  L <- n.pareto[b+1]
  if(is.null(L)) # one dimensional array 
  {
#    print("L=1 !!! ")
    L = 1
    loss.vec = loss_PS(new.X, loss.type, loss.params)
  }  else
  {
    #    loss.vec <- rep(0, L)
    #    for(i in 1:L)
    #      loss.vec[i] <- loss_PS(new.X[i,], loss.type, loss.params)
    loss.vec <- loss_PS_mat(new.X, loss.type, loss.params)
  }
  i.min <- which.min(loss.vec) # find vector minimizing loss 
  
  loss.eval.time <- difftime(Sys.time() , start.time, units="secs") - bb.time - merge.time 
#  print(paste0("loss eval time (sec.):", loss.eval.time))
  

  if(L == 1)  # output one vector at the end 
  {
#    print("RETURN L=1")
#    return(list(opt.X = new.X, opt.c = new.c, opt.loss = min(loss.vec)))
  }
  else
  {
    new.X <- new.X[i.min,]
    new.c <- new.c[i.min,]
  }
  return(list(opt.X = new.X, opt.c = new.c, opt.loss = min(loss.vec), L.vec = L.vec,
        bb.time = bb.time, merge.time = merge.time))
  
#  return(list(opt.X = new.X[i.min,], opt.c = new.c[i.min,], opt.loss = min(loss.vec), L.vec = L.vec,
#              bb.time = bb.time, merge.time = merge.time))
  
}  





###############################################################
# Filter solutions that can't be grown into optimal solution for monotone loss  
# Input: 
# sol - vector of solutions for partial problems (can be of different sizes)
# loss.type - what loss to compute 
# loss.params - parameters of loss 
#
# Output: 
# sol - new with reduced number of solutions
###############################################################
filter_solutions <- function(sol, loss.type, loss.params)
{
  if(!is_monotone_loss(loss.type))
  {
    print("Error! Can't filter solutions for non-monotone losses!!!")
    return(NULL)
  }
  n.blocks <- length(sol) # number of blocks 
  if(is.vector(sol[[1]]$pareto.opt.X))
    T <- length(sol[[1]]$pareto.opt.X)
  else
    T <- dim(sol[[1]]$pareto.opt.X)[2]

  opt.X <- max.X <- rep(0, T)
  n.pareto <- rep(0, n.blocks)
  for(b in 1:n.blocks) # loop on Multi-Blocks
  {
    n.pareto[b] <- length(sol[[b]]$loss.vec)  # number of vectors in each block
    max.X <- max.X + sol[[b]]$max.X # colMaxs(sol[[b]]$pareto.opt.X, value = TRUE) # Get minimum at each coordinate 
    opt.X <- opt.X + sol[[b]]$opt.X # Get best vector for each sub-problem (better than getting coordinate max) 
  }
  n.pareto.new <- n.pareto
  upper = loss_PS(opt.X, loss.type, loss.params)
  for(b in 1:n.blocks)
    if(n.pareto[b]>1)
    {
      cur.max.X <- max.X - sol[[b]]$max.X #   cur.opt.X <- opt.X - sol[[b]]$max.X 
#      print("INPUT TO LOSS:")
#      print(sol[[b]]$pareto.opt.X + t(replicate(n.pareto[b], cur.max.X)))
#      print("n.pareto")
      lower <- loss_PS_mat(sol[[b]]$pareto.opt.X + t(replicate(n.pareto[b], cur.max.X)), loss.type, loss.params)  # Vectorized version 
#      print("OUTPUT OF LOSS, lower:")
#      print(lower)
      cur.good.inds <- which(lower <= upper)  # Filter all vectors .. 
      n.pareto.new[b] <- length(cur.good.inds)
      
#      print("Filter: cur-good-inds:")
#      print(cur.good.inds)
#      print("Paret-X:")
#      print(sol[[b]]$pareto.opt.X)
      sol[[b]]$pareto.opt.X <- sol[[b]]$pareto.opt.X[cur.good.inds,]
      if(n.pareto.new[b]>1)
      {
        if(is.matrix(sol[[b]]$pareto.opt.c)) 
        {
          if(dim(sol[[b]]$pareto.opt.c)[2]==1)
            sol[[b]]$pareto.opt.c <- t(t(sol[[b]]$pareto.opt.c[cur.good.inds,]))
          else
            sol[[b]]$pareto.opt.c <- sol[[b]]$pareto.opt.c[cur.good.inds,]
        }
      }
      sol[[b]]$loss.vec <- sol[[b]]$loss.vec[cur.good.inds]
      if(n.pareto.new[b]==1)
        sol[[b]]$max.X  <- sol[[b]]$pareto.opt.X # update also max (one vector) 
      else
        sol[[b]]$max.X  <- colMaxs(sol[[b]]$pareto.opt.X, value = TRUE) # update also max 
    }
  
#  print("AFTER FILTERING NEW # VECTORS:")
#  print(n.pareto.new)
  return(list(sol=sol, n.pareto=n.pareto.new))  # updated array 
}  


###############################################################
# A closed-form solution for the relaxation version case of stabilizing selection 
#
# Input: 
# X - tensor of polygenic scores 
# loss.type - string signifying loss type
# loss.params - parameters of the loss function
#
# Output: 
# List with the following 
# opt.X - optimal X 
# opt.loss - loss of optimal X
# opt.c - optimal value of the loss 
# C.mat - matrix of selector variables 
# loss.mat -  
# Big.A - matrix of 
# b - 
###############################################################
optimize_C_stabilizing_exact <- function(X, loss.type, loss.params)
{
  if(!("eta" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$eta <- 0 
  
  print("Start optimize stabilizing")
  M <- dim(X)[1];   C <- dim(X)[2];  T <- dim(X)[3]
  #  A <- grad_loss_PS(X, C, "stabilizing", loss.params)
  
  A <- -2 * loss.params$eta * eye(M*C) # new: add regularization # matrix(0, nrow=M*C, ncol=M*C)
  for(k in c(1:T))
    A <- A + 2 * loss.params$theta[k] *  as.vector(t(X[,,k])) %*% t(as.vector(t(X[,,k])))
  E <- matrix(0, nrow=M, ncol=M*C)
  for(i in c(1:M))
    E[i,((i-1)*C+1):(i*C)] <- 1
  b <- c(rep(0, M*C), rep(1, M)) # free vector for linear system   
  Big.A <- rbind(cbind(A, t(E)), cbind(E, matrix(0, nrow=M, ncol=M)))
  
  #  return(list(  Big.A=Big.A, b=b)) # temp debug
  debug.print <- 0
  if(debug.print)
  {
    print("A:")
    print(A)
    print("Big.A:")
    print(Big.A)
    print("Dim(Big.A):")
    print(dim(Big.A))
    print("Rank(Big.A):")
    print(qr(Big.A)$rank)
  }  
  
  if(T+2*M >= (C+1)*M)  # unique solution 
    v <- solve(Big.A, b) # Solve system
  else  
    v <- pinv(Big.A) %*% b  # infinite solutions. Use pseudo-inverse
  C.mat <- matrix(v[1:(M*C)], nrow=M, ncol=C, byrow = TRUE)
  
#  v2 <- v + (eye(M*(C+1)) - pinv(Big.A)%*%Big.A) %*% rnorm(M*(C+1))
#  C.mat <- matrix(v2[1:(M*C)], nrow=M, ncol=C, byrow = TRUE)
  
  loss.mat <- loss_PS(compute_X_C_mat(X, C.mat), loss.type, loss.params) #  c.p.v <- compute_X_C_mat(X, C.mat)

#  print("TOTAL LOSS IS: ")
#  print(loss.mat - loss.params$eta * sum(C.mat**2))
#  c.stacked <- as.vector(C.mat)
#  c.stacked <- as.vector(t(C.mat))
#  print("TOTAL LOSS VEC IS:")
#  print(0.5 * t(c.stacked) %*% A %*% c.stacked)
#  print(0.5 * t(c.stacked) %*% A.clean %*% c.stacked)
#  
#  c.stacked.better <- as.vector((sol.rel.reg$C.mat))
#  print(0.5 * t(c.stacked.better) %*% A %*% c.stacked.better)
#  
#  c.stacked.better <- as.vector((t(sol.rel.reg$C.mat)))
#  A.clean <- A + 2 * loss.params$eta * eye(M*C) 
#  print(0.5 * t(c.stacked.better) %*% A.clean %*% c.stacked.better)
#  A.dirty <- -2 * loss.params$eta * eye(M*C)
#  print(0.5 * t(c.stacked.better) %*% A.dirty %*% c.stacked.better)
  
  
  # NEW: Do quadratic programming including inequality constraints (not closed-form)
#  if(quad.prog)
#  {
#    sol.qp <- quadprog(A.clean, rep(0, params$M*params$C), A = NULL, b = NULL,  # use pracma function
#             Aeq = E, beq = rep(1, params$M), lb = rep(0, params$M*params$C), ub = rep(1, params$M*params$C))  # WORKS ONLY FOR POS.DEF!!

#    sol.qp <- solveqp(A.clean, h=NULL, lb = rep(0, params$M*params$C), ub = rep(1, params$M*params$C), 
#                      A = E, Alb = rep(1, params$M), Aub = rep(1, params$M)) # equality constraint
#    A.clean <- diag( params$M*params$C)
#    sol.qp <- QPmin(G = A.clean, g = rep(0, params$M*params$C), A = t(E), b = rep(1, params$M), neq = params$M,
#                    Lb = rep(0, params$M*params$C), Ub = NULL)
#                       Lb = rep(0, params$M*params$C), Ub = rep(1, params$M*params$C), tol = 1e-06)  # WORKS ONLY FOR POS.DEF!!

#      }
  
  
#  c.vec <- max.col(C.mat) # Convert to matrix and take max of each row 
#  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.type, loss.params)  # cost of the rounded solution 
#  opt.X <- compute_X_c_vec(X, c.vec)

  ret <- real_to_integer_solution(X, C.mat, loss.type, loss.params)
  ret$C.mat <- C.mat
  ret$loss.mat <- loss.mat
  ret$Big.A <- Big.A
  ret$b <- ret$b

  print("Relax Loss:")
  print(ret$loss.mat)
  
    return(ret)
  
    
#  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, C.mat=C.mat, loss.mat=loss.mat, 
#              Big.A=Big.A, b=b))
}


###############################################################
# A closed-form solution for the relaxation version case of stabilizing selection 
#
# Input: 
# X - tensor of polygenic scores 
# loss.type - string signifying loss type
# loss.params - parameters of the loss function
#
# Output: 
# List with the following 
# opt.X - optimal X 
# opt.loss - loss of optimal X
# opt.c - optimal value of the loss 
# C.mat - matrix of selector variables 
# loss.mat -  
# Big.A - matrix of 
# b - 
###############################################################
optimize_C_stabilizing_SDR_exact <- function(X, loss.type, loss.params)
{
  if(!("eta" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$eta <- 0 
  
  print("Start optimize stabilizing Semidefinite Relaxation")
  M <- dim(X)[1];   C <- dim(X)[2];  T <- dim(X)[3]

  # Old stuff below: replace by SDR: 
  A <- -2 * loss.params$eta * eye(M*C) # new: add regularization # matrix(0, nrow=M*C, ncol=M*C)
  for(k in c(1:T))
    A <- A + 2 * loss.params$theta[k] *  as.vector(t(X[,,k])) %*% t(as.vector(t(X[,,k])))
  E <- matrix(0, nrow=M, ncol=M*C)
  for(i in c(1:M))
    E[i,((i-1)*C+1):(i*C)] <- 1
  print("Dims:")
  print(dim(A))
  print(dim( t(A %*% rep(1, M*C))     ))
  A.hom <- cbind( rbind(A, t(A %*% rep(1, M*C))), rbind(A %*% rep(1, M*C), 0) )
  print("Did A.hom")
  b <- c(rep(1, M*C+1), rep(2-C, M), (M*(2-C)+1)**2  ) # free vector for linear system   

  H <-  vector("list", length = M*C+M+2)  # list of pos-def matrices for the constraints 
  for(i in 1:(M*C+1))
  {
    H[[i]] <- list(matrix(0, nrow = M*C+1, ncol = M*C+1))
    H[[i]][[1]][i,i] = 1
  }
  for(i in 1:M)
  {
    H[[i+M*C+1]] <- list(matrix(0, nrow = M*C+1, ncol = M*C+1)) # set as zeros 
    H[[i+M*C+1]][[1]][M*C+1,1:(M*C)] <- 0.5*E[i,]  # factor 2 correlation
    H[[i+M*C+1]][[1]][1:(M*C),M*C+1] <- 0.5*t(E[i,])
  }
  H[[M*C+M+2]] <- list(matrix(1.0, nrow = M*C+1, ncol = M*C+1)) # set as zeros 
  
    print("Did H")
  
  K <- c()
  K$type = "s"  # positive semidefinite
  K$size = M*C+1
  SDR.ret <- csdp(list(-A.hom), H, b, K) # , control=csdp.control()) # package maximizes, need to take minus
  print("Finished SDP, now get integer solution")
  ret <- SDP_to_integer_solution(X, SDR.ret$X[[1]], loss.type, loss.params, method = loss.params$sdr_to_int)  # randomization")  # "svd"
  print("Got integer solution")
  
#  SDR.svd <- svd(SDR.ret$X[[1]], 1, 1) # take the best rank-1 approximation
#  heatmap.2(A.hom, scale = "none", col = bluered(100), dendrogram = "none",
#            trace = "none", density.info = "none", Rowv=NA, Colv=NA)
#  heatmap.2(SDR.ret$X[[1]], scale = "none", col = bluered(100), dendrogram = "none",
#            trace = "none", density.info = "none", Rowv=NA, Colv=NA)
#  correct.c.vec01 <- c(0,0,1,0,1,0,1)
#    correct.c.vec <-correct.c.vec01*2-1
#  correct.C.mat <- correct.c.vec %*% t(correct.c.vec)
#  correct.C.mat01 <- (correct.c.vec01 %*% t(correct.c.vec01))[1:(M*C),1:(M*C)]
#  rank1.approx <- SDR.svd$u %*% t(SDR.svd$u) * ((M*C+1)/sum(SDR.svd$u**2))  # normalize
#  heatmap.2(  correct.C.mat, scale = "none", col = bluered(100), dendrogram = "none",
#            trace = "none", density.info = "none", Rowv=NA, Colv=NA)
#
#  opt.C.relax.linear <- optimize_C_stabilizing_exact(X, loss.type, loss.params)
#  
#  writeLines(paste0("Opt SDR COST:", 0.25 * (sum(diag(A.hom %*% SDR.ret$X[[1]])) + sum(A)), 
#                    "\nOpt rank-1 approx COST:",    0.25 * (sum(diag(A.hom %*% rank1.approx)) + sum(A)), 
#                    "\nOpt linear-relaxation COST:",   opt.C.relax.linear$opt.loss, 
#                    "\nOpt orig problem COST:",  0.25 * ( sum(diag(A.hom %*%correct.C.mat)) + sum(A))  ))
#  sum(diag(A %*% correct.C.mat01)) - sum(A)

# constraint.vec.ret <- rep(0, length(H))
#  constraint.vec.correct <- rep(0, length(H))
#  for(i in 1:length(H))
#  {
#    constraint.vec.ret[i] <-sum(diag(H[[i]][[1]] %*% rank1.approx ))
#    constraint.vec.correct[i] <-sum(diag(H[[i]][[1]] %*% correct.C.mat ))
#  }    
#  # check trace of the solution: 
#  diag(SDR.ret$X[[1]])  # should be all ones 
#  sum(diag(H[[8]][[1]] %*% SDR.ret$X[[1]])) # should be C-2 = -1
#  sum(diag(H[[9]][[1]] %*% SDR.ret$X[[1]])) # should be C-2 = -1
  
#  ggplot(SDR.ret$X[[1]], aes(X, Y, fill= Z)) +     geom_tile()
##  heatmap(SDR.ret$X[[1]], Rowv=NA, Colv=NA)
##  plot(SDR.svd$d)
###  ret <- c()
#  ret$C.mat <- (sign(SDR.svd$u)+1)/2 # take first eigenvector 
  ###   ret$c.vec <- apply(matrix(SDR.svd$u[-1], nrow=M), 1, FUN = which.max) 
  ### ret$C.mat <- matrix(0, nrow=M, ncol=C)
  ### ret$C.mat[cbind(1:M, ret$c.vec)] <- 1
  # Specialized rounding: Take for each block of M the maximum value:
  
  ret$loss.mat <- loss_PS(compute_X_C_mat(X, ret$C.mat), loss.type, loss.params) #  c.p.v <- compute_X_C_mat(X, C.mat)
  ### print("SDP Loss:")
  ### print(ret$loss.mat)
  
  ### ret$opt.loss <- loss_PS(compute_X_c_vec(X, ret$c.vec), loss.type, loss.params)  # cost of the rounded solution 
  ### ret$opt.X <- compute_X_c_vec(X, ret$c.vec) #  opt.X <- compute_X_C_mat(X, C.cur)
  
  return(ret)  

#  ret <- real_to_integer_solution(X, C.mat, loss.type, loss.params)
#  ret$C.mat <- C.mat
#  ret$loss.mat <- loss.mat
#  ret$Big.A <- Big.A
#  ret$b <- ret$b
#  return(ret)
  
}



###############################################################
# A function for embryo selection algorithm  
# Input: 
# X - tensor of polygenic scores 
# loss.type - string signifying loss type
# loss.params - parameters of the loss function
#
# Output: 
# List with the following for the best embryo 
# opt.X - optimal X 
# opt.loss - 
# opt.c - optimal value of the loss 
###############################################################
optimize_C_embryo <- function(X, loss.type, loss.params)
{
  C = dim(X)[2]
  loss.vec <- rep(0, C)
  X.block.sum <- apply(X, 2, colSums) # sum over blocks 
  for(i in 1:C)
  {
    #    print("Compute i=")
    #    print(i)
    #    print("X vec col:")
    #    print(X.block.sum[,i])
    loss.vec[i] = loss_PS(X.block.sum[,i], loss.type, loss.params)  
  }
  
  opt.c <- which.min(loss.vec)
  #  print("Opt c=")
  #  print(opt.c)
  opt.loss <- min(loss.vec)
  #  print("Dim X block:")
  #  print(dim(X.block.sum))
  opt.X <- X.block.sum[,opt.c]
  
  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=opt.c)) # return identity of optimal embryo 
}


###############################################################
# A wrapper function for all optimization algorithms
# Find c vector minimizing: loss (\sum_i X[i,c.vec[i],])
#
# Input: 
# X - tensor of polygenic scores 
# loss.type - string signifying loss type
# loss.params - parameters of the loss function
# alg.str - which algorithm to use for optimization
#
# Output: 
# Whatever the algorithm called to returns
# opt.loss - 
# .c.opt - optimal value of the loss 
###############################################################
optimize_C <- function(X, loss.type, loss.params, alg.str)
{
  start.time <- Sys.time()
  M <- dim(X)[1]; C <- dim(X)[2];  T <- dim(X)[3]
  
  ret = switch(alg.str,
    "embryo" = optimize_C_embryo(X, loss.type, loss.params), 
    "quant" = optimize_C_quant(X, loss.type, loss.params),
    "closed_form" = optimize_C_stabilizing_exact(X, loss.type, loss.params),
    "relax" = optimize_C_relax(X, loss.params$C.init, loss.type, loss.params),
    "SDR" = optimize_C_SDR(X, loss.params$C.init, loss.type, loss.params),
    "SDR_closed_form" = optimize_C_stabilizing_SDR_exact(X, loss.type, loss.params),
    "branch_and_bound" = optimize_C_branch_and_bound(X, loss.type, loss.params),
    "branch_and_bound_lipschitz" = optimize_C_branch_and_bound_lipschitz(X, loss.type, loss.params),
    "branch_and_bound_divide_and_conquer" = optimize_C_branch_and_bound_divide_and_conquer(X, loss.type, loss.params)
  )
  ret$run.time <- difftime(Sys.time(), start.time, units = "secs")[[1]]  # add running time as output
  return(ret)
  #  if(alg.str == "embryo") # take best embryo (no separation to chromosomes)  
#    return(optimize_C_embryo(X, loss.type, loss.params))
#  if(loss.type == "quant") # easy optimization for quantitative traits 
#    return(optimize_C_quant(X, loss.type, loss.params))
#  if((alg.str == "closed_form") && (loss.type == "stabilizing"))
#    return(optimize_C_stabilizing_exact(X, loss.type, loss.params))
#  if(alg.str == "relax")  # here we need to set init
#    return(optimize_C_relax(X, loss.params$C.init, loss.type, loss.params))  
#  if(alg.str == "branch_and_bound")
#    return(optimize_C_branch_and_bound(X, loss.type, loss.params))
#  if(alg.str == "branch_and_bound_lipschitz")
#    return(optimize_C_branch_and_bound_lipschitz(X, loss.type, loss.params))
#  if(alg.str == "branch_and_bound_divide_and_conquer")
#    return(optimize_C_branch_and_bound_divide_and_conquer(X, loss.type, loss.params))
  
}  



###############################################################
# Compute average gain using simulations 
# Input: 
# params - dimensions (M, C and T)
# loss.type - cost function
# loss.params - cost function parameters 
#
# Output: 
# A list with the following items: 
# - rand.mat - matrix with values for random selection
# - gain.tensor - a 3rd-order table with the gain for each simulation
# - gain.mat - matrix of differences (algs * C)
###############################################################
compute_gain_sim <- function(params, loss.type, loss.params)
{
  if(!("do.checks" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$do.checks <- 0 
  
  n.algs <- length(params$alg.str)
  n.c <- length(params$c.vec)
  gain.mat <- matrix(rep(0, n.c*n.algs), nrow=n.c)
  gain.tensor <- array(0, dim=c(params$iters, n.c, n.algs))
  rand.mat <- matrix(rep(0, params$iters*n.c), nrow=params$iters)  # a matrix (one column for each c value) 
  runs.tensor <- vector("list", length = params$iters)
  
  for (t in 1:params$iters)
  {
    print(paste0("Compute gain, iter=", t, " out of ", params$iters))
    runs.tensor[[t]] = vector("list", length = n.c)
    # New: Set Covariance matrices inside function for each iteration ! 
    Sigma <- 0.5*diag(params$T) + matrix(0.5, nrow=params$T, ncol=params$T)   # Fixed trait-correlations matrix 
    df <- params$T # For wishart distribution
    Sigma.T <- rWishart(1, df, Sigma)[,,1]  # Randomize traits correlation matrix 
    params$max.C <- max(params$c.vec)
    Sigma.K <- 0.5*diag(params$max.C) + matrix(0.5, nrow=params$max.C, ncol=params$max.C)   # kinship-correlations matrix 
    
#    print(paste0(params$alg.str, ": Iter=", t, ", Dim: (M, C, T)=", params$M, " ", params$max.C, " ", params$T))
    X = simulate_PS_chrom_disease_risk(params$M, params$max.C, params$T, Sigma.T, Sigma.K, params$sigma.blocks, rep(0.5, k))

    # New: loop on all methods (Use same input X to reduce variance) 
    for(i.c in 1:n.c)  # loop on C, take partial data     
    {
      runs.tensor[[t]][[i.c]] = vector("list", length = n.algs)
      c <- params$c.vec[i.c]
      # Compute also score without selection: 
      rand.mat[t,i.c] <- loss_PS(compute_X_c_vec(X[,1:c,], rep(1, params$M)), loss.type, loss.params)
      for(a in 1:n.algs)
      {
        sol <- optimize_C(X[,1:c,], loss.type, loss.params, params$alg.str[a])
        save("X", "c", "loss.type", "loss.params", "params", "sol", file="temp_bad_loss.Rdata")
        if(loss.params$do.checks)
        {
          sol2 <- optimize_C(X[,1:2,], loss.type, loss.params, params$alg.str)
          runs.tensor[[t]][[i.c]][[a]] <- sol2
          if(sol$opt.loss > sol2$opt.loss)
            print("Error! Adding C increased error!")
          sol.e <- optimize_C(X, loss.type, loss.params, "embryo")
          if(sol$opt.loss > sol.e$opt.loss)
            print("Error! embryo selection has lower loss!")
          
          if(params$alg.str[a] != "embryo")
          {
            sol.bb <- optimize_C_branch_and_bound(X, loss.type, loss.params)
            if(abs(sol$opt.loss - sol.bb$opt.loss) > 0.000000001)
              print("Error! bb has different loss!")
            if(sol$opt.c != sol.bb$opt.c)
              print("Error! bb has different c!")
            if(max(abs(sol$opt.X - sol.bb$opt.X)) > 0.000000001)
              print("Error! bb has different X!")
          }
        }
        
        # Next compute average gain vs. random: 
#        print("dim: ")
#        print(dim(gain.tensor))
#        print(c(t, i.c, a))
#        print("Alg:")
#        print(params$alg.str[a])
#        print("opt loss:")
#        print(sol$opt.loss)
        gain.tensor[t,i.c,a] <- sol$opt.loss
        #    if("loss.mat" %in% names(sol))
        #      gain.mat[t] <- sol$loss.mat
        #    else
        #      gain.mat[t] <- 0 
      } # loop on algorithm 
    } # loop on C 
  } # loop on iters

  for(a in 1:n.algs)
    gain.mat[,a] <- t(colMeans(rand.mat) - colMeans(gain.tensor[,,a])) # compute mean random loss minus optimal loss.
  return(list(gain.mat=gain.mat, gain.tensor=gain.tensor, rand.mat=rand.mat, runs.tensor=runs.tensor)) # Need to reduce the mean gain without selection 
}



