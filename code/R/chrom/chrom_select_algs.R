# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(pracma)
library(tensor)
library(olpsR) # for projection onto the simplex remotes::install_github("ngloe/olpsR")

Rcpp::sourceCpp("cpp/chrom_funcs.cpp")  # fast functions  
source('chrom_select_funcs.R')



# a function for optimizing the selection of chromosomes:
# Run the optimization to find the optimal C:
optimize_C_relax <- function(X, C.init, loss.C, loss.params)
{
  print("Start optimize relax")
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
    C.cur <- project_stochastic_matrix(C.next) # Project onto the simplex 
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
    X = simulate_PS_chrom_disease_risk(M, C, k, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
    sol.bb <- optimize_C_branch_and_bound(X, "quant", loss.params) # run B&B. Loss at the end doesn't matter. 
    n.pareto <- n.pareto + length(sol.bb$loss.vec)
  }
  return(n.pareto / (iters*n))
}



# Choose in a greedy manner the best X. Mininize loss   
optimize_C_quant <- function(X, loss.C, loss.params)
{
  M <- dim(X)[1]
  T <- dim(X)[3]
  c.vec <- rep(0, M)
  opt.x <- rep(0, T)
  for( i in c(1:M))
  {
    v <- sum(X[i,,] * loss.params$theta)  
    c.vec[i] <- which.min(rowSums(X[i,,] * loss.params$theta)  )
    opt.x <- opt.x + X[i,c.vec[i],]
  }
  return(list(opt.X = opt.x, opt.c = c.vec, opt.loss = sum(opt.x*loss.params$theta)))
}



# A branch and bound algorithm for finding the X combination with minimal loss 
optimize_C_branch_and_bound <- function(X, loss.C, loss.params)
{
# print("Start optimize B&B") 
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  
  par.X <- get_pareto_optimal_vecs_rcpp(X[1,,]) # Save only Pareto-optimal vectors . Needs fixing 
#  par.X.R <- get_pareto_optimal_vecs(X[1,,])
  cur.c <- t(t(par.X$pareto.inds))
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  L.vec <- rep(0, M)
  L.vec[1] = L
#  print(paste("L=", L, " start loop"))
  for(i in 2:M) # loop on block
  {
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    print(paste0("B&B i=", i, " L=", L))
    new.X <- c()
    new.c <- c()
    
    # We know that the first vectors pareto optimal
    if(L>1)
    {
      new.X <- t(t(cur.X) + X[i,1,])  # may need to transpose here !
    } else
      new.X <- matrix(cur.X + X[i,1,], nrow=1)
    new.c <- cbind(cur.c, rep(1, L) )

    # new version: create sums and take union
    for(c in 2:C)
    {
      temp.X <- t(t(cur.X) + X[i,c,])
      temp.X <- get_pareto_optimal_vecs_rcpp(temp.X)
#      print("new X:")
#      print(new.X)
#      print(dim(new.X))
      union.X <- union_pareto_optimal_vecs_rcpp(new.X, temp.X$pareto.X)
      new.X <- union.X$pareto.X
      add.c <- cbind(cur.c, rep(c, L) )
      new.c <- rbind(new.c[union.X$pareto.inds1,], add.c[union.X$pareto.inds2,]) # need to modify here indices 
    }
    cur.X <- new.X
    cur.c <- new.c
    if(is.null(dim(new.X)[1]))
      L.vec[i] = 1
    else
      L.vec[i] = dim(new.X)[1]  
        
    #    if(i == M)
#      print(paste0("B&B C=", C, " i=", i, " out of ", M, " Stack Size:", dim(new.X)[1]))
  } # end loop on blocks 
  
  # Finally find the cost-minimizer out of the Pareto-optimal vectors
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
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), 
              loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X, pareto.opt.c = cur.c))
}



# A branch and bound algorithm 
optimize_C_branch_and_bound_lipschitz_middle <- function(X, loss.type, loss.params)
{
  # Add timing: 
  start.time <- Sys.time()
  M <- dim(X)[1];   C <- dim(X)[2];  T <- dim(X)[3]
  if(!("n.blocks" %in% names(loss.params)))  
    loss.params$n.blocks <- 2  # default: divide to two blocks 
  lip.alpha <- compute_lipschitz_const(loss.type, loss.params)
  lip.tensor <- get_tensor_lipshitz_params(X, loss.type, loss.params)  # get loss and bounds for individual vectors 

  # Divide to blocks: 
  M.vec <- rep(floor(M/loss.params$n.blocks), loss.params$n.blocks)
  if( mod(M, loss.params$n.blocks)>0 )
    M.vec[1:mod(M, loss.params$n.blocks)] <- M.vec[1:mod(M, loss.params$n.blocks)] + 1 # correct to sum to M
  M.vec.cum <- c(0, cumsum(M.vec))
  B <- vector("list", loss.params$n.blocks) 
  n.pareto <- rep(0,  loss.params$n.blocks)
  opt.X.upperbound <- rep(0, T)
  opt.c.upperbound <- c()
  for(b in 1:loss.params$n.blocks) # get pareto-optimal vectors for each block 
  {
    B[[b]] <- optimize_C_branch_and_bound(X[(M.vec.cum[b]+1):M.vec.cum[b+1],,], loss.type, loss.params)  # compute pareto optimal vectors for block
    opt.X.upperbound <- opt.X.upperbound + B[[b]]$opt.X
    opt.c.upperbound <- c(opt.c.upperbound,  B[[b]]$opt.c)
    B[[b]]$max.X <- colMaxs(B[[b]]$pareto.opt.X, value = TRUE) # Get maximum at each coordinate 
    n.pareto[b] <-  length(B[[b]]$loss.vec)
#    print("BB Pareto X:")
#    print(B[[b]]$pareto.opt.X)
#    print("BB Pareto c:")
#    print(B[[b]]$pareto.opt.c)
  }
  L.upperbound <- loss_PS(opt.X.upperbound, loss.type, loss.params) + 0.00000000001
  bb.time <- difftime(Sys.time() , start.time, units="secs") 
  
  print(paste0("b&b time (sec.): ", bb.time))

    
  ##############################################
  # Next loop from one side and merge blocks:
  ##############################################
  if(loss.params$n.blocks == 1)
  {
    run.blocks <- c()
    new.X <- B[[b]]$pareto.opt.X
    new.c <- B[[b]]$pareto.opt.c
  } else
    run.blocks <- 1:(loss.params$n.blocks-1)
  
  for(b in run.blocks)
  {
    max.X <- rep(0, T)
    for(b2 in c((b+1):loss.params$n.blocks))
      max.X <- max.X + B[[b2]]$max.X # this reduction can be applied on both sides 
    B[[b]]$L.lowerbound.vec <- rep(0, n.pareto[b])
    for(i in 1:n.pareto[b])
      B[[b]]$L.lowerbound.vec[i] = loss_PS(B[[b]]$pareto.opt.X[i,] + max.X, loss.type, loss.params)  # here loss_PS should be vectorized 
    cur.good.inds <- which(B[[b]]$L.lowerbound.vec <= L.upperbound)
    
    if(b == run.blocks) # last layer
    {
      min.loss <- 999999999999
      min.b <- 0
      for(j in cur.good.inds)  # loop over all vectors in the current stack   
      {
        new.v <- matrix(rep(B[[b]]$pareto.opt.X[j,], n.pareto[b+1]), nrow=n.pareto[b+1], byrow=TRUE) + B[[b+1]]$pareto.opt.X # B[[b]]$pareto.opt.X[j,] + B[[b+1]]$pareto.opt.X  # take all vectors together
        new.loss.vec <- loss_PS_mat(new.v, loss.type, loss.params)
        i.min <- which.min(new.loss.vec)
        new.min.loss <- min(new.loss.vec)
        if(new.min.loss < min.loss)
        {
          min.loss <- new.min.loss
          min.X <- new.v[i.min,]
          print(paste("j=", j))
          print("dim new.c:")
          print(dim(new.c))
          print(paste0("i.min=", i.min))
          print("Dim b+1 pareto:")
          print(dim(B[[b+1]]$pareto.opt.c))
          
          min.c <- c(new.c[j,],  B[[b+1]]$pareto.opt.c[i.min,])  # need to modify here 
        }
      }
      
      return(list(opt.X = min.X, opt.c =min.c, opt.loss = min.loss, bb.time = bb.time, merge.time = merge.time))
    }
    
    
#    print(paste0("num. good inds: ", length(cur.good.inds), " out of: ", length(B[[b]]$L.lowerbound.vec)))
    new.X <- c()
    new.c <- c()
    ctr <- 0
    for(j in cur.good.inds)  # loop over all vectors in the current stack   
    {
      ctr <- ctr + 1
      new.v <- matrix(rep(B[[b]]$pareto.opt.X[j,], n.pareto[b+1]), nrow=n.pareto[b+1], byrow=TRUE) + B[[b+1]]$pareto.opt.X # B[[b]]$pareto.opt.X[j,] + B[[b+1]]$pareto.opt.X  # take all vectors together
      new.v <- get_pareto_optimal_vecs_rcpp(new.v)

      # Next merge the two 
      new.X <- rbind(new.X, new.v$pareto.X)
      if( length(new.v$pareto.inds)<=1 )
      {
        print("1: ")
        print(cbind(matrix(rep(B[[b]]$pareto.opt.c[j,], length(new.v$pareto.inds)), nrow=length(new.v$pareto.inds), byrow=TRUE)    ) )
        print("2: ")
        print( B[[b+1]]$pareto.opt.c[new.v$pareto.inds,] )
      }
      new.c <- rbind(new.c, cbind(matrix(rep(B[[b]]$pareto.opt.c[j,], length(new.v$pareto.inds)), nrow=length(new.v$pareto.inds), byrow=TRUE), 
                                  matrix(B[[b+1]]$pareto.opt.c[new.v$pareto.inds,], nrow= length(new.v$pareto.inds), byrow=TRUE)) )
      new.X <- get_pareto_optimal_vecs_rcpp(new.X) # could be costly to run again all new.X against themselves 
      if( length(new.X$pareto.inds)<=1 )
        print(paste0("Num pareto union: ", length(new.X$pareto.inds)))
      new.c <- new.c[new.X$pareto.inds,]
      new.X <- new.X$pareto.X
    }  # loop on cur good inds 
    print(paste0("Merge i=", b, "L=", dim(new.c)[1]))
  # update next layer: 
    B[[b+1]]$pareto.opt.X <- new.X
    B[[b+1]]$pareto.opt.c <- new.c
    B[[b+1]]$max.X <- colMaxs(B[[b+1]]$pareto.opt.X, value = TRUE) # Update: Get maximum at each coordinate 
    if(is.matrix(B[[b+1]]$pareto.opt.X))
      n.pareto[b+1] <- dim(B[[b+1]]$pareto.opt.X)[1] # update number of vectors in next layer
    else
    {
      print("Only one!")
      print(B[[b+1]]$pareto.opt.X)
      n.pareto[b+1] <- 1 # only one vector 
    }
    
#    print(paste0("Num. Next layer: ", dim(new.X)[1]))
    
#        cur.X <- new.X
#    cur.c <- new.c
#    L.vec[i] = dim(new.X)[1]  
    if(b+1 == loss.params$n.blocks)
      print(paste0("Finished B&B Mid. Stack Size:", dim(new.X)[1]))
  } # end loop on blocks 
  
  print("merge!")
  merge.time <- difftime(Sys.time() , start.time, units="secs") - bb.time
  print(paste0("merge time (sec.):", merge.time))


  
  # Finally find the cost-minimizer out of the Pareto-optimal vectors
  L <- dim(new.X)[1]
  if(is.null(L)) # one dimensional array 
  {
    print("L=1 !!! ")
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
  print(paste0("loss eval time (sec.):", loss.eval.time))
  
  

#  print("loss.vec:")
#  print(loss.vec)
#  print("new c:")
#  print(new.c)
#  print("i min:")
#  print(i.min)
  
  if(L == 1)
  {
    print("i.min:")
    print(i.min)
    return(list(opt.X = new.X, opt.c = new.c, opt.loss = min(loss.vec)))
  }
  else
    return(list(opt.X = new.X[i.min,], opt.c = new.c[i.min,], opt.loss = min(loss.vec), bb.time = bb.time, merge.time = merge.time))
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
  if(!("eta" %in% names(loss.params)))   # negative L2 regularizer
    loss.params$eta <- 0 
  
  print("Start optimize stabilizing")
  M <- dim(X)[1];   C <- dim(X)[2];  T <- dim(X)[3]
  #  A <- grad_loss_PS(X, C, "stabilizing", loss.params)
  
  A <- -loss.params$eta * eye(M*C) # new: add regularization # matrix(0, nrow=M*C, ncol=M*C)
  for(k in c(1:T))
    #    A <- A + 2 * loss.params$theta[k] *  as.vector((X[,,k])) %*% t(as.vector((X[,,k])))
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

  c.p.v <- compute_X_C_mat(X, C.mat)
  loss.mat <- loss_PS(compute_X_C_mat(X, C.mat), loss.C, loss.params)
  c.vec <- max.col(C.mat) # Convert to matrix and take max of each row 

  opt.loss <- loss_PS(compute_X_c_vec(X, c.vec), loss.C, loss.params)  # cost of the rounded solution 
  opt.X <- compute_X_c_vec(X, c.vec)

  return(list(opt.X=opt.X, opt.loss=opt.loss, opt.c=c.vec, C.mat=C.mat, loss.mat=loss.mat, 
              Big.A=Big.A, b=b))
}


# A wrapper function for all optimizations
optimize_C_embryo <- function(X, loss.C, loss.params)
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
    loss.vec[i] = loss_PS(X.block.sum[,i], loss.C, loss.params)  
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
  
  


# A wrapper function for all optimizations
optimize_C <- function(X, loss.C, loss.params, alg.str)
{
  M <- dim(X)[1]; C <- dim(X)[2];  T <- dim(X)[3]

  if(alg.str == "embryo") # take best embryo (no separation to chromosomes)  
    return(optimize_C_embryo(X, loss.C, loss.params))
  
  
  if(loss.C == "quant") # easy optimization for quantitative traits 
    return(optimize_C_quant(X, loss.C, loss.params))
  if(alg.str == "relax")  # here we need to set init
    return(optimize_C_relax(X, loss.params$C.init, loss.C, loss.params))  
  if(loss.C == "stabilizing")
    return(optimize_C_stabilizing_exact(X, loss.C, loss.params))
  if(alg.str == "branch_and_bound")
    return(optimize_C_branch_and_bound(X, loss.C, loss.params))
  if(alg.str == "branch_and_bound_lipschitz")
    return(optimize_C_branch_and_bound_lipschitz(X, loss.C, loss.params))
  if(alg.str == "branch_and_bound_lipschitz_middle")
    return(optimize_C_branch_and_bound_lipschitz_middle(X, loss.C, loss.params))
}  
  


# Compute average gain using simulations 
compute_gain_sim <- function(params, loss.C, loss.params)
{
  n.algs <- length(params$alg.str)
  gain.vec <- rep(0, params$iters)
  gain.mat <- matrix(rep(0, params$iters*n.algs), nrow=params$iters)
  rand.vec <- rep(0, params$iters)
  for (t in 1:params$iters)
  {
    print(paste0(params$alg.str, ": Iter=", t, ", Dim: (M, C, T)=", params$M, " ", params$C, " ", params$T))
    X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))

    #    print("Solve:")
    
    # New: loop on all methods (same X to reduce variance) 
    
    # Compute also score without selection: 
    rand.vec[t] <- loss_PS(compute_X_c_vec(X, rep(1, params$M)), loss.C, loss.params)
    for(a in 1:n.algs)
    {
      sol <- optimize_C(X, loss.C, loss.params, params$alg.str[a])
#    sol2 <- optimize_C(X[,1:2,], loss.C, loss.params, params$alg.str)
#    if(sol$opt.loss > sol2$opt.loss)
#      print("Error! Adding C increased error!")
#    sol.e <- optimize_C(X, loss.C, loss.params, "embryo")
#    if(sol$opt.loss > sol.e$opt.loss)
#      print("Error! embryo selection has lower loss!")
#    

    # Next compute average gain vs. random: 
      gain.mat[t,a] <- sol$opt.loss
#    if("loss.mat" %in% names(sol))
#      gain.mat[t] <- sol$loss.mat
#    else
#      gain.mat[t] <- 0 
    } # loop on algorithm 
  } # loop on iters
  gain <- colMeans(gain.mat) - mean(rand.vec) # compute optimal loss. Should subtract mean loss  
#  gain.mat <- mean(gain.mat) # gain.mat=gain.mat,
  return(list(gain=gain,  gain.mat=gain.mat, rand.vec=rand.vec)) # Need to reduce the mean gain without selection 
}


