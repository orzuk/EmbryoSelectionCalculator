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



# A branch and bound algorithm 
optimize_C_branch_and_bound <- function(X, loss.C, loss.params)
{
 print("Start optimize B&B") 
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  
  par.X <- get_pareto_optimal_vecs(X[1,,]) # Save only Pareto-optimal vectors . Needs fixing 
  cur.c <- t(t(par.X$pareto.inds))
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  L.vec <- rep(0, M)
  L.vec[1] = L
  print(paste("L=", L, " start loop"))
  for(i in 2:M)
  {  
    print(paste0("i=", i))
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1
    new.X <- c()
    new.c <- c()
    for(j in 1:L)  # loop over all vectors in the current stack      
      for(c in 1:C)  # loop over possible vectors to add 
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
optimize_C_branch_and_bound_lipschitz <- function(X, loss.type, loss.params)
{
  print("Start optimize B&B Lipschitz") 
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  lip <- get_tensor_lipshitz_params(X, loss.type, loss.params)  # get loss and bounds for individual vectors 
  
  lip$max.pos <- apply(lip$lip.pos.mat, 1, max)
  lip$max.neg <- apply(lip$lip.neg.mat, 1, max)
  # Need to save also c-vec for each branch

  par.X <- get_pareto_optimal_vecs(X[1,,]) # Save only Pareto-optimal vectors . Needs fixing 
  print("par x:")
  print(par.X)
  
  cur.c <- t(t(par.X$pareto.inds))
  cur.X <- par.X$pareto.X
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
    L = 1
  
  L.vec <- rep(0, M)
  L.vec[1] = L
  print(paste("L=", L, " start loop"))
  for(i in c(2:M))
  {  
#    print(paste0("i=", i))
    L <- dim(cur.X)[1]
    if(is.null(L)) # one dimensional array 
      L = 1

    # New: exclude vectors with high loss 
    if(i < M)
    {
#      print(paste0("L is:", L))
#      print("Dim cur.X:")
#      print(dim(cur.X))
      cur.scores <- rep(0, L)
      for(j in 1:L)
        cur.scores[j] = loss_PS(cur.X[j,], loss.type, loss.params) # New part: here filter pareto-optimal vectors with too low values 
      min.score <- min(cur.scores) # - sum(lip$max.pos[(i+1):M]) # Next, exclude     
      # Next, exclude all solutions with too high scores: 
#      print("Score bound:")
#      print(sum(lip$max.pos[(i+1):M]) +  sum(lip$max.neg[(i+1):M]))
#      print("Max cur score:")
#      print(max(cur.scores))
      
      
      good.inds <- which(cur.scores <= min.score + sum(lip$max.pos[(i+1):M]) +  sum(lip$max.neg[(i+1):M]))
      if(length(good.inds) < L)
      {
        print("Saved: Good inds:")
        print(length(good.inds))
        print(L)
      }
      cur.X <- cur.X[good.inds,] # take only good ones 
      cur.c <- cur.c[good.inds,]
    }  
    
    
    new.X <- c()
    new.c <- c()
    L <- length(good.inds)
#    print(dim(cur.X))
    for(j in 1:L)  # loop over all vectors in the current stack      
      for(c in 1:C)  # loop over possible vectors to add 
      {
        if(is.null(dim(cur.X)))
          v = cur.X+X[i,c,]
        else
          v = cur.X[j,]+X[i,c,]
        if(is_pareto_optimal(v, new.X))  # first check if pareto optimal 
        {
          new.X <- rbind(new.X, v)
          if(is.null(dim(cur.c)))
            new.c <- rbind(new.c, c(cur.c[j], c) )
          else
            new.c <- rbind(new.c, c(cur.c[j,], c) )
        }
      }
    L.vec[i] = dim(new.X)[1]  
    cur.X <- new.X
    cur.c <- new.c
    if(i == M)
      print(paste0("B&B C=", C, " i=", i, " out of ", M, " Stack Size:", dim(new.X)[1]))
  }
  
  # Finally find the cost-minimizer out ofthe Pareto-optimal vectors
  L <- dim(cur.X)[1]
  if(is.null(L)) # one dimensional array 
  {
    L = 1
    loss.vec = loss_PS(cur.X, loss.type, loss.params)
  }  else
  {
    loss.vec <- rep(0, L)
    for(i in 1:L)
      loss.vec[i] <- loss_PS(cur.X[i,], loss.type, loss.params)
  }
  i.min <- which.min(loss.vec) # find vector minimizing loss 
  
  return(list(opt.X = cur.X[i.min,], opt.c = cur.c[i.min,], opt.loss = min(loss.vec), 
              loss.vec = loss.vec, L.vec = L.vec, pareto.opt.X= cur.X, pareto.opt.c = cur.c))
}


# A branch and bound algorithm 
optimize_C_branch_and_bound_lipschitz_middle <- function(X, loss.type, loss.params)
{
  # Compute pareto optimal vectors for each half and then merge them 
  M1 <- floor(M/2)
  M2 <- M-M1
  
  X1 <- optimize_C_branch_and_bound(X[1:M1,,], loss.type, loss.params)
  X2 <- optimize_C_branch_and_bound(X[(M1+1):M,,], loss.type, loss.params)

  print("Lengths:")
  print(length(X1$loss.vec))
  print(length(X2$loss.vec))
  
  lip <- compute_lipschitz_const(loss.type, loss.params)
  
  lip1 <- c()
  lip2 <- c()
  lip1$pos.mat <- rep(0, length(X1$loss.vec))
  lip1$neg.mat <- rep(0, length(X1$loss.vec))
  lip2$pos.mat <- rep(0, length(X2$loss.vec))
  lip2$neg.mat <- rep(0, length(X2$loss.vec))
  for(i in 1:length(X1$loss.vec))
  {
    lip1$pos.mat[i] <-  pmax(X1$pareto.opt.X[i,], 0) %*% lip
    lip1$neg.mat[i] <- -pmin(X1$pareto.opt.X[i,], 0) %*% lip
  }
  for(i in 1:length(X2$loss.vec))
  {
    lip2$pos.mat[i] <-  pmax(X2$pareto.opt.X[i,], 0) %*% lip
    lip2$neg.mat[i] <- -pmin(X2$pareto.opt.X[i,], 0) %*% lip
  }
  print("Upperbounds:")
  L.lowerbound <- min(X1$opt.loss - max(lip2$pos.mat), X2$opt.loss - max(lip1$pos.mat))  
  L.upperbound <- max(X1$opt.loss + max(lip2$neg.mat), X2$opt.loss + max(lip1$neg.mat))  
  
  # Next exclude all vectors exceeding the upperbound
  good.inds1 <- which(X1$loss.vec - max(lip2$pos.mat) <= L.upperbound)
  good.inds2 <- which(X2$loss.vec - max(lip1$pos.mat) <= L.upperbound)
#  good.inds1 <- 1:length(X1$loss.vec) # TEMP DEBUG!
#  good.inds2 <- 1:length(X2$loss.vec) # TEMP DEBUG!

  # Try another upperbound: 
  X1$max.vec = rowMaxs(X1$pareto.opt.X) # take maximum value 
  X2$max.vec = rowMaxs(X2$pareto.opt.X) # take maximum value 
  X1$L.lowerbound.vec <- rep(0, length(X1$loss.vec))
  for(i in 1:length(X1$loss.vec))
    X1$L.lowerbound.vec[i] = loss_PS(X1$pareto.opt.X[i,] + X2$max.vec, loss.type, loss.params)
  good.inds1 <- intersect(good.inds1, which(X1$L.lowerbound.vec <= L.upperbound))
  
  min.loss <- L.upperbound
  
  print("good inds: ")
  print(paste0("1:", length(good.inds1), " of ", length(X1$loss.vec)))
  print(paste0("2:", length(good.inds2), " of ", length(X2$loss.vec)))
  print(good.inds1)
  print(good.inds2)
  
  new.good.inds1 <- good.inds1
  new.good.inds2 <- good.inds2
  
  print(paste0("Saved: ",  1-length(good.inds1)*length(good.inds2) / (length(X1$loss.vec)*length(X1$loss.vec))))
  for(i1 in good.inds1)
  {
    if(!(i1 %in% new.good.inds1)) # skip
    {
      print("SKIP i1!")
      next
    }
    
    print("i1:")
    print(i1)
    if(loss_PS(X1$pareto.opt.X[i1,] + X2$max.vec, loss.type, loss.params) > L.upperbound)
    {
      print("SKIP i1 with max!")
      next
    }
#    print("loop over i2:")
    
    for(i2 in good.inds2)
    {
#      print("i1, i2:")
#      print(i1)
#      print(i2)
      cur.loss <- loss_PS(X1$pareto.opt.X[i1,] + X2$pareto.opt.X[i2,], loss.type, loss.params)
      if(cur.loss < min.loss)
      {
        min.loss <- cur.loss
        opt.X <- X1$pareto.opt.X[i1,] + X2$pareto.opt.X[i2,]
        opt.c <- c(X1$pareto.opt.c[i1,], X2$pareto.opt.c[i2,])
        
        # update criteria for testing: 
        L.upperbound <- min.loss
        print("L Upper:")
        print(L.upperbound)
        print("L Upper + Lip:")
        print(L.upperbound + max(lip2$pos.mat) )
        print(L.upperbound + max(lip1$pos.mat) )
        
        new.good.inds1 <- which(X1$loss.vec - max(lip2$pos.mat) <= L.upperbound)
        new.good.inds2 <- which(X2$loss.vec - max(lip1$pos.mat) <= L.upperbound)
      }
      if(!(i2 %in% new.good.inds2)) # skip
      {
        print("SKIP i2!")
        next
      }
      
    }
  }
  print("Now return:")
  return(list(opt.X=opt.X, opt.c=opt.c, opt.loss = min.loss))      
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
optimize_C <- function(X, loss.C, loss.params, alg.str)
{
  M <- dim(X)[1]; C <- dim(X)[2];  T <- dim(X)[3]

  if(loss.C == "quant") # easy optimization for quantitative traits 
    return(optimize_C_quant(X, loss.C, loss.params))
  if(alg.str == "branch_and_bound")
    return(optimize_C_branch_and_bound(X, loss.C, loss.params))
  if(alg.str == "branch_and_bound_lipschitz")
    return(optimize_C_branch_and_bound_lipschitz(X, loss.C, loss.params))
  if(alg.str == "relax")  # here we need to set init
    return(optimize_C_relax(X, loss.params$C.init, loss.C, loss.params))  
  if(loss.C == "stabilizing")
    return(optimize_C_stabilizing_exact(X, loss.C, loss.params))
  
}  
  


# Compute average gain using simulations 
compute_gain_sim <- function(params, loss.C, loss.params)
{
  gain.vec <- rep(0, params$iters)
  gain.mat <- rep(0, params$iters)
  for (t in 1:params$iters)
  {
    X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
    sol <- optimize_C(X, loss.C, loss.params, params$alg.str)
    
    # Next compute average gain vs. random: 
    gain.vec[t] <- sol$opt.los
    if("loss.mat" %in% names(sol))
      gain.mat[t] <- sol$loss.mat
    else
      gain.mat[t] <- 0 
  }
  gain <- mean(gain.vec) # compute optimal loss. Should subtract mean loss  
  gain.mat <- mean(gain.mat)
  return(list(gain=gain, gain.mat=gain.mat))
}


