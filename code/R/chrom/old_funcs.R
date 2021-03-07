
# A branch and bound algorithm 
optimize_C_branch_and_bound_lipschitz_old <- function(X, loss.type, loss.params)
{
  print("Start optimize B&B Lipschitz") 
  M <- dim(X)[1]; C <- dim(X)[2]; T <- dim(X)[3]
  lip <- get_tensor_lipshitz_params(X, loss.type, loss.params)  # get loss and bounds for individual vectors 
  
  lip$max.pos <- apply(lip$lip.pos.mat, 1, max)
  lip$max.neg <- apply(lip$lip.neg.mat, 1, max)
  # Need to save also c-vec for each branch
  
  par.X <- get_pareto_optimal_vecs(X[1,,]) # Save only Pareto-optimal vectors . Needs fixing 
  
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


old_optimize_C_branch_and_bound_lipschitz_middle <- function(X, loss.type, loss.params)
{
  # Compute pareto optimal vectors for each half and then merge them 
  M1 <- floor(M/2)
  M2 <- M-M1
  
  X1 <- optimize_C_branch_and_bound(X[1:M1,,], loss.type, loss.params)
  X2 <- optimize_C_branch_and_bound(X[(M1+1):M,,], loss.type, loss.params)
  
  print("Lengths:")
  print(length(X1$loss.vec))
  print(length(X2$loss.vec))
  
  
  max.c <- rowMins(lip.tensor$X.loss.mat)
  max.X <- rep(0, T)
  for(i in 1:M)
    max.X <- max.X + X[i,max.c[i],]    
  L.upperbound <- loss_PS(max.X, loss.type, loss.params)
  print("Greedy sol: ")
  print(L.upperbound)
  
  
  lip1 <- c()
  lip2 <- c()
  lip1$pos.mat <- rep(0, length(X1$loss.vec))
  lip1$neg.mat <- rep(0, length(X1$loss.vec))
  lip2$pos.mat <- rep(0, length(X2$loss.vec))
  lip2$neg.mat <- rep(0, length(X2$loss.vec))
  for(i in 1:length(X1$loss.vec))
  {
    lip1$pos.mat[i] <-  pmax(X1$pareto.opt.X[i,], 0) %*% lip.alpha
    lip1$neg.mat[i] <- -pmin(X1$pareto.opt.X[i,], 0) %*% lip.alpha
  }
  for(i in 1:length(X2$loss.vec))
  {
    lip2$pos.mat[i] <-  pmax(X2$pareto.opt.X[i,], 0) %*% lip.alpha
    lip2$neg.mat[i] <- -pmin(X2$pareto.opt.X[i,], 0) %*% lip.alpha
  }
  print("Upperbounds:")
  L.lowerbound <- min(X1$opt.loss - max(lip2$pos.mat), X2$opt.loss - max(lip1$pos.mat))  
  L.upperbound <- max(X1$opt.loss + max(lip2$neg.mat), X2$opt.loss + max(lip1$neg.mat))  
  
  # greedy   
  L.upperbound <- loss_PS(X1$pareto.opt.X[which.min(X1$loss.vec),] + X2$pareto.opt.X[which.min(X2$loss.vec),], loss.type, loss.params)
  
  # Next exclude all vectors exceeding the upperbound
  good.inds1 <- which(X1$loss.vec - max(lip2$pos.mat) <= L.upperbound)
  good.inds2 <- which(X2$loss.vec - max(lip1$pos.mat) <= L.upperbound)
  #  good.inds1 <- 1:length(X1$loss.vec) # TEMP DEBUG!
  #  good.inds2 <- 1:length(X2$loss.vec) # TEMP DEBUG!
  
  # Try another upperbound: 
  X1$max.vec = colMaxs(X1$pareto.opt.X, value=TRUE) # take maximum value 
  X2$max.vec = colMaxs(X2$pareto.opt.X, value=TRUE) # take maximum value 
  X1$L.lowerbound.vec <- rep(0, length(X1$loss.vec))
  for(i in 1:length(X1$loss.vec))
    X1$L.lowerbound.vec[i] = loss_PS(X1$pareto.opt.X[i,] + X2$max.vec, loss.type, loss.params)
  good.inds1 <- intersect(good.inds1, which(X1$L.lowerbound.vec <= L.upperbound))
  
  min.loss <- L.upperbound + 0.0000000000001
  
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
    
    #    print("i1:")
    #    print(i1)
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
