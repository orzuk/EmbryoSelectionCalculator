library(cubature) # for multi-dimensional integration
library(iterpc)
library(Rmpfr) # arbitrary precision
library(pracma)
# library(VeryLargeIntegers)


# Function for counting pareto-optimal vectors 
# Compute Pareto optimal probability under independence with simulations 
pareto_P_sim <- function(n, k, iters=1000)
{
  n.pareto <- rep(0, iters)
  for(i in 1:iters)
    n.pareto[i] <- length(get_pareto_optimal_vecs(matrix(runif(n*k), nrow=n, ncol=k))$pareto.inds) # Simulate vectors 
  return(list(n.pareto=n.pareto, p.pareto=mean(n.pareto) / n))
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
  
  load("p_k_n_big_tab.Rdata")  # pre-computed (for speed)
  if((dim(p_k_n_big_tab)[1] >= n) & (dim(p_k_n_big_tab)[2] >= k))
    return(  p_k_n_big_tab[n,k] )
  
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

# Comptue the matrix of alpha_j,k coefficients and the vector of m_k coefficients
pareto_alpha_mat_m_vec <- function(max.k, log.flag = FALSE)
{
  log.m.k = c(0, (1:(max.k-1)) * ( log(1:(max.k-1)) - 1) - lfactorial(1:(max.k-1)))
  m.k <- exp(log.m.k)
  alpha.mat <- matrix(0, max.k, max.k)
  alpha.mat[1,1] <- 1
  for(k in c(2:max.k))
  {
    alpha.mat[2:k, k] <- alpha.mat[1:(k-1), k-1]
    alpha.mat[1,k] <- 2 * (sum(m.k[1:(k-1)] * alpha.mat[1:(k-1),k-1]) + m.k[k])
  }
  
  return(list(m.k=m.k, alpha.mat=alpha.mat))  
}


pareto_P_var_integrand_n2 <- function(x)
{
  
  k <- length(x) / 2 # no need to give as input
  x1 <- x[1:k]
  x2 <- x[(k+1):(2*k)]
  return(  (1 - all(x1 < x2) - all(x2 < x1))  )
}  


pareto_P_var_integrand_all <- function(x)
{
  k <- length(x) / 2 # no need to give as input
  x1 <- x[1:k]
  x2 <- x[(k+1):(2*k)]
  return(  all(x1 < x2)  )
}  


# a function of x a vector in R^2k, and k and n 
pareto_P_var_integrand <- function(x)
{
  
  #  print("Dim x start:")
  #  print(dim(x))
  #  print("X start:")
  #  print(x)
  k <- length(x) / 2 # no need to give as input
  x1 <- x[1:k]
  x2 <- x[(k+1):(2*k)]
  return(  (1 - all(x1 < x2) - all(x2 < x1)) * 
             (1 - prod(x1) - prod(x2) + prod(pmin(x1, x2)))^(my.n-2) )
  
  #  return(  (1 - all(x[1:k] < x[(k+1):(2*k)]) - all(x[(k+1):(2*k)] < x[1:k])) * 
  #             (1 - prod(x[1:k]) - prod(x[(k+1):(2*k)]) + prod(pmin(x[1:k], x[(k+1):(2*k)])))^(my.n-2) )
  
  
  #  print("Dim x:")
  #  print(dim(x))
  
  #  yy <- x1 < x2
  #  print(yy)
  #  pp <- apply(x1 < x2, 2, prod)
  #  print("pp:")
  #  print(pp)
  
  
  #  pm <- pmin(x1,x2)
  #  print("Dim pm:")
  #  print(dim(pm))
  #  r1 <- (1 - apply(x1 < x2, 2, prod) - apply(x2 < x1, 2, prod))
  #  print("r1:")
  #  print(r1)
  #  r2 <- (1 - apply(x1, 2, prod) - apply(x2, 2, prod)) #  + apply(pmin(x1,x2), 2, prod))
  #  print("r2:")
  #  print(r2)
  #  r <- (1 - apply(x1 < x2, 2, prod) - apply(x2 < x1, 2, prod)) * 
  #          (1 - apply(x1, 2, prod) - apply(x2, 2, prod) + apply(pmin(x1,x2), 2, prod))^(my.n-2)
  #  print("r:")
  #  print(r)
  #  print("Dim r:")
  #  print(dim(r))
  
  #  return(r)  
}

# Vectorized version 
pareto_P_var_integrand_v <- function(x)
{
  
  print("MAT!")
  k <- dim(x)[1] / 2 # no need to give as input
  x1 <- x[1:k,]
  x2 <- x[(k+1):(2*k),]
  return( matrix( (1 - apply(x1 < x2, 2, all) - apply(x2 < x1, 2, all)) * 
                    (1 - apply(x1, 2, prod) - apply(x2, 2, prod) + apply(pmin(x1,x2), 2, prod))^(my.n-2) , ncol(x)) )
}


my.sum <- function(x)
{
  #  print("dim x:")
  #  print(dim(x))
  #  print("x:")
  #  print(x)
  #  print("Dim output sum:")
  #  print(dim(sum(x)))
  if(is.vector(x))
    return(sum(x))
  else
    return(rowSums(x))
}

# Compute the variance of the Pareto-front size for n vectors in R^k
pareto_P_var <- function(k, n, integral_method = "cuhre", vectorize = FALSE) # "hcubature"
{
  start_time <- Sys.time()
  my.n <- n # give from outside
#  if(vectorize)
#  {
#    print("Vectorize")
#    V <- cubintegrate(f = pareto_P_var_integrand_v, lower = rep(0, 2*k), upper = rep(1, 2*k), method = integral_method, nVec = 128) # vectorInterface=TRUE "hcubature") # integrate   
#  }  else
#    V <- cubintegrate(f = pareto_P_var_integrand, lower = rep(0, 2*k), upper = rep(1, 2*k), method = integral_method) # "hcubature") # integrate   
  
  
  # New: Use combinatorial   
  print("Compute e_n_k")
  load("e_k_n_big_tab.Rdata")
  if(e_k_n_big_tab[k,n] >= 0)
    e_k_n <- e_k_n_big_tab[k,n]
  else
  {  
    e_k_n <- as.numeric(as.character(pareto_E_Z1Z2_python(as.integer(k), as.integer(n))))
    e_k_n_big_tab[k,n] <- e_k_n
    save(e_k_n_big_tab, file = "e_k_n_big_tab.Rdata") # update file 
  }
  print("Compute p_n_k")
  p_n_k <- pareto_P2(n, k)
  print(paste0("e_{n,k}=", round(e_k_n, 6), " , p_{n,k}=", round(p_n_k, 6), 
               " , p_{n,k}^2=", round(p_n_k^2, 6), " , COV_{n,k}=", round(e_k_n - p_n_k^2, 6)))
  
  return(list(V = n * p_n_k * (1-p_n_k) + nchoosek(n, 2) * (e_k_n - p_n_k^2), run.time = Sys.time() - start_time))
}


#V1 <- pareto_P_var(2, 3)
#V2 <- pareto_P_var(2, 3, "cuhre")
#cubintegrate(f = my.sum, lower = rep(0, 2*k), upper = rep(1, 2*k), method = "hcubature") # integrate   
#cubintegrate(f = pareto_P_var_integrand, lower = rep(0, 2*k), upper = rep(1, 2*k), method = "hcubature") # integrate   

# Probability that the two first vectors are in teh Pareto-front
pareto_E_Z1Z2 <- function(k, n, log.flag = FALSE, alternate.flag = FALSE, mulit.prec=FALSE, dig=128)
{
  start.time <- Sys.time()
  max.val <- 0
  e_k_n <- 0
  S1 <- 0
  S2 <- 0 
  S3 <- 0
  
  if(mulit.prec==TRUE)
  {
#    dig <- 128
    e_k_n <- mpfr(0, dig)
    kk <- mpfr(k, dig)
    nn <- mpfr(n, dig)
    run.vec <- mpfr(c(0:(n-2)), dig)
  }
  
  e_k_n_vec <- rep(0, (n-2)^3)
  ctr <- 1
  if(alternate.flag)  # sum (a+1,b,c,d) and (a,b,c+1,d) together
  {
    for(a in seq(0, n-3, 2)) # take only even values !!! 0:(n-3))  # don't take the last a !!!
      for(b in 0:(n-3-a))
        for(c in 0:(n-3-a-b))
        {
          d <- n-3-a-b-c  # (take a+1,b,c,d)
          
          log.frac <- log( (a+b+2*c+3)^k - (a+c+2)^k - (b+c+1)^k )  - k * (log(a+c+2) + log(b+c+1) + log(a+b+c+3))  # log-based implementation
          log.fact <- lfactorial(n-2) -sum(lfactorial(c(a+1,b,c,d)))
          log.ratio <- log(c+1) - log(a+1) + k*(log(b+c+2) - log(b+c+1) ) + log( (a+b+2*c+3)^k - (a+c+2)^k - (b+c+1)^k  ) - log( (a+b+2*c+4)^k - (a+c+2)^k - (b+c+2)^k  )
          
          
          #          print(paste0("Alt1: (a,b,c,d)= ", paste0(c(a+1,b,c,d), collapse=","), ", Val=", round( (-1)^(a+b+1) * exp(log.fact + log.frac), 3)))
          #          print(paste0("Alt2: (a,b,c,d)= ", paste0(c(a,b,c+1,d), collapse=","), ", Val=", round( (-1)^(a+b) * exp(log.fact + log.frac - log.ratio), 3)))
          
          #          print(paste0("log.ratio=", round(log.ratio, 5)))
          #          print(c(a,b,c,d,k))
          
          e_k_n <- e_k_n + (-1)^(a+b+1) * exp(log.fact + log.frac) * (1 - exp(-log.ratio))
          e_k_n_vec[ctr] <- (-1)^(a+b+1) * exp(log.fact + log.frac) * (1 - exp(-log.ratio))
          ctr <- ctr+1
        }
    
    # add boundaries  
    for(a in seq(0, n-2, 2)) # take only even values !!! 0:(n-3))  # don't take the last a !!!
      for(b in 0:(n-2-a))
      {
        #        a <- 0
        c <- 0
        d <- n-2-a-b-c  # (take a+1,b,c,d)
        log.frac <- log( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k )  - k * (log(a+c+1) + log(b+c+1) + log(a+b+c+2))  # log-based implementation
        log.fact <- lfactorial(n-2) -sum(lfactorial(c(a,b,c,d)))
        e_k_n <- e_k_n + (-1)^(a+b) * exp(log.fact + log.frac)
        e_k_n_vec[ctr] <- (-1)^(a+b) * exp(log.fact + log.frac)
        ctr <- ctr + 1
        #        print(paste0("Alt: (a,b,c,d)= ", paste0(c(a,b,c,d), collapse=","), ", Val=", round( (-1)^(a+b) * exp(log.fact + log.frac), 3)))
      }
    
    print("TOTAL VEC SUM:")
    print(sum(e_k_n_vec))
    
  } else
  {
    if(mulit.prec==TRUE)
    {
#      for(a in run.vec)
#        for(b in run.vec[1:(n-1-as.integer(a))]) 
#            for(c in run.vec[1:(n-1-as.integer(a+b))])
              
      for(a in 0:(n-2))
      {
        aa <- run.vec[a+1]
          for(b in 0:(n-2-a))
          {
            bb <- run.vec[b+1]
            for(c in 0:(n-2-a-b))              
            {
              cc <- run.vec[c+1]
              dd <- nn-2-aa-bb-cc
              log.frac <- log( (aa+bb+2*cc+2)^k - (aa+cc+1)^k - (bb+cc+1)^k )  - k * (log(aa+cc+1) + log(bb+cc+1) + log(aa+bb+cc+2))  # log-based implementation
              log.fact <- lfactorial(nn-2) -sum(lfactorial(c(aa,bb,cc,dd)))
              e_k_n <- e_k_n + (-1)^(a+b) * exp(log.fact + log.frac)

#                            d <- n-2-a-b-c
#              log.frac <- log( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k )  - k * (log(a+c+1) + log(b+c+1) + log(a+b+c+2))  # log-based implementation
#              log.fact <- lfactorial(nn-2) -sum(lfactorial(c(a,b,c,d)))
#              e_k_n <- e_k_n + (-1)^(a+b) * exp(log.fact + log.frac)
            }
          }
      }
    } else # don't use high-precision
    {
      for(a in 0:(n-2))
        for(b in 0:(n-2-a))
          for(c in 0:(n-2-a-b))
          {
            d <- n-2-a-b-c
            if(log.flag == FALSE)
              e_k_n <- e_k_n + (-1)^(a+b) * multichoose(c(a,b,c,d)) * 
                ( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k ) / ( (a+c+1)*(b+c+1)*(a+b+c+2) )^k  # naive implementation. Overflow for large n
            else
            {
              #            if(mulit.prec==TRUE)
              #            {
              #              aa <- mpfr(a, dig)
              #              bb <- mpfr(b, dig)
              #              cc <- mpfr(c, dig)
              #              dd <- mpfr(d, dig)
              #              log.frac <- log( (aa+bb+2*cc+2)^k - (aa+cc+1)^k - (bb+cc+1)^k )  - kk * (log(aa+cc+1) + log(bb+cc+1) + log(aa+bb+cc+2))  # log-based implementation
              #              log.fact <- lfactorial(nn-2) -sum(lfactorial(c(aa,bb,cc,dd)))
              #              e_k_n <- e_k_n + (-1)^(a+b) * exp(log.fact + log.frac)
              #            } else
              #            {

              # New: Split into three terms
              log.frac <- log( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k )  - k * (log(a+c+1) + log(b+c+1) + log(a+b+c+2))  # log-based implementation
              L1 <- k * (log(a+b+2*c+2) - log(a+c+1) - log(b+c+1) - log(a+b+c+2))  # log-based implementation
              L2 <- k * ( - log(a+c+1) - log(a+b+c+2))  # log-based implementation
              L3 <- k * ( - log(b+c+1) - log(a+b+c+2))  # log-based implementation
              
              log.fact <- lfactorial(n-2) -sum(lfactorial(c(a,b,c,d)))
              e_k_n_vec[ctr] <- (-1)^(a+b) * exp(log.fact + L2) # set S2
              ctr <- ctr + 1
              S1 <- S1 + (-1)^(a+b) * exp(log.fact + L1)
              S2 <- S2 + (-1)^(a+b) * exp(log.fact + L2)
              S3 <- S3 + (-1)^(a+b) * exp(log.fact + L3)
              
              
              e_k_n <- e_k_n + (-1)^(a+b) * exp(log.fact + log.frac)
              #              print(paste0("Add: (a,b,c,d)= ", paste0(c(a,b,c,d), collapse=","), ", Val=", round( (-1)^(a+b) * exp(log.fact + log.frac), 5)))
              #            }
              #            max.val <- max(max.val, exp(log.fact + log.frac))
            }
          } # end for loop
    } # end else high-precision
  }
  #  print("Max term:")
  #  print(max.val)
#  e_k_n_vec <- 0
  
  
  
  print(paste0("S1=", round(S1, 5), " ; S2=", round(S2, 5))) # , " ; S3=", round(S3, 5)))
  print(paste0("SUM=", round(S1-2*S2, 5), " ; e_k_n=", round(e_k_n, 5)))
  
  return(list(e_k_n=e_k_n, e_k_n_vec=e_k_n_vec[1:ctr], run.time = Sys.time()-start.time))
}        



# Get the lcm of 1,2,..,n
lcm_first_large_int <- function(n)
{
  first_primes <- as.numeric(VeryLargeIntegers::primes(n, test = "MR", iter = 10, bar=FALSE)) # get all primes 
  
  lcm <- as.vli(1)
  for(i in 1:length(first_primes))
  {
    lcm <- lcm * as.vli(first_primes[i]) ^ as.vli((log(n)%/%log(first_primes[i])))
  }
  return(lcm)
}


# Get float a/b for a,b big integers
divide_big_integers <- function(a, b)
{
  log.diff <- as.numeric(as.character(log10(a) - log10(b)))
  
  return( 10^log.diff * as.numeric(paste0("0.", as.character(a))) / as.numeric(paste0("0.", as.character(b))))
  
}



# Do computation with large integers package
pareto_E_Z1Z2_large_int <- function(k, n, dig=128)
{
  start.time <- Sys.time()
  max.val <- 0

  print("Start:")
  n_factorial_vec <- vli(n-1) # prepare in advance
  for(i in c(1:(n-1)))
    n_factorial_vec[[i]] <- factvli(max(1, i-1)) # 0!=1!=1
  max_denom_lcm = lcm_first_large_int(n)^(3*k)
  
  e_k_n_int <- as.vli(0)
  for(a in 0:(n-2))
  {
    if(a%%10 == 0)
      print(a)
    for(b in 0:(n-2-a))
    {
      for(c in 0:(n-2-a-b))              
      {
        d <- n-2-a-b-c
#        numerator <- (n_factorial_vec[[n-1]] / (n_factorial_vec[[a+1]] * n_factorial_vec[[b+1]] * n_factorial_vec[[c+1]] * n_factorial_vec[[d+1]])) * 
#          (-1)^(a+b) * as.vli( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k ) 
        numerator <- (n_factorial_vec[[n-1]] / (n_factorial_vec[[a+1]] * n_factorial_vec[[b+1]] * n_factorial_vec[[c+1]] * n_factorial_vec[[d+1]])) * 
           ( (a+b+2*c+2)^k - (a+c+1)^k - (b+c+1)^k ) 
        denominator <- as.vli(( (a+c+1)*(b+c+1)*(a+b+c+2) )^k)  # naive implementation. Overflow for large n
        
        if( (a+b)%%2 == 0)
          e_k_n_int <- e_k_n_int + numerator * (max_denom_lcm / denominator)
        else
          e_k_n_int <- e_k_n_int - numerator * (max_denom_lcm / denominator)
        
      }
    }
  }
  
  print(e_k_n_int)
  print(max_denom_lcm)
  
  e_k_n <- divide_big_integers(e_k_n_int, max_denom_lcm)
  
  return(list(e_k_n=e_k_n, run.time = Sys.time()-start.time))
}        


