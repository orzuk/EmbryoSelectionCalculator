# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
#setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") # Or

root_dir = "C:/Code/GitHub/EmbryoSelectionCalculator/code/R/chrom" # Or. Change to your path
figs_dir = paste0(root_dir, "/Figures/")
sim_res_dir = paste(root_dir, "/sim_res/")

setwd(root_dir) 
source("chrom_select_funcs.R")
source("chrom_select_algs.R")
source("chrom_select_plots.R")

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)
C <- 10 # number of chromosomal copies
T <- 4 # number of traits
M <- 23  # number of blocks 

df <- T # For Wishart distribution

h.ps <- 0.3  # Variance explained by the polygenic score 
prev <- logspace(-4, -0.5, n=T) # c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease . Should match T 
theta <- rep(1, T) # c(1, 1, 1, 1, 1)  # importance of each disease . Should match T 
sigma.blocks = chr.lengths * h.ps
Sigma.K <- 0.5*diag(C) + matrix(0.5, nrow=C, ncol=C)   # kinship-correlations matrix 
is.positive.definite(Sigma.K)
Sigma <- 0.5*diag(T) + matrix(0.5, nrow=T, ncol=T)   # trait-correlations matrix 
df = T
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 

X = simulate_PS_chrom_disease_risk(M, C, T, Sigma.T, Sigma.K, sigma.blocks[1:M], prev)

# New: test pareto optimal 
X.mat <- X[1,,]

par.R <- get_pareto_optimal_vecs(X.mat)  
par.CPP <- get_pareto_optimal_vecs_rcpp(X.mat)  

if(!all(par.R$pareto.inds == par.CPP$pareto.inds+1))
  print("errr!!!! Pareto!!!")


loss.C <- "disease"
loss.params <- c()
loss.params$K <- prev
loss.params$h.ps <- rep(h.ps, T)
loss.params$theta <- theta
loss.params$eta <- 0 # negative L2 regularization 
loss.params$n.blocks <- 4


run.plots = TRUE
old.runs = FALSE
pareto.tests = FALSE



if(run.plots)  # New: plots for paper: 
{
  plot_BB_num_vectors_errorbars(params, time.iters = 10, save.figs = TRUE, force.rerun = FALSE)
    
  
  
  params$c.vec <- 2:5
  params$iters <- 5
  loss.params$n.blocks = 4 
  params$loss.type <- "disease"
  params$M <- 8  # reduce to run fast !!  
  
  
  params$c.vec <- 2:10
  params$iters <- 100
  loss.params$n.blocks = 13
  loss.params$eta <- 0.0
  loss.params$sdr_to_int <- "randomization" # randomization"  # "svd"
  params$loss.type <- "stabilizing"
  params$M <- 23  # reduce to run fast !! 
  
  plot_BB_accuracy(params, save.figs = TRUE, force.rerun = FALSE)
}



if(old.runs)
{
  sol.bb.mid <- optimize_C_branch_and_bound_divide_and_conquer(X, loss.C, loss.params)
  
  # Only if block size small enough
  sol.bb <- optimize_C_branch_and_bound(X, loss.C, loss.params)
  print(sol.bb.mid$opt.loss-sol.bb$opt.loss)
  #sol.bb.lip <- optimize_C_branch_and_bound_lipschitz(X, loss.C, loss.params)
  
  
  sol.quant <- optimize_C_quant(X, "quant", loss.params)
  
  # Example of choosing the index for each block
  c.vec = sample(C, M, replace=TRUE)
  C.mat = matrix(rexp(M*C), nrow=M, ncol=C)
  C.mat = C.mat / rowSums(C.mat)
  
  X.c = compute_X_c_vec(X, c.vec)
  X.c2 = compute_X_C_mat(X, C.mat)
  
  
  loss.C <- "stabilizing"
  bal.bb <- optimize_C_branch_and_bound(X, loss.C, loss.params)
  bal.relax <- optimize_C_relax(X, c(), loss.C, loss.params)
  
  
  loss_PS(compute_X_C_mat(X, C.mat), loss.C, loss.params)
  g.d = grad_loss_PS(X, C.mat, "disease", loss.params)
  g.b = grad_loss_PS(X, C.mat, "stabilizing", loss.params)
  
  H.d = hessian_loss_PS(X, C.mat, "disease", loss.params)
  H.b = hessian_loss_PS(X, C.mat, "stabilizing", loss.params)
  
  
  
  # Try many times: 
  diff.vec <- rep(0, 100)
  for(i in 1:100)
  {
    X = simulate_PS_chrom_disease_risk(M, C, T, Sigma.T, Sigma.K, sigma.blocks[1:M], prev)
    bal.exact <- optimize_C_stabilizing_exact(X, loss.C, loss.params)
    bal.exact2 <- optimize_C_stabilizing_exact(X[,1:2,], loss.C, loss.params) # take subset of C. Loss should be higher
    diff.vec[i] <- bal.exact$loss.mat - bal.exact2$loss.mat
    #  bal.exact$loss.mat
    #  bal.exact2$loss.mat
    if(bal.exact2$loss.mat < bal.exact$loss.mat)
      print("Error! loss got smaller!")
  }
  
  
  bal.relax2 <- optimize_C_relax(X, bal.exact$C.mat, loss.C, loss.params)
  plot(bal.relax2$loss.vec)
  
  grad_loss_PS(X, bal.exact$C.mat, loss.C, loss.params)
  
  C.mat <- bal.exact$C.mat
  #C.mat <- matrix(runif(M*C), nrow=M, ncol=C)
  #C.mat <- C.mat / rowSums(C.mat)
  X.c <- compute_X_C_mat(X, C.mat)
  LOSS1 <- loss_PS(X.c, loss.C, loss.params)
  C.mat.to.vec <- as.vector(t(C.mat))
  A <- bal.exact$Big.A[1:(M*C),1:(M*C)]
  LOSS2 <- 0.5 * t(C.mat.to.vec) %*% A %*% (C.mat.to.vec)
  print(LOSS1)
  print(LOSS2)
  
  
  average.disease.loss = sum(loss.params$theta * loss.params$K)
  print(paste("Optimal disease loss: ", sol$opt.loss, " Average disease loss: ", average.disease.loss))
}









if(pareto.tests)
{
  #po = 0
  #for(i in 1:730)
  #  po = po + is_pareto_optimal(sol$pareto.opt.X[i,], sol$pareto.opt.X)
  #print(po)
  
  #colMeans(sol$pareto.opt.X)
  
  sol$L.vec
  
  
  # Check projection: 
  C.mat <- V / sum(V)
  C.proj <- project_stochastic_matrix(C.mat)
  C.proj  
  
  
  
  # Test Pareto-optimal function for k=2:
  n <- 50
  k <- 2
  X.mat <- matrix(runif(n*k), nrow=n, ncol=k)
  P <- get_pareto_optimal_vecs(X.mat)
  
  plot(X.mat[,1], X.mat[,2])
  points(P$pareto.X[,1], P$pareto.X[,2], pch=20, col="red", cex=2)
  
  
  # Try worst case vectors. The goal is to find a set of vectors such that all C^K vectors will be Pareto optimal 
  M <- 6
  C <- 2 # number of vectors
  T <- 4 # dimension
  
  u <- matrix(runif(M*C), nrow=M, ncol=C)
  W = array(0, dim=c(M, C, T))
  for(i in 1:M)
    for(j in 1:C)
    {
      W[i,j,] = u[i,j]
      W[i,j,T] = (1-T)*u[i,j]
    }
  
  W[1,1,] <- c(0,1)
  W[1,2,] <- c(1,0)
  W[2,1,] <- c(0.5, sqrt(0.5))
  W[2,2,] <- c(sqrt(0.5), 0.5)
  sol.W <- optimize_C_branch_and_bound(W, "quant", loss.params)
  
  plot(sol.W$pareto.opt.X[,1], sol.W$pareto.opt.X[,2])
}
