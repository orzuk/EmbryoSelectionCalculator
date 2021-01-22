# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") # Or
source("chrom_select_funcs.R")
source("chrom_select_algs.R")

# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)
C <- 3 # number of chromosomal copies
T <- 4 # number of traits
M <- 12  # number of blocks 

df <- T # For wishart distribution
k <- 4

h.ps <- 0.3  # variane explained by the polygenic score 
prev <- logspace(-4, -0.5, n=T) # c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease . Should match T 
theta <- rep(1, T) # c(1, 1, 1, 1, 1)  # importance of each disease . Should match T 
sigma.blocks = chr.lengths * h.ps
Sigma.K <- 0.5*diag(C) + matrix(0.5, nrow=C, ncol=C)   # kinship-correlations matrix 
is.positive.definite(Sigma.K)
Sigma <- 0.5*diag(T) + matrix(0.5, nrow=T, ncol=T)   # trait-correlations matrix 
df = T
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 

X = simulate_PS_chrom_disease_risk(M, C, T, Sigma.T, Sigma.K, sigma.blocks[1:M], prev)


loss.C <- "disease"
loss.params <- c()
loss.params$K <- prev
loss.params$h.ps <- rep(h.ps, T)
loss.params$theta <- theta
loss.params$eta <- 0 # negative L2 regularization 
loss.params$n.blocks <- 4


sol.bb <- optimize_C_branch_and_bound(X, loss.C, loss.params)
sol.bb.mid <- optimize_C_branch_and_bound_lipschitz_middle(X, loss.C, loss.params)

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


V = t(array(c(c(7,0,2), c(1,2,3), c(4,5,4), c(2,4,1), c(-4,-1,6), c(3,2,1)), dim=c(3,5)))
V
V.pareto = get_pareto_optimal_vecs(V)
V.pareto

max_n <- 50
n.vec <- seq(5, max_n, 5)
compute.pareto.p.stat <- 0
if(compute.pareto.p.stat)
{
  n.vec <- C^(2:8)
  iters <- 50
  par <- compare_pareto_P(n.vec, k, C, iters)
  
  max.p.k.n <- max(max(par$p.k*n.vec), max(par$p.k.blocks*n.vec), max(par$p.k.asymptotic2*n.vec)) * 1.01
  
  plot(n.vec, par$p.k*n.vec, type="l", xlab="n", ylab="p_k(n) n", ylim = c(0, max.p.k.n))  # exact 
  lines(n.vec, par$p.k.asymptotic*n.vec, col="green")  # asymptotic 
  lines(n.vec, par$p.k.asymptotic2*n.vec, col="blue")  # asymptotic 
  lines(n.vec, par$p.k.sim*n.vec, col="red")  # simulation 
  lines(n.vec, par$p.k.blocks*n.vec, col="cyan")  # simulation 
  legend(0.8 * max(n.vec), 0.3*max.p.k.n,   lwd=c(2,2,2,2), 
         c("exact", "approx", "approx2", "sim", "block"), col=c("black", "green", "blue", "red", "cyan"), cex=0.75) #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  
  plot(n.vec, par$p.k.asymptotic / par$p.k, xlab="n", ylab="ratio")
  print(max(par$p.k.asymptotic / par$p.k))
}
#p_k <- rep(0, max_n)
#for(n in n.vec)
#{
#  p_k[n] <- pareto_P2(n, k)
#  # add also simulation  
#}
#plot(n.vec, p_k*n.vec, xlab="n", ylab="p_k(n) n")  # analytic 


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

  

# Test pareto-optimal function for k=2:
n <- 50
k <- 2
X.mat <- matrix(runif(n*k), nrow=n, ncol=k)
P <- get_pareto_optimal_vecs(X.mat)

plot(X.mat[,1], X.mat[,2])
points(P$pareto.X[,1], P$pareto.X[,2], pch=20, col="red", cex=2)


# Try worst case vectors. The goal is to find a set of vectors such that all C^K vectors will be pareto optimal 
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
