# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") # Or
source("chrom_select.R")


# (sum(sqrt(chr.lengths)) + sum(sqrt(chr.lengths[1:22]))) / sqrt(2*pi)
C <- 2 # number of chromosomal copies
T <- 5 # number of traits
M <- 22 # number of blocks 

df <- 5 # For wishart distirbution
k <- 1
max_n <- 50
n_vec <- c(1:max_n)
p_k <- rep(0, max_n)
for(n in n_vec)
  p_k[n] <- pareto_P(n, k)

plot(n_vec, p_k*n_vec)

h.ps <- 0.3  # variane explained by the polygenic score 
prev <- c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease 

sigma.blocks = chr.lengths * h.ps
Sigma.K <- 0.5*diag(C) + matrix(0.5, nrow=C, ncol=C)   # kinship-correlations matrix 
is.positive.definite(Sigma.K)
Sigma <- 0.5*diag(T) + matrix(0.5, nrow=T, ncol=T)   # trait-correlations matrix 
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 
X = simulate_PS_chrom_disease_risk(M, C, T, Sigma.T, sigma.blocks, prev)

# Example of choosing the index for each block
c.vec = sample(C, M, replace=TRUE)
C.mat = matrix(rexp(M*C), nrow=M, ncol=C)
C.mat = C.mat / rowSums(C.mat)

X.c = compute_X_c_vec(X, c.vec)
X.c2 = compute_X_C_mat(X, C.mat)


V = t(array(c(c(7,0,2), c(1,2,3), c(4,5,4), c(2,4,1), c(-4,-1,6), c(3,2,1)), dim=c(3,5)))
V
V.pareto = get_pareto_optimal_vecs(V)
V.pareto

loss.C <- "disease"
loss.params <- c()
loss.params$K <- prev
loss.params$h.ps <- rep(h.ps, T)
x.opt <- optimize_C_branch_and_bound(X, loss.C, loss.params)


