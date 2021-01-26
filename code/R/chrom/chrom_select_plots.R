# Plot performance
# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") # Or
source("chrom_select_funcs.R")
source("chrom_select_algs.R")


# SEt all parameters
params <- c()
params$M <- 6 # 22 
params$c.vec <- 2:5
params$T <- 5
params$iters <- 100
df <- 5 # For wishart distribution

h.ps <- 0.3  # variane explained by the polygenic score 

sigma.blocks = chr.lengths * h.ps
Sigma <- 0.5*diag(params$T) + matrix(0.5, nrow=params$T, ncol=params$T)   # trait-correlations matrix 
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 
# Loss parameters 
loss.C <- "disease" # stabilizing" # "disease"
loss.params <- c()
loss.params$K <- c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease 
loss.params$h.ps <- rep(h.ps, params$T)
loss.params$theta <- c(1, 1, 1, 1, 1)  # importance of each disease 
loss.params$eta <- 0 # negative L2 regularization 


gain.vec <- rep(0, length(params$C.vec))
gain.mat <- rep(0, length(params$C.vec))
run.plots <- 1
params$alg.str <- "branch_and_bound_lipschitz_middle" # "exact" # "branch_and_bound"
if(run.plots)
  for(i in 1:length(params$c.vec))
  {
    params$C <- params$c.vec[i]
    Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
    is.positive.definite(Sigma.K)
    L  <- compute_gain_sim(params, loss.C, loss.params)
    gain.vec[i] <- L$gain
#    gain.mat[i] <- L$gain.mat
  }
    
# Plot: 
plot(params$c.vec, gain.vec, xlab="C", ylab="Gain", ylim = c(min(0, min(gain.vec)), max(0, max(gain.vec))), main=paste0("Gain for ", loss.C, " loss"))
# points(params$c.vec, gain.mat, col="red", xlab="C", ylab="Gain Relaxed")



    