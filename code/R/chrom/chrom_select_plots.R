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
params$M <- 22
params$C.vec <- 1:10
params$T <- 5
params$iters <- 10
df <- 5 # For wishart distribution

h.ps <- 0.3  # variane explained by the polygenic score 
prev <- c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease 
theta <- c(1, 1, 1, 1, 1)  # importance of each disease 
sigma.blocks = chr.lengths * h.ps
Sigma <- 0.5*diag(params$T) + matrix(0.5, nrow=params$T, ncol=params$T)   # trait-correlations matrix 
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 



gain.vec <- rep(0, length(params$C.vec))
run.plots <- 1
if(run.plots)
  for(i in 1:length(params$C.vec))
  {
    params$C <- params$C.vec[i]
    Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
    is.positive.definite(Sigma.K)
    gain.vec[i] <- compute_gain_sim(params, loss.C, loss.params)
    
  }
    
    