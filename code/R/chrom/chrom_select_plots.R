# Plot performance
# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") 
source("chrom_select_funcs.R")
source("chrom_select_algs.R")


figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\"


# SEt all parameters
params <- c()
params$M <- 12 # try full chromosomes  
params$c.vec <- 2:10
params$T <- 5 # number of traits 
params$iters <- 3
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
loss.params$n.blocks <- 2

gain.embryo.vec <- bb.gain.vec <- gain.vec <- rep(0, length(params$c.vec))  # the gain when selecting embryos (no chromosomes)
run.plots <- 1
params$alg.str <- c("embryo", "branch_and_bound_lipschitz_middle") # ) "branch_and_bound") # "exact" # "branch_and_bound"
n.methods <- length(params$alg.str)
gain.mat <- matrix(rep(0, length(params$c.vec)*n.methods), ncol = n.methods)
#embryo.loss.params = loss.params
#embryo.loss.params$alg.str = "embryo"
#bb.loss.params = loss.params
#bb.loss.params$alg.str = "branch_and_bound"
if(run.plots)
  for(i in 1:length(params$c.vec))
  {
    params$C <- params$c.vec[i]
    Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
    is.positive.definite(Sigma.K)
    gain.mat[i,]  <- compute_gain_sim(params, loss.C, loss.params)$gain # chromosomal selection
#    bb.gain.vec[i]  <- compute_gain_sim(params, loss.C, bb.loss.params)$gain # chromosomal selection
#    gain.embryo.vec[i] <- compute_gain_sim(params, loss.C, embryo.loss.params)$gain # embryo selection   multi.trait.gain.mean
#    gain.mat[i] <- L$gain.mat
  }
    
# Plot: 

jpeg(paste0(figs_dir, 'diseaes_gain_chrom.jpg'))
plot(params$c.vec, gain.mat[,1], xlab="C", ylab="Gain", type="b", ylim = c(1.5*min(gain.mat), max(0, max(gain.mat))), main=paste0("Gain for ", loss.C, " loss"))
lines(params$c.vec, gain.mat[,2], type="b", col="red") # compare to gain just form embryo selection 
legend(0.8 * max(params$c.vec), 0,   lwd=c(2,2), 
       c( "embryo", "chrom"), col=c("black", "red"), cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
dev.off()
#points(params$c.vec, gain.mat[,3], col="blue", pch=3) # compare to gain just form embryo selection 
#legend(0.7 * max(params$c.vec), 0,   lwd=c(2,2), 
#       c( "embryo", "chrom", "bb"), col=c("black", "red", "blue"), pch=c(4,1,3), cex=0.75) #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),


# points(params$c.vec, gain.mat, col="red", xlab="C", ylab="Gain Relaxed")



    