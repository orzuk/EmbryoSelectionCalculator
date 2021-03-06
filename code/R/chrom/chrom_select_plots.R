# Plot performance
# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") 
source("chrom_select_funcs.R")
source("chrom_select_algs.R")


figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\"
start.time <- Sys.time()


# SEt all parameters
params <- c()
params$M <- 22 # try full chromosomes  
params$c.vec <- 2:10
params$T <- 5 # number of traits 
params$iters <- 10
df <- 5 # For wishart distribution

h.ps <- 0.3  # variane explained by the polygenic score 

sigma.blocks = chr.lengths * h.ps
params$sigma.blocks = sigma.blocks
Sigma <- 0.5*diag(params$T) + matrix(0.5, nrow=params$T, ncol=params$T)   # trait-correlations matrix 
Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 
# Loss parameters 
loss.type <- "disease" # stabilizing" # "disease"
loss.params <- c()
loss.params$K <- c(0.01, 0.05, 0.1, 0.2, 0.3) # prevalence of each disease 
loss.params$h.ps <- rep(h.ps, params$T)
loss.params$theta <- c(1, 1, 1, 1, 1)  # importance of each disease 
loss.params$eta <- 0 # negative L2 regularization 
loss.params$n.blocks <- 2
loss.params$cpp <- TRUE  #  TRUE # run in cpp 

gain.embryo.vec <- bb.gain.vec <- gain.vec <- rep(0, length(params$c.vec))  # the gain when selecting embryos (no chromosomes)
run.plots <- 1
params$alg.str <- c("embryo", "branch_and_bound_lipschitz_middle") # ) "branch_and_bound") # "exact" # "branch_and_bound"
n.methods <- length(params$alg.str)
gain.mat <- matrix(rep(0, length(params$c.vec)*n.methods), ncol = n.methods)
#embryo.loss.params = loss.params
#embryo.loss.params$alg.str = "embryo"
#bb.loss.params = loss.params
#bb.loss.params$alg.str = "branch_and_bound"
loss.params$do.checks = 0
if(run.plots)
  for(i in 1:length(params$c.vec))
  {
    params$C <- params$c.vec[i]
    Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
    is.positive.definite(Sigma.K)
    gain.mat[i,]  <- compute_gain_sim(params, loss.type, loss.params)$gain # chromosomal selection
#    bb.gain.vec[i]  <- compute_gain_sim(params, loss.type, bb.loss.params)$gain # chromosomal selection
#    gain.embryo.vec[i] <- compute_gain_sim(params, loss.type, embryo.loss.params)$gain # embryo selection   multi.trait.gain.mean
#    gain.mat[i] <- L$gain.mat
  }


overall.plot.time <- difftime(Sys.time() , start.time, units="secs")
print(paste0("Overall Running Time for Plots (sec.):", overall.plot.time))

# Save results to file: 
save(params, loss.type, loss.params, gain.mat, overall.plot.time, file="disease_gain_chrom.Rdata")

# Plot: 
jpeg(paste0(figs_dir, 'diseaes_gain_chrom.jpg'))
plot(params$c.vec, gain.mat[,1], xlab="C", ylab="Gain", type="b", ylim = c(1.5*min(gain.mat), max(0, max(gain.mat))), main=paste0("Gain for ", loss.type, " loss"))
lines(params$c.vec, gain.mat[,2], type="b", col="red") # compare to gain just form embryo selection 
legend(0.8 * max(params$c.vec), 0,   lwd=c(2,2), 
       c( "embryo", "chrom"), col=c("black", "red"), cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
grid(NULL, NULL, lwd = 2)
dev.off()
#points(params$c.vec, gain.mat[,3], col="blue", pch=3) # compare to gain just form embryo selection 
#legend(0.7 * max(params$c.vec), 0,   lwd=c(2,2), 
#       c( "embryo", "chrom", "bb"), col=c("black", "red", "blue"), pch=c(4,1,3), cex=0.75) #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),


# points(params$c.vec, gain.mat, col="red", xlab="C", ylab="Gain Relaxed")

# Temp debug one simulation: 
##X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
# PRoblem: not solved optimally 
##loss.params$cpp = FALSE
##sol.bb.R <- optimize_C(X, loss.type, loss.params, "branch_and_bound")
##loss.params$cpp = TRUE
##sol.bb.cpp <- optimize_C(X, loss.type, loss.params, "branch_and_bound")

##loss_PS_mat(X[1,,], loss.type, loss.params)
##loss_PS_mat_rcpp(X[1,,], loss.type, loss.params)

#M = matrix(rnorm(100000), ncol=5)
#t = Sys.time(); P = get_pareto_optimal_vecs_rcpp(M); 
#get.pareto.time <- difftime(Sys.time() , t, units="secs") 
#print(paste0("Num pareto: ", length(P$pareto.inds), ", get pareto time: (sec.): ", get.pareto.time))
#t = Sys.time(); P.ecr = get_pareto_optimal_vecs(M)
#get.pareto.ecr.time <- difftime(Sys.time() , t, units="secs") 
#print(paste0("Num pareto ecr: ", length(P.ecr$pareto.inds), ", get pareto ecr time: (sec.): ", get.pareto.ecr.time))



