# Plot performance
# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
library(stringr)

# root_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom"
# figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\"
#figs_dir = paste0(substr(root_dir, 1, tail(str_locate_all(root_dir, "/")[[1]][,1], 3)[1]), "Figures/chrom/")
root_dir = "C:/Code/GitHub/EmbryoSelectionCalculator/code/R/chrom"
figs_dir = paste0(root_dir, "/Figures/")

setwd(root_dir)
source("chrom_select_funcs.R")
source("chrom_select_algs.R")

start.time <- Sys.time()

# For plotting
col.vec <- c("red", "green", "blue", "orange", "purple", "pink")
pal <- colorRamp(c("red", "blue"))
num.k <- 5 # length(k.vec)
num.c <- 5 # length(c.vec)
c.col.vec <- matrix(0, num.c, 3) # rep('', num.k)
for(k in c(1:num.c))
{
  c.col.vec[k,] <- pal((k-1) / (num.c-1))
}
c.col.vec <- c.col.vec/ 255
chr.c.col.vec <- rep("", num.c)
for(k in c(1:num.c))
{
  chr.c.col.vec[k] <- rgb(c.col.vec[k,1], c.col.vec[k,2], c.col.vec[k,3])
}


# Set all parameters
params <- c()
params$M <-23 # 10 for fast running 22 # try full chromosomes  
params$T <- 5 # number of traits 
df <- 5 # For Wishart distribution

h.ps <- 0.3  # variance explained by the polygenic score 

params$max.C = max(params$c.vec)
sigma.blocks = c(chr.lengths, chr.lengths) * h.ps  # allow up to 46 
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
loss.params$n.blocks <- 8
loss.params$cpp <- TRUE  # run in cpp (faster)
loss.params$max.L <- 10**6 # maximal number to take in B&B algorithm


gain.embryo.vec <- bb.gain.vec <- gain.vec <- rep(0, length(params$c.vec))  # the gain when selecting embryos (no chromosomes)
run.plots <- 1
n.methods <- length(params$alg.str)
gain.mat <- matrix(rep(0, length(params$c.vec)*n.methods), ncol = n.methods)
#embryo.loss.params = loss.params
#embryo.loss.params$alg.str = "embryo"
#bb.loss.params = loss.params
#bb.loss.params$alg.str = "branch_and_bound"
loss.params$do.checks = 0
save.figs <- FALSE


###############################################################
# Figure 1: Size of tree for Branch-and-Bound for chromosomal selection
###############################################################
params$C = 2 # max(params$c.vec)
Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
loss.params$lipschitz <- FALSE
bb.start.time <- Sys.time()
sol.bb <- optimize_C(X, loss.type, loss.params, "branch_and_bound")
print(paste0("Overall B&B Running Time (sec.):", difftime(Sys.time() , bb.start.time, units="secs")))
loss.params$lipschitz <- TRUE
loss.params$lipschitz.alpha <- lipschitz_loss_PS(loss.type, loss.params)  
sol.bb.lip <- optimize_C(X, loss.type, loss.params, "branch_and_bound_divide_and_conquer") #optimize_C_branch_and_bound_divide_and_conquer

if(any(sol.bb.lip$opt.c != sol.bb$opt.c)) # check that we got the same solution!
{
  print("ERROR!!! TWO ALGORITHMS GAVE DIFFERENT RESULTS!!!")
} else
  print("GOOD!!! TWO ALGORITHMS GAVE THE SAME RESULTS!!!")

if(save.figs)
  jpeg(paste0(figs_dir, 'bb_runtime_chrom.jpg'))
plot(1:params$M, sol.bb$L.vec, xlab="Chrom.", ylab="Num. Vectors", type="l", log='y', ylim = c(1, params$C**params$M), 
      col="red", main=paste0("Number of vectors considered "))
lines(1:params$M, sol.bb.lip$L.vec, type="l", col="blue") # compare to gain just form embryo selection 
lines(1:params$M, params$C ** c(1:params$M), type="l", col="black") # compare to gain just form embryo selection 
grid(NULL, NULL, lwd = 2)
legend(1, 0.8*params$C**params$M,   lwd=c(2,2), 
       c(  "Exp.", "Branch&Bound", "Div&Conq"), col=c("black", "red", "blue"), 
       cex=0.75, box.lwd = 0, box.col = "white", bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
if(save.figs)
  dev.off()


# Figure 1.b.: average many runs 
time.iters <- 3
params$C = 2 # max(params$c.vec)

params$max.M <- 10
num.m <- params$max.M - 1
dc.num.vecs <- bb.num.vecs <- matrix(0, nrow=time.iters, ncol = params$max.M)
for(i in 1:time.iters)
{
  for(j in 2:params$max.M)
  {
    params$M <- j
    print(paste0("Run iter=", i, ", ", j, " out of ", time.iters, " , ", params$max.M))
    Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
    X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
    loss.params$lipschitz <- FALSE
    bb.start.time <- Sys.time()
    sol.bb <- optimize_C(X, loss.type, loss.params, "branch_and_bound")
    print(paste0("Overall B&B Running Time (sec.):", difftime(Sys.time() , bb.start.time, units="secs")))
    loss.params$lipschitz <- TRUE
    loss.params$lipschitz.alpha <- lipschitz_loss_PS(loss.type, loss.params)  
    sol.dc <- optimize_C(X, loss.type, loss.params, "branch_and_bound_divide_and_conquer") #optimize_C_branch_and_bound_divide_and_conquer
    
    bb.num.vecs[i,j] <- max(sol.bb$L.vec)
    dc.num.vecs[i,j] <- max(sol.dc$L.vec)
  } 
}

###############################################################
# Figure 2: gain as function of copies, for embryo and chromosomal selection
###############################################################
params$c.vec <- 2:10
params$iters <- 10
loss.params$n.blocks = 4 
params$M <- 23  # reduce to run fast !! 
params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer", "relax") # ) "branch_and_bound") # "exact" # "branch_and_bound"
params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer") # ) "branch_and_bound") # "exact" # "branch_and_bound"
params$alg.str <- c("embryo", "relax") # ) "branch_and_bound") # "exact" # "branch_and_bound"
if(run.plots)
  gain.res <- compute_gain_sim(params, loss.type, loss.params) # chromosomal selection
#  for(i in 1:length(params$c.vec))
#  {
#    params$C <- params$c.vec[i]
#    bb.gain.vec[i]  <- compute_gain_sim(params, loss.type, bb.loss.params)$gain # chromosomal selection
#    gain.embryo.vec[i] <- compute_gain_sim(params, loss.type, embryo.loss.params)$gain # embryo selection  multi.trait.gain.mean
#    gain.mat[i] <- L$gain.mat
#  }

overall.plot.time <- difftime(Sys.time() , start.time, units="secs")
print(paste0("Overall Running Time for Plots (sec.):", overall.plot.time))


# Plot: 
if(save.figs)
{
  # Save results to file: 
  save(params, loss.type, loss.params, gain.res, overall.plot.time, file="disease_gain_chrom.Rdata")
  jpeg(paste0(figs_dir, 'diseaes_gain_chrom.jpg'))
}
n.algs <- length(params$alg.str)
plot(params$c.vec, gain.res$gain.mat[,1], xlab="C", ylab="Gain", type="b", col=col.vec[1],
     ylim = c(min(0,min(gain.res$gain.mat)), max(0, max(gain.res$gain.mat))), 
     main=paste0("Gain for ", loss.type, " loss, M=", params$M, " C=", params$C, " T=", params$T))
for(j in 2:n.algs)
  lines(params$c.vec, gain.res$gain.mat[,j], type="b", col=col.vec[j]) # compare to gain just form embryo selection 
grid(NULL, NULL, lwd = 2)
legend(0.8 * max(params$c.vec), 0.2*max(gain.res$gain.mat),   lwd=c(2,2), 
       c( "embryo", "chrom", "relax"), col=col.vec[1:n.algs], cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
if(save.figs)
  dev.off()
###############################################################


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



