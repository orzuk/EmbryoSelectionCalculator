# Plot performance
# Functions for chromosomal selection 
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
library(stringr)
library(Hmisc)

# root_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom"
# figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\"
#figs_dir = paste0(substr(root_dir, 1, tail(str_locate_all(root_dir, "/")[[1]][,1], 3)[1]), "Figures/chrom/")
root_dir = "C:/Code/GitHub/EmbryoSelectionCalculator/code/R/chrom"
figs_dir = paste0(root_dir, "/Figures/")
sim_res_dir = paste(root_dir, "/sim_res/")


setwd(root_dir)
source("chrom_select_funcs.R")
source("chrom_select_algs.R")
# Read all functions for embryo-selection
setwd('../')
embryo.file.sources = list.files(pattern="*.R$", full.names = TRUE)
sapply(embryo.file.sources,source,.GlobalEnv)
setwd(root_dir)


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

params$c.vec <- 2:5
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
# Figure 2.a.: Size of tree for Branch-and-Bound for chromosomal selection
###############################################################
plot_one_sim_BB_num_vectors <- function(params, save.figs = TRUE)
{
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
  
  bb.runtime <- sol.bb$run.time
  sol.bb <- sol.bb$alg.ouptut
  bb.lip.runtime <- sol.bb.lip$run.time
  sol.bb.lip <- sol.bb.lip$alg.ouptut
  
  
  if(any(sol.bb.lip$opt.c != sol.bb$opt.c)) # check that we got the same solution!
  {
    print("ERROR!!! TWO ALGORITHMS GAVE DIFFERENT RESULTS!!!")
  } else
    print("GOOD!!! TWO ALGORITHMS GAVE THE SAME RESULTS!!!")
  
  if(save.figs)
    jpeg(paste0(figs_dir, 'bb_runtime_chrom.jpg'),  res=300) # increase resolution
  
  plot(1:params$M, sol.bb$L.vec, xlab="Chrom.", ylab="Num. Vectors", type="l", log='y', ylim = c(1, params$C**params$M), 
       col="red", yaxt="n", xaxp = c(1, 23, 11), cex=1.25)  # , main=paste0("Number of vectors considered "))
  #pow.lab <- apply( rbind(rep("10^", 8), as.character(c(0:7))), 2, paste, collapse="")
  pow.lab <- parse(text=paste("10^", abs(seq(0, 7, 1)), sep=""))
  axis(2, 10^seq(0L,7L,1L), cex.axis=1.25, labels=pow.lab)
  
  
  lines(1:params$M, sol.bb.lip$L.vec, type="l", col="blue") # compare to gain just form embryo selection 
  lines(1:params$M, params$C ** c(1:params$M), type="l", col="black") # compare to gain just form embryo selection 
  grid(NULL, NULL, lwd = 2)
  legend(1, 0.8*params$C**params$M,   lwd=c(2,2), 
         c(  "Exp.", "Branch&Bound", "Div&Conq"), col=c("black", "red", "blue"), 
         cex=1.25, box.lwd = 0, box.col = "white", bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  if(save.figs)
    dev.off()
}

###############################################################
# Figure 2.b.: average many runs 
###############################################################
plot_BB_num_vectors_errorbars <- function(params, time.iters = 100, save.figs = TRUE, force.rerun = FALSE)
{
  time.iters <- 5 # 100
  params$C = 2 # max(params$c.vec)
  
  params.C.vec <- c(2,4) # 6
  params$max.M <- 10 # 23
  num.m <- params$max.M - 1
  
  sim.res.file <- ''
  
  
  if(force.rerun | !file.exists(sim.res.file))   # Run analysis (could be heavy)  
  {
    dc.num.vecs <- bb.num.vecs <- bb.time.mat <- dc.time.mat <- vector("list", length(params.C.vec))
    for(k in 1:length(params.C.vec))
      dc.num.vecs[[k]] <- bb.num.vecs[[k]] <- bb.time.mat[[k]] <- dc.time.mat[[k]] <- matrix(0, nrow=time.iters, ncol = params$max.M)
    
    for(i in 1:time.iters)
    {
      for(j in 2:params$max.M)
      {
        params$M <- j
        print(paste0("Run iter=", i, ", ", j, " out of ", time.iters, " , ", params$max.M))
        
        for(k in 1:length(params.C.vec))
        {    
          params$C <- params.C.vec[k]
          Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
          X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
          loss.params$lipschitz <- FALSE
          bb.start.time <- Sys.time()
          sol.bb <- optimize_C(X, loss.type, loss.params, "branch_and_bound")
          print(paste0("Overall B&B Running Time (sec.):", difftime(Sys.time() , bb.start.time, units="secs")))
          loss.params$lipschitz <- TRUE
          loss.params$lipschitz.alpha <- lipschitz_loss_PS(loss.type, loss.params)  
          sol.dc <- optimize_C(X, loss.type, loss.params, "branch_and_bound_divide_and_conquer") #optimize_C_branch_and_bound_divide_and_conquer
          
          bb.num.vecs[[k]][i,j] <- max(sol.bb$alg.ouptut$L.vec)
          dc.num.vecs[[k]][i,j] <- max(sol.dc$alg.ouptut$L.vec)
          bb.time.mat[[k]][i,j] <- sol.bb$run.time
          dc.time.mat[[k]][i,j] <- sol.dc$run.time
        }
      } 
    }  # end loop on iters
  } else # here load results of previous run
    load(sim.res.file)
  print("Finished loops!! ")
  
  if(save.figs)
    jpeg(paste0(figs_dir, 'bb_runtime_average.jpg'))
  errbar(1:params$max.M, log10(pmax(1, colMeans(bb.num.vecs))), 
         log10(pmax(1, colMeans(bb.num.vecs)+sqrt(colVars(bb.num.vecs)))), 
         log10(pmax(1, colMeans(bb.num.vecs)-sqrt(colVars(bb.num.vecs)))),
         type='b', main="Num. Vectors", xlab="M", ylab="Num. vectors", 
         col="red", errbar.col="red", yaxt="n", xaxp = c(1, 23, 11), cex=1.25) # , log="y")
  # errbar(1:params$max.M, colMeans(dc.num.vecs), colMeans(dc.num.vecs)+sqrt(colVars(dc.num.vecs)),colMeans(dc.num.vecs)-sqrt(colVars(dc.num.vecs)),
  errbar(1:params$max.M, log10(pmax(1, colMeans(dc.num.vecs))), 
         log10(pmax(1, colMeans(dc.num.vecs)+sqrt(colVars(dc.num.vecs)))), 
         log10(pmax(1, colMeans(dc.num.vecs)-sqrt(colVars(dc.num.vecs)))),
         type='b', col="blue", errbar.col="blue", add="TRUE")
  lines(1:params$max.M, log10(params$C ** c(1:params$M)), type="l", col="black") # compare to gain just form embryo selection 
  axis(2, seq(0L,4L,1L), cex.axis=1.25, labels=parse(text=paste("10^", abs(seq(0, 4, 1)), sep="")))
  
  print("Finished first fig!! ")
  
  
  #plot(2:params$max.M, bb.num.vecs, xlab="Chrom.", ylab="Num. Vectors", type="l", log='y', ylim = c(1, params$C**params$M), 
  #     col="red", main=paste0("Number of vectors considered "))
  #lines(1:params$M, sol.bb.lip$L.vec, type="l", col="blue") # compare to gain just form embryo selection 
  
  grid(NULL, NULL, lwd = 2)
  #legend(1, 0.8*params$C**params$max.M,   lwd=c(2,2), # legend appears in other plot
  #       c(  "Exp.", "Branch&Bound", "Div&Conq"), col=c("black", "red", "blue"), 
  #       cex=0.75, box.lwd = 0, box.col = "white", bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  if(save.figs)
    dev.off()
  
  
  ###############################################################
  # Figure 2.NEW.: Running time
  ###############################################################
  print("Another plot!! ")
  
  if(plot.runtime)
  {
    max.y <- 0
    for(k in 1:length(params.C.vec))
    {
      max.y <- max( max.y, max(colMeans(bb.time.mat[[k]]) + sqrt(colVars(bb.time.mat[[k]]))), 
                    max(colMeans(dc.time.mat[[k]]) + sqrt(colVars(dc.time.mat[[k]]))) )
    }  
    
    print("Another plot max.y!! ")
    
    
    plot(1:params$max.M, log10(pmax(1, colMeans(bb.time.mat[[1]]))), ylim = c(-1, log10(max.y*1.05)), 
         xlab="M", ylab="Run. time (sec.)", cex=1.25) 
    for(k in 1:length(params.C.vec))
    {
      errbar(1:params$max.M, log10(pmax(1, colMeans(bb.time.mat[[k]]))), 
             log10(pmax(1, colMeans(bb.time.mat[[k]])+sqrt(colVars(bb.time.mat[[k]])))), 
             log10(pmax(1, colMeans(bb.time.mat[[k]])-sqrt(colVars(bb.time.mat[[k]])))),
             type='b', main="Num. Vectors", xlab="M", ylab="Run. time (sec.)", 
             col="red", errbar.col="red",  cex=1.25, add="TRUE") # , log="y") # yaxt="n", xaxp = c(1, 23, 11), 
      errbar(1:params$max.M, log10(pmax(1, colMeans(dc.time.mat[[k]]))), 
             log10(pmax(1, colMeans(dc.time.mat[[k]])+sqrt(colVars(dc.time.mat[[k]])))), 
             log10(pmax(1, colMeans(dc.time.mat[[k]])-sqrt(colVars(dc.time.mat[[k]])))),
             type='b', col="blue", errbar.col="blue", add="TRUE")
    } 
    legend(1, 2*max( colMeans(bb.time.mat)),   lwd=c(2,2), 
           c(  "Branch&Bound", "Div&Conq"), col=c( "red", "blue"), 
           cex=1.25, box.lwd = 0, box.col = "white", bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  } # end if plot run time
}  # end function for plotting 


###############################################################
# Figure 2.c.: gain as function of copies, for embryo and chromosomal selection
###############################################################
plot_BB_accuracy <- function(params, time.iters = 100, save.figs = TRUE, force.rerun = FALSE)
{
  if(!("loss.type" %in% names(params)))
    params$loss.type <- "disease" # default
  if(params$loss.type == "disease")
    params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer", "relax", "naive_block_by_block") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  if(params$loss.type == "stabilizing")
    params$alg.str <- c("embryo", "closed_form", "SDR_closed_form", "naive_block_by_block") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  legend.vec <- params$alg.str
  for(i in 1:length(legend.vec))
    legend.vec[i] <- strsplit(params$alg.str[i], "_")[[1]][1]
  
  #params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  #params$alg.str <- c("embryo", "relax") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  if(run.plots)
    gain.res <- compute_gain_sim(params, params$loss.type, loss.params) # chromosomal selection
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
    save(params, loss.type, loss.params, gain.res, overall.plot.time, file=paste0(params$loss.type, "_gain_chrom.Rdata"))
    jpeg(paste0(figs_dir, params$loss.type, '_gain_chrom.jpg'))
  }
  n.algs <- length(params$alg.str)
  plot(params$c.vec, gain.res$gain.mat[,1], xlab="C", ylab="Gain", type="b", col=col.vec[1],
       xaxt="n", cex=1.25, ylim = c(min(0,min(gain.res$gain.mat)), max(0, max(gain.res$gain.mat)))) # 
  #     main=paste0("Gain for ", loss.type, " loss, M=", params$M, " C=", params$C, " T=", params$T))
  for(j in 2:n.algs)
    lines(params$c.vec, gain.res$gain.mat[,j], type="b", col=col.vec[j]) # compare to gain just form embryo selection 
  grid(NULL, NULL, lwd = 2)
  legend(0.7 * max(params$c.vec), 0.275*max(gain.res$gain.mat),   lwd=c(2,2), 
         legend.vec, col=col.vec[1:n.algs], cex=1.25, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  axis(1, 2:5, cex.axis=1.25, labels=c(2:5))
  
  print(legend.vec)
  if(save.figs)
    dev.off()
  
  # save also running times: 
  
} # end plot accuracy
###############################################################


###############################################################
# Figure 2.d.: gain as function of ???, for stabilizing selection loss
###############################################################
# Similar to Figure 2.c. No need for a new function
run.2c = FALSE
if(run.2c)
{
  params$c.vec <- 2:10
  params$iters <- 50
  loss.params$n.blocks = 23
  loss.params$eta <- 0.0
  
  loss.type <- "stabilizing"
  params$M <- 23  # reduce to run fast !! 
  #params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer", "relax") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  #params$alg.str <- c("embryo", "branch_and_bound_divide_and_conquer") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  params$alg.str <- c("embryo", "closed_form") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  if(run.plots)  
    gain.stab <- compute_gain_sim(params, loss.type, loss.params) # chromosomal selection
  
  overall.plot.time <- difftime(Sys.time() , start.time, units="secs")
  print(paste0("Overall Running Time for Plots (sec.):", overall.plot.time))
  if(save.figs) # Plot: 
  {
    # Save results to file: 
    save(params, loss.type, loss.params, gain.stab, overall.plot.time, file="stab_gain_chrom.Rdata")
    jpeg(paste0(figs_dir, 'stab_gain_chrom.jpg'))
  }
  n.algs <- length(params$alg.str)
  plot(params$c.vec, gain.stab$gain.mat[,1], xlab="C", ylab="Gain", type="b", col=col.vec[1],
       ylim = c(min(0,min(gain.stab$gain.mat)), max(0, max(gain.stab$gain.mat)))) 
  #     main=paste0("Gain for ", loss.type, " loss, M=", params$M, " C=", params$C, " T=", params$T))
  for(j in 2:n.algs)
    lines(params$c.vec, gain.stab$gain.mat[,j], type="b", col=col.vec[j+1]) # compare to gain just form embryo selection 
  grid(NULL, NULL, lwd = 2)
  #legend(0.8 * max(params$c.vec), 0.2*max(gain.res$gain.mat),   lwd=c(2,2),  # use legend of previous figure 
  #       c(params$alg.str[1], "relax"), col=col.vec[c(1, 3:(n.algs+1))], cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  if(save.figs)
    dev.off()
  
  
  # New: another plot? save table with running time 
  save("gain.stab", "params", "loss.params", file="temp_stabilizing_results")
}

#########################################################################
# Debug: exact stabilizing worse than embryo?? (for large C)
#########################################################################
debug.code = FALSE
if(debug.code)
{
  params$C <- 6
  Sigma.K <- 0.5*diag(params$C) + matrix(0.5, nrow=params$C, ncol=params$C)   # kinship-correlations matrix 
  X = simulate_PS_chrom_disease_risk(params$M, params$C, params$T, Sigma.T, Sigma.K, sigma.blocks, rep(0.5, k))
  loss.params$C.init <- NULL
  sol.stab <- optimize_C(X, loss.type, loss.params, "closed_form")
  sol.rel <- optimize_C(X, loss.type, loss.params, "relax")
  sol.emb <- optimize_C(X, loss.type, loss.params, "embryo")
  loss.params$eta <- 0.03
  sol.stab.reg <- optimize_C(X, loss.type, loss.params, "closed_form")
  sol.rel.reg <- optimize_C(X, loss.type, loss.params, "relax")
  
  loss_PS(compute_X_C_mat(X, C.mat), loss.type, loss.params) - loss.params$eta * sum(C.mat**2) # TEMP
  loss_PS(compute_X_C_mat(X, sol.stab.reg$C.mat), loss.type, loss.params) - loss.params$eta * sum(sol.stab.reg$C.mat**2)
  loss_PS(compute_X_C_mat(X, sol.rel.reg$C.mat), loss.type, loss.params) - loss.params$eta * sum(sol.rel.reg$C.mat**2)
  loss_PS(compute_X_C_mat(X, sol.stab.reg$C.mat), loss.type, loss.params) - 2*loss.params$eta * sum(sol.stab.reg$C.mat**2)
  loss_PS(compute_X_C_mat(X, sol.rel.reg$C.mat), loss.type, loss.params) - 2*loss.params$eta * sum(sol.rel.reg$C.mat**2)
  
  loss.params$C.init <- sol.stab.reg$C.mat
  sol.rel.reg.mixed <- optimize_C(X, loss.type, loss.params, "relax")
  loss_PS(compute_X_C_mat(X, sol.rel.reg.mixed$C.mat), loss.type, loss.params) - loss.params$eta * sum(sol.rel.reg.mixed$C.mat**2)
  
  
  C.rand <- matrix(0, nrow = params$M, ncol = params$C)
  for(i in 1:params$M)
  {
    C.rand[i, sample(1:params$C, 1)] <- 1
  }
  loss_PS(compute_X_C_mat(X, C.rand), loss.type, loss.params) - loss.params$eta * sum(C.rand**2)
  loss.str <- paste0("L embryo: ", round(sol.emb$opt.loss, 3), 
                     " exact: ", round(sol.stab$opt.loss, 3), 
                     " grad: ", round(sol.rel$opt.loss, 3))
  print(loss.str)
  
  n.rand <- 500
  c.rand <- matrix(sample(1:params$C, n.rand*params$M, replace=TRUE), nrow=n.rand)
  
  X.c.mat <- compute_X_c_vecs(X, c.rand)
  i = 19
  compute_X_c_vec(X, c.rand[i,]) - X.c[i,]
  loss.dist <- loss_PS_mat(X.c.mat, loss.type, loss.params)
  hist(loss.dist, breaks=50, main=loss.str)
  
  
  # Debug SDP:
  # Find brute-force the optimal C
  loss.vec <- rep(0, C**M)
  ctr = 1
  for(a in 1:C)
    for(b in 1:C)
    {
      c.vec <- c(a,b)
      loss.vec[ctr] <- loss_PS(compute_X_c_vec(X, c.vec), loss.type, loss.params)
      ctr <- ctr + 1
    }
  opt.c.vec <- c(3,2)  
  loss_PS(compute_X_c_vec(X, c(3,2)), loss.type, loss.params)
  ret <- SDP_to_integer_solution(X, SDR.ret$X[[1]], loss.type, loss.params, method = "svd")  # randomization")  # "svd"
#  > ret
# (1,1)  
  

}
#########################################################################
# End debug 
#########################################################################




###############################################################
# Figure 3.: linear gain asymptotic vs. simulation 
###############################################################
plot_linear_asymptotic_vs_sim <- function(params, save.figs = TRUE)
{
  
  params$c.vec <- 2:10
  params$iters <- 50
  loss.params$n.blocks = 45
  loss.params$eta <- 0.0
  Sigma.T <- rWishart(1, df, Sigma)[,,1]  # traits correlation matrix 
  
  loss.type <- "quant"  # linear quantitative 
  params$M <- 23  # reduce to run fast !! 
  params$alg.str <- c("embryo", "closed_form") # ) "branch_and_bound") # "exact" # "branch_and_bound"
  
  if(run.plots)  
    gain.linear <- compute_gain_sim(params, loss.type, loss.params) # chromosomal selection
  
  if(save.figs) # Plot: 
  {
    jpeg(paste0(figs_dir, 'gain_linear_quant.jpg'))
  }
  n.algs <- length(params$alg.str)
  plot(params$c.vec, gain.linear$gain.mat[,1], xlab="C", ylab="Gain", type="p", col=col.vec[1],
       ylim = c(min(0,min(gain.linear$gain.mat)), max(0, max(gain.linear$gain.mat)))) 
  for(j in 2:n.algs)
    points(params$c.vec, gain.linear$gain.mat[,j], type="p", col=col.vec[j+1]) # compare to gain just form embryo selection 
  grid(NULL, NULL, lwd = 2)
  legend(0.8 * max(params$c.vec), 0.2*max(gain.linear$gain.mat),   lwd=c(2,2),  # use legend of previous figure 
         c(params$alg.str[1], "chrom"), col=col.vec[c(1, 3:(n.algs+1))], cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
  
  # Add asymptotic formula (depends heavily on Sigma.t)
  gain.linear.approx <- multi_trait_gain_mean(params$c.vec, Sigma.T, loss.params$theta, 'approx')
  lines(params$c.vec, gain.linear.approx, type="l", col=col.vec[j+1]) # compare to gain just from embryo selection 
  
  gain.gamete.linear.approx <- multi_trait_gain_gamete_mean(params$c.vec, Sigma.T, loss.params$theta, 'approx')
  lines(params$c.vec, gain.gamete.linear.approx, type="l", col=col.vec[j+2]) # compare to gain just from gamete selection 
  
  if(save.figs)
    dev.off()
}



