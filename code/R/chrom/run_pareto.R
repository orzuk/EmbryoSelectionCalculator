# Temp script for computing asymptotic behaviour of pareto function s
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
library(latex2exp)
library(pracma)

Sys.setenv(RETICULATE_PYTHON = "C:\\ProgramData\\Anaconda3")  # SEt python path 
library(reticulate) # calling python function
source_python("pareto.py")



main_dir = 'C:\\Users\\Or Zuk\\Documents\\GitHub\\EmbryoSelectionCalculator\\code\\R\\chrom\\'
setwd(main_dir) 
source("pareto_funcs.R")
source("chrom_select_algs.R")

figs_dir = paste0(main_dir, 'Pareto\\') 

max.n <- 10^4
max.k <- 100

q.mat <- pareto_P_binary_mat(max.n, max.k, FALSE)
r.mat <- pareto_P_binary_mat(max.n, max.k, TRUE)


n.vec <- round(logspace(1, 10, 1000))
# Try different values for k: 
k.vec <- ceil(2.5*log(n.vec))
# k.vec <- round(sqrt(n.vec))
#k.vec <- round(n.vec^(1/4))
# k.vec <- round(log(n.vec)*log(log(n.vec)))
num.n <- length(n.vec)

q.vec <- r.vec <- rep(0, num.n)
for(i in 1:num.n)
{
  q.vec[i] <- pareto_P_binary_log(n.vec[i], k.vec[i], FALSE)
  r.vec[i] <- pareto_P_binary_log(n.vec[i], k.vec[i], TRUE)
}

y.lim <- range((q.vec), (r.vec))
plot(log(n.vec), (q.vec), ylim = y.lim)
points(log(n.vec), (r.vec), col="red", pch="+")
legend(log(n.vec[1]), y.lim[2], lwd=rep(1, 2),  c("q", "r"), col=c("black", "red"), 
       pch = rep(10, 2), lty=rep(1, 2), cex=1, box.lwd = 0,box.col = "white",bg = "white") 


min(q.vec - r.vec)

  
  
  
  
  