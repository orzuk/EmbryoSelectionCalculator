# Temp script for computing asymptotic behavior of pareto function s
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



main_dir = 'C:\\Users\\Or Zuk\\Documents\\GitHub\\EmbryoSelectionCalculator\\code\\R\\chrom\\'
setwd(main_dir) 
source_python("pareto.py")
source("pareto_funcs.R")
source("chrom_select_algs.R")

figs_dir = paste0(main_dir, 'Pareto\\') 


col.vec <- c("red", "green", "blue", "orange", "purple", "pink")
pal <- colorRamp(c("red", "blue"))

max.n <- 10^4
max.k <- 100

q.mat <- pareto_P_binary_mat(max.n, max.k, FALSE)
r.mat <- pareto_P_binary_mat(max.n, max.k, TRUE)  # r.mat (strong-front) < q.mat (weak-front)


n.vec <- round(logspace(1, 60, 1000))

k.vec <- 1:100 # fix k 




# Try different values for k: 
#k.vec <- ceil(0.5*log(n.vec))
# k.vec <- round(sqrt(n.vec))
#k.vec <- round(n.vec^(1/4))
# k.vec <- round(log(n.vec)*log(log(n.vec)))
png(paste0(figs_dir, "DiscretePhaseTransition.png"), res=300,  width = 300, height=300, units = "mm")
p.ber <- 0.5
gamma.ber <- -p.ber*log(p.ber)
c.vec <- round(gamma.ber*100)/100 + (-3:2)/100 # c(0.32, 0.33, 0.34, 0.35, 0.36, 0.37)
num.k <- length(k.vec)
num.c <- length(c.vec)
q.vec <- r.vec <- rep(0, num.k)

c.col.vec <- matrix(0, num.c, 3) # rep('', num.k)
for(k in c(1:num.c))
{
  c.col.vec[k,] <- pal((k-1) / (num.c-1))
}
c.col.vec <- c.col.vec/ 255

chr.c.col.vec <- rep("", num.c)
for(k in c(1:num.c))
{
  chr.c.col.vec[k] <-  rgb(c.col.vec[k,1], c.col.vec[k,2], c.col.vec[k,3])
}


q.mat <- r.mat <-  matrix(0, nrow=num.c, ncol=num.k)

for(c in 1:length(c.vec))
{
  n.vec <- ceil(exp(k.vec * c.vec[c]))
  for(i in 1:num.k)
  {
    q.mat[c,i] <- pareto_P_binary_log(n.vec[i], k.vec[i], FALSE, p.ber)$log.p.n.k
    r.mat[c,i] <- pareto_P_binary_log(n.vec[i], k.vec[i], TRUE, p.ber)$log.p.n.k
  }
}
y.lim <- range(r.mat) # y.lim <- range(q.vec)
par(mar=c(5,6,4,1)+.1) # increase left margin
for(c in 1:length(c.vec))
{
  cur.col <- rgb(c.col.vec[c,1], c.col.vec[c,2], c.col.vec[c,3])
  if(c == 1)    
  {
#    y.lim[2] <- max(0, y.lim[2])
    plot(k.vec, r.mat[c,], ylim = y.lim, col=cur.col, pch=1, cex=1, xlab="k", ylab=TeX("$log(p_{k,n})$"), cex.lab=2, cex.axis=2)
  }
  else
    points(k.vec, r.mat[c,], pch=1, cex=1,  col=cur.col)
  # Add q: 
  points(k.vec, q.mat[c,], pch=4, cex=1,  col=cur.col)
}
c.legend <- as.character(c.vec)
for(i in 1:length(c.vec))
  c.legend[i] <- paste0("c=", c.legend[i])

grid(col = "darkgray", lwd=1.5)
legend(0.04*max(k.vec), 0.8*y.lim[1], lwd=rep(1, length(c.vec)),  
       c.legend, col=chr.c.col.vec, 
       pch = rep(10, num.c), lty=rep(1, num.c), cex=2, box.lwd = 0,box.col = "white",bg = "white") 
mtext("c.", side=3, adj=0, line=1.2, cex=4, font=2); 
dev.off()


#y.lim <- range((q.vec), (r.vec)) # y.lim <- range(q.vec)
#plot(log(n.vec), q.vec, ylim = y.lim)
#points(log(n.vec), r.vec, col="red", pch="+")
#legend(log(n.vec[1]), y.lim[2], lwd=rep(1, 2),  c("q", "r"), col=c("black", "red"), 
#       pch = rep(10, 2), lty=rep(1, 2), cex=1, box.lwd = 0,box.col = "white",bg = "white") 


min(q.vec - r.vec)

  
  
  
  
  
#######################################  NEW STUFF ##########################

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

