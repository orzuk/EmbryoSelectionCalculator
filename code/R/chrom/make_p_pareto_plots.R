library(ggplot2)
library(reshape2)
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
library(latex2exp)


main_dir = 'C:\\Users\\Or Zuk\\Documents\\GitHub\\EmbryoSelectionCalculator\\code\\R\\chrom\\'

setwd(main_dir) 
source("chrom_select_funcs.R")
source("chrom_select_algs.R")


figs_dir = paste0(main_dir, 'Pareto\\') #  C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\Pareto\\"
start.time <- Sys.time()


col.vec <- c("red", "green", "blue", "orange", "purple", "pink", "yellow", "cyan", "black", "brown")

# Plot p_{k,n} using asymptotic 

max.n <- 10000000 # 00

c.vec <- c(0.2, 0.5, 0.8, 1, 1.3)
num.c <- length(c.vec)
max.k <- ceil(log(max.n) * max(c.vec))

p_k_n_mat <- pareto_P_mat(max.n, max.k)
mu_k_n_mat <- p_k_n_mat * c(1:max.n) # from probabilities to expectations of counts 
q_k_n_mat <- pareto_P_binary_mat(max.n, max.k)
r_k_n_mat <- pareto_P_binary_mat(max.n, max.k, TRUE)

#binary_flag <- 1
#if(binary_flag)
#  p_k_n_mat <- q_k_n_mat


# 1. Fixed k plot: 

# pal <- colorRamp(c("black", "lightgray"))
pal <- colorRamp(c("red", "blue"))
num.k <- 10
# cont.col.vec <- c('black', 'brown4',  'red',  'orange', 'yellow3')

cont.col.vec <- matrix(0, num.k, 3) # rep('', num.k)
for(k in c(1:num.k))
{
  cont.col.vec[k,] <- pal((k-1) / (num.k-1))
}
cont.col.vec <- cont.col.vec/ 255

n.vec <- round(10^seq(0, 5, 0.1))
# jpeg(paste0(figs_dir, 'p_k_n_fixed_k.jpg'))
jpeg('p_k_n_fixed_k.jpg')
#png('p_k_n_fixed_k.png')
plot(log(n.vec), log(mu_k_n_mat[n.vec,k]), col="white", ylim=c(-0.0, 1.01*max(log(mu_k_n_mat[,1:num.k]))), xlab=TeX("log(n)"), ylab=TeX("$log(\\mu_{k,n})$"))
grid()
chr.col.vec <- rep("", num.k)
for(k in c(1:num.k))
{
  chr.col.vec[k] <-  rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3])
  points(log(n.vec), log(mu_k_n_mat[n.vec,k]), pch=10, col = rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3]))  # exact values 
  lines(log(n.vec), (k-1) * log(log(n.vec)) - lfactorial(k-1) ,  col =  rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3]))  # exact values 
}
#cbind(rep("k=", 10), as.character(1:10))
legend(0, 0.99*max(log(mu_k_n_mat[,1:num.k])), lwd=rep(1, num.k),  paste0(rep("k=", num.k), as.character(1:num.k)), col=chr.col.vec, 
       cex=1, box.lwd = 0,box.col = "white",bg = "white") 
dev.off()
# Problem: this creates steps due to rounding 
#old_plot = 0
#if(old_plot)
#{
#  p_k_n <- matrix(0, nrow=max.n, ncol=num.c)
#  ctr <- 1
#  for(c in c.vec)
#  {
#    k.vec <- pmax(1, ceil(c * log(c(1:max.n))))
#    
#    for(n in 1:max.n)
#      p_k_n[n, ctr] <- p_k_n_mat[n, k.vec[n]]
#    ctr <- ctr + 1
#  }
#  
#  p_k_n <- as.data.frame(p_k_n)
#  names(p_k_n) <- as.character(c.vec)
#  p_k_n$n <- c(1:max.n)
#  
#  
#  
#  p_k_n_melt <- melt(p_k_n , id.vars = 'n', variable.name = 'series')
#  
#  ggplot(p_k_n_melt, aes(n,value)) + geom_line(aes(colour = series)) + 
#    scale_x_continuous(trans='log') +
#    scale_y_continuous(trans='log')
#}



# 2. Plot for k_n = c * log(n). 
# Try again: 
p_k_n <- matrix(1, nrow=max.k, ncol=num.c) # 1 is unfilled 
n_k_n <- matrix(0, nrow=max.k, ncol=num.c)
for(k in c(1:max.k))
{
  cur.n.vec <- round(  exp(k / c.vec))
  n_k_n[k,] <- cur.n.vec
  
  for(c in c(1:num.c))
    if(cur.n.vec[c] <= max.n)
      p_k_n[k,c] <- p_k_n_mat[cur.n.vec[c], k]
}
p_k_n <- as.data.frame(p_k_n)
names(p_k_n) <- as.character(c.vec)


# Plot for multiple n values 
# jpeg(paste0(figs_dir, 'p_c_log_n_vs_n.jpg'))
mu_k_n <- p_k_n * n_k_n
plot(log(n_k_n[,1]), log(mu_k_n[,1]), col="white", xlim=c(0, log(max.n)), ylim=c(-0.2, log(max.n)), 
     xlab=TeX("$log(n)$"), ylab=TeX("$log(\\mu_{c log n,n})$"), main=TeX("$\\mu_{c log n, n}$ vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
for(c in c(1:num.c))
{
  max.n.ind <- max(which(n_k_n[,c] <= max.n))
  if(max.n.ind>0)
  { # * log(n_k_n[1:max.n.ind,c])
    points(log(n_k_n[1:max.n.ind,c]), log(mu_k_n[1:max.n.ind,c]), col=col.vec[c])
    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
  }
}
grid()
legend(0, 10, lwd=rep(1, num.c),  as.character(c.vec), col=col.vec[1:num.c], 
       cex=1, box.lwd = 0,box.col = "white",bg = "white") 


# Try ratio
# jpeg(paste0(figs_dir, 'p_c_log_n_vs_n.jpg'))
plot(log(n_k_n[,1]), rep(1, length(n_k_n[,1])), col="white", xlim=c(0, log(max.n)), ylim=c(-0.2, 5), 
     xlab=TeX("$log(n)$"), ylab=TeX("$ratio$"), main=TeX("$\\frac{\\mu_{c log n, n}}{n f_{c log n}(n)}$  vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
for(c in c(1:num.c))
{
  max.n.ind <- max(which(n_k_n[,c] <= max.n))
  if(max.n.ind>0)
  { # * log(n_k_n[1:max.n.ind,c])
    y_vec <- n_k_n[1:max.n.ind,c] ^ (c.vec[c]*(1-log(c.vec[c]))) * sqrt(c.vec[c]) / sqrt(2*pi * log( n_k_n[1:max.n.ind,c] ))
    points(log(n_k_n[1:(max.n.ind),c]), mu_k_n[1:max.n.ind,c] / y_vec[1:(max.n.ind)], col=col.vec[c])
#    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
  }
}
grid()
legend(0, 5, lwd=rep(1, num.c),  as.character(c.vec), col=col.vec[1:num.c], 
       cex=1, box.lwd = 0,box.col = "white",bg = "white") 
# dev.off()


# dev.off()  
  

# Plot p_k_n
plot(log(n_k_n[,1]), rep(1, length(n_k_n[,1])), col="white", xlim=c(0, log(max.n)), ylim=c(-10, 0.1), 
     xlab=TeX("$log(n)$"), ylab=TeX("$log-prob.$"), main=TeX("$log(p_{c log n, n})$  vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
for(c in c(1:num.c))
{
  max.n.ind <- max(which(n_k_n[,c] <= max.n))
  if(max.n.ind>0)
  { # * log(n_k_n[1:max.n.ind,c])
    points(log(n_k_n[1:(max.n.ind),c]), log(p_k_n[1:max.n.ind,c]) , col=col.vec[c])
    #    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
  }
}
legend(0, -6, lwd=rep(1, num.c),  as.character(c.vec), col=col.vec[1:num.c], 
       cex=1, box.lwd = 0,box.col = "white",bg = "white") 



# Figure 2: Comparison of the continuous and binary probability 
jpeg(paste0(figs_dir, 'p_k_n_comparison.jpg'))
col.vec <- c('black', 'brown4',  'red',  'orange', 'yellow3')
# Compare pareto for continuous and binary: 
k.plt <- 5
plot(log(1:max.n), rep(0, max.n), col="white", ylim=c( min(log(p_k_n_mat)), 0.1), 
     xlab="log(n)", ylab=TeX("$log(P_{Pareto})$"))
for(k in c(1:k.plt))
{
  lines(log(1:max.n), log(p_k_n_mat[,k]), type="l", col=col.vec[k])
  lines(log(1:max.n), log(q_k_n_mat[,k]), col=col.vec[k], lty=5)
  lines(log(1:max.n), log(r_k_n_mat[,k]), col=col.vec[k], lty=3)
}
legend(0, -7, lwd=rep(1, k.plt),  as.character(c(1:k.plt)), col=col.vec[1:k.plt], 
       cex=1, box.lwd = 0,box.col = "white",bg = "white") 
# cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
dev.off()



# Compute upper-bound m_k
library(pracma)
k = 1000
log_m_k = (0:(k-1)) * ( log(0:(k-1)) - 1) - lfactorial(0:(k-1))
log_m_k_tag = log( 0:(k-1)) - psi(0, 1:k)




plot(1:k, log_m_k)

plot(1:k, log_m_k_tag)


max.k <- 10
cc <- pareto_alpha_mat_m_vec(max.k)

plot(log(cc$alpha.mat[1,]))

plot(log(cc$alpha.mat[1,]) / (1:max.k))



# Check maxima of geometric
n <- 10000
alpha <- 0.01
iters <- 100

M <- matrix(0, iters, n)
for(i in 1:iters)
{
  G <- rgeom(n, alpha)
  M[i, ] <- cummax(G)
}
# M <- M / iters

med <- apply(M, 2, median)

plot( -log(1-alpha) * med / log(c(1:n)))










