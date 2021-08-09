library(ggplot2)
library(reshape2)
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
setwd("C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\chrom") 
source("chrom_select_funcs.R")
source("chrom_select_algs.R")


figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\Pareto\\"
start.time <- Sys.time()


col.vec <- c("red", "green", "blue", "orange", "purple", "pink", "yellow", "cyan", "black", "brown")

# Plot p_{k,n} using asymptotic 

max.n <- 10000000

c.vec <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9)
num.c <- length(c.vec)
max.k <- ceil(log(max.n) * max(c.vec))

p_k_n_mat <- pareto_P_mat(max.n, max.k)
q_k_n_mat <- pareto_P_binary_mat(max.n, max.k)
r_k_n_mat <- pareto_P_binary_mat(max.n, max.k, TRUE)


binary_flag <- 1
if(binary_flag)
  p_k_n_mat <- q_k_n_mat


# Problem: this creates steps due to rounding 
p_k_n <- matrix(0, nrow=max.n, ncol=num.c)
ctr <- 1
for(c in c.vec)
{
  k.vec <- pmax(1, ceil(c * log(c(1:max.n))))
  
  for(n in 1:max.n)
    p_k_n[n, ctr] <- p_k_n_mat[n, k.vec[n]]
  ctr <- ctr + 1
}

p_k_n <- as.data.frame(p_k_n)
names(p_k_n) <- as.character(c.vec)
p_k_n$n <- c(1:max.n)



p_k_n_melt <- melt(p_k_n , id.vars = 'n', variable.name = 'series')

ggplot(p_k_n_melt, aes(n,value)) + geom_line(aes(colour = series)) + 
  scale_x_continuous(trans='log') +
  scale_y_continuous(trans='log')




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
plot(log(n_k_n[,1]), log(p_k_n[,1]), col=col.vec[1], xlim=c(0, log(max.n)), ylim=c(-0.2, log(max.n)), 
     xlab="log(n)", ylab="log(mu_{k,n})", main="Prob vs. asymptote") # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
for(c in c(2:num.c))
{
  max.n.ind <- max(which(n_k_n[,c] <= max.n))
  if(max.n.ind>0)
  {
    points(log(n_k_n[1:max.n.ind,c]), log(mu_k_n[1:max.n.ind,c]), col=col.vec[c])
    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]), col=col.vec[c]) # Add approximation line 
  }
}
#legend(0, 15, lwd=rep(1, k.plt),  as.character(c(1:k.plt)), col=col.vec[1:k.plt], 
#       cex=1, box.lwd = 0,box.col = "white",bg = "white") 

# dev.off()  
  


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



