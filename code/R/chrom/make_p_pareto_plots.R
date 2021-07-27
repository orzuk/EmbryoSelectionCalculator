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


figs_dir = "C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\Pareto\\"
start.time <- Sys.time()


col.vec <- c("red", "green", "blue", "orange", "purple", "pink", "yellow", "cyan", "black", "brown")

# Plot p_{k,n} using asymptotic 

max.n <- 10000

c.vec <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9)
num.c <- length(c.vec)
max.k <- ceil(log(max.n) * max(c.vec))

p_k_n_mat <- pareto_P_mat(max.n, max.k)
p_k_n_binary_mat <- pareto_P_binary_mat(max.n, max.k)

binary_flag <- 1
if(binary_flag)
  p_k_n_mat <- p_k_n_binary_mat


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




# try again: 
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
plot(log(n_k_n[,1]), log(p_k_n[,1]), col=col.vec[1], type="l", xlim=c(0, 100), ylim=c(-5,0.2))
for(c in c(2:num.c))
{
  lines(log(n_k_n[,c]), log(p_k_n[,c]), col=col.vec[c])
}
  
  




# Compare pareto for continuous and binary: 




plot(log(p_k_n_mat[,3]))
points(log(p_k_n_binary_mat[,3]), col="red")



