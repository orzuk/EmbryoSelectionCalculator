library(ggplot2)
library(reshape2)
library(mvtnorm)
library(rWishart)
library(matrixcalc)
library(matrixNormal)
library(Rfast)
library(ecr)
library(latex2exp)

Sys.setenv(RETICULATE_PYTHON = "C:\\ProgramData\\Anaconda3")  # SEt python path 
library(reticulate) # calling python function
source_python("pareto.py")



main_dir = 'C:\\Users\\Or Zuk\\Documents\\GitHub\\EmbryoSelectionCalculator\\code\\R\\chrom\\'
setwd(main_dir) 
source("pareto_funcs.R")
source("chrom_select_algs.R")

figs_dir = paste0(main_dir, 'Pareto\\') #  C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Figures\\chrom\\Pareto\\"
start.time <- Sys.time()

col.vec <- c("red", "green", "blue", "orange", "purple", "pink", "yellow", "cyan", "black", "brown")

max.n <- 10000000 # 00

# c.vec <- c(0.2, 0.5, 0.8, 1, 1.3)
c.vec <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
num.c <- length(c.vec)
max.k <- ceil(log(max.n) * max(c.vec))

p_k_n_mat <- pareto_P_mat(max.n, max.k)
mu_k_n_mat <- p_k_n_mat * c(1:max.n) # from probabilities to expectations of counts 
q_k_n_mat <- pareto_P_binary_mat(max.n, max.k)
r_k_n_mat <- pareto_P_binary_mat(max.n, max.k, TRUE)

#binary_flag <- 1
#if(binary_flag)
#  p_k_n_mat <- q_k_n_mat


# Figure 1: Fixed k plot: 
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

n.vec <- round(10^seq(0, 7, 0.1))
max.y <- max(log(mu_k_n_mat[,1:num.k]))
png(paste0(figs_dir, 'p_k_n_fixed_k.png')  , res=300,  width = 300, height=300, units="mm")
par(mar=c(5,6,4,1)+.1) # increase left margin
plot(log(n.vec), log(mu_k_n_mat[n.vec,k]), col="white", ylim=c(-0.0, 1.01*max.y), 
     xlab=TeX("log(n)"), ylab=TeX("$log(\\mu_{k,n})$"), cex.lab=2, cex.axis=2)
grid()
chr.col.vec <- rep("", num.k)
for(k in c(1:num.k))
{
  chr.col.vec[k] <-  rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3])
  points(log(n.vec), log(mu_k_n_mat[n.vec,k]), pch=10, col = rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3]))  # exact values 
  lines(log(n.vec), (k-1) * log(log(n.vec)) - lfactorial(k-1) ,  col =  rgb(cont.col.vec[k,1], cont.col.vec[k,2], cont.col.vec[k,3]))  # exact values 
}
legend(0, 0.99*max.y, lwd=rep(1, num.k),  paste0(rep("k=", num.k), as.character(1:num.k)), col=chr.col.vec, 
       cex=2, box.lwd = 0,box.col = "white",bg = "white") 
dev.off()


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


# # Plot for multiple n values 
# # jpeg(paste0(figs_dir, 'p_c_log_n_vs_n.jpg'))
# mu_k_n <- p_k_n * n_k_n
# plot(log(n_k_n[,1]), log(mu_k_n[,1]), col="white", xlim=c(0, log(max.n)), ylim=c(-0.2, log(max.n)), 
#      xlab=TeX("$log(n)$"), ylab=TeX("$log(\\mu_{c log n,n})$"), main=TeX("$\\mu_{c log n, n}$ vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
# for(c in c(1:num.c))
# {
#   max.n.ind <- max(which(n_k_n[,c] <= max.n))
#  if(max.n.ind>0)
#  { # * log(n_k_n[1:max.n.ind,c])
#    points(log(n_k_n[1:max.n.ind,c]), log(mu_k_n[1:max.n.ind,c]), col=col.vec[c])
#    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
#  }
# }
# grid()
# legend(0, 10, lwd=rep(1, num.c),  as.character(c.vec), col=col.vec[1:num.c], 
#       cex=1, box.lwd = 0,box.col = "white",bg = "white") 
#
#
# # Try ratio
# # jpeg(paste0(figs_dir, 'p_c_log_n_vs_n.jpg'))
# plot(log(n_k_n[,1]), rep(1, length(n_k_n[,1])), col="white", xlim=c(0, log(max.n)), ylim=c(-0.2, 5), 
#     xlab=TeX("$log(n)$"), ylab=TeX("$ratio$"), main=TeX("$\\frac{\\mu_{c log n, n}}{n f_{c log n}(n)}$  vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
# for(c in c(1:num.c))
# {
#  max.n.ind <- max(which(n_k_n[,c] <= max.n))
#  if(max.n.ind>0)
#  { # * log(n_k_n[1:max.n.ind,c])
#    y_vec <- n_k_n[1:max.n.ind,c] ^ (c.vec[c]*(1-log(c.vec[c]))) * sqrt(c.vec[c]) / sqrt(2*pi * log( n_k_n[1:max.n.ind,c] ))
#    points(log(n_k_n[1:(max.n.ind),c]), mu_k_n[1:max.n.ind,c] / y_vec[1:(max.n.ind)], col=col.vec[c])
##    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
#  }
# }
# grid()
# legend(0, 5, lwd=rep(1, num.c),  as.character(c.vec), col=col.vec[1:num.c], 
#       cex=1, box.lwd = 0,box.col = "white",bg = "white") 
# # dev.off()



# Figure 2: Plot p_k_n for k_n ~ c * log(n)
c.col.vec <- matrix(0, num.c, 3) # rep('', num.k)
chr.c.col.vec <- rep("", num.c)
for(c in c(1:num.c))
{
  c.col.vec[c,] <- (pal((c-1) / (num.c-1))) / 255
  chr.c.col.vec[c] <-  rgb(c.col.vec[c,1], c.col.vec[c,2], c.col.vec[c,3])
}

# jpeg(paste0(figs_dir, 'p_clogn_n.jpg'))
png(paste0(figs_dir, 'p_clogn_n.png')  , res=300,  width = 300, height=300, units="mm")
par(mar=c(5,6,4,1)+.1) # increase left margin
plot(log(n_k_n[,1]), rep(1, length(n_k_n[,1])), col="white", xlim=c(0, log(max.n)), ylim=c(-6, 0.01), 
     xlab=TeX("$log(n)$"), ylab=TeX("$log(p_{c log n, n})$"), cex.lab=2, cex.axis=2) # , main=TeX("$log(p_{c log n, n})$  vs. $n$")) # type="l", # log(max(mu_k_n[!is.infinite(mu_k_n)]))
grid()
for(c in c(1:num.c))
{
  max.n.ind <- max(which(n_k_n[,c] <= max.n))
  if(max.n.ind>0)
  { # * log(n_k_n[1:max.n.ind,c])
    points(log(n_k_n[1:(max.n.ind),c]), log(p_k_n[1:max.n.ind,c]) , col=chr.c.col.vec[c], pch=10)
    lines(log(n_k_n[1:(max.n.ind),c]), log(p_k_n[1:max.n.ind,c]) , col=chr.c.col.vec[c])
    #    lines(log(n_k_n[1:max.n.ind,c]), c.vec[c]*(1-log(c.vec[c])) * log(n_k_n[1:max.n.ind,c]) - 0.5 * log (2*pi * log(n_k_n[1:max.n.ind,c]) ), col=col.vec[c]) # Add approximation line 
  }
}
legend(0, -4, lwd=rep(1, num.c),  paste0(rep("c=", num.c), as.character(c.vec)), col=chr.c.col.vec[1:num.c], 
       pch = rep(10, num.c), lty=rep(1, num.c), cex=2, box.lwd = 0,box.col = "white",bg = "white") 
dev.off()


# Figure 3: Comparison of the continuous and binary probability 
#jpeg(paste0(figs_dir, 'p_k_n_comparison.jpg'))
min.y <- min(log(p_k_n_mat))
png(paste0(figs_dir, 'p_k_n_comparison.png')  , res=300,  width = 300, height=300, units="mm")
par(mar=c(5,6,4,1)+.1) # increase left margin
col.vec <- c('black', 'brown4',  'red',  'orange', 'yellow3')
# Compare pareto for continuous and binary: 
k.plt <- 5
k.col.vec <- matrix(0, k.plt, 3) # rep('', num.k)
chr.k.col.vec <- rep("", k.plt)
for(k in c(1:k.plt))
{
  k.col.vec[k,] <- (pal((k-1) / (k.plt-1))) / 255
  chr.k.col.vec[k] <-  rgb(k.col.vec[k,1], k.col.vec[k,2], k.col.vec[k,3])
}
plot(log(n.vec), rep(0, length(n.vec)), col="white", ylim=c(min.y, 0.1), 
     xlab="log(n)", ylab=TeX("$log(P_{Pareto})$"), cex.lab=2, cex.axis=2)
grid()
for(k in c(1:k.plt))
{
  lines(log(n.vec), log(p_k_n_mat[n.vec,k]), type="l", col=chr.k.col.vec[k])
  lines(log(n.vec), log(q_k_n_mat[n.vec,k]), col=chr.k.col.vec[k], lty=5)
  lines(log(n.vec), log(r_k_n_mat[n.vec,k]), col=chr.k.col.vec[k], lty=3)
}
legend(-0.5, -12, lwd=rep(1, k.plt),  paste0(rep("k=", 5), as.character(c(1:k.plt))), 
       col=chr.k.col.vec[1:k.plt], cex=2, box.lwd = 0,box.col = "white",bg = "white") 
# cex=0.75, box.lwd = 0,box.col = "white",bg = "white") #  y.intersp=0.8, cex=0.6) #  lwd=c(2,2),
dev.off()



# 4.Figure 4: Compare distribution to Poisson
n.vec <- c(100, 200, 400) # Computing var takes a really long time 
# n.vec <- c(20, 40, 60) # A small example 
k.vec <- c(3, 5, 7)
iters <- 100000
#n.vec <- c(100)
#k.vec <- c(3)



# e_k_n_big_tab = matrix(-1, 20, 1000)  # Store all pre-computed values 

# p_k_n_big_tab = matrix(-1, 20, 1000)  # Store all pre-computed values 
# p_k_n_big_tab = pareto_P_mat(10000, 200)
# save(p_k_n_big_tab, file = "p_k_n_big_tab.Rdata")
# save(e_k_n_big_tab, file = "e_k_n_big_tab.Rdata")
# load("e_k_n_big_tab.Rdata")

png(paste0(figs_dir, 'p_k_n_dist.png'), 
            res=300,  width = 300, height=300, units="mm")
par(mfrow=c(3,3), mai = c(0.0, 0.0, 0.0, 0.0), mar=c(4,6,2,0)+.1)  # mar=c(5,6,4,1)+.1)  # For subplots  # mai = c(1, 0.1, 0.1, 0.1)
# par(mfrow=c(3,3))  # For subplots  # mai = c(1, 0.1, 0.1, 0.1)
for(n in n.vec)
  for(k in k.vec)
  {
    print(paste0("Run n=", n, " k=", k))
    sim.file.name <- paste0("sim_n_", n, "_k_", k, "_iters_", as.integer(iters), ".Rdata")
    
    if(file.exists(sim.file.name))
      load(sim.file.name)
    else
    {
      simpar <- pareto_P_sim(n, k, iters)
      simrange <- min(simpar$n.pareto):max(simpar$n.pareto)
      save(simpar, simrange, n, k, iters, 
           file = paste0("sim_n_", n, "_k_", k, "_iters_", as.integer(iters), ".Rdata"))
    }
    
    
    print("Compute Var combinatorial (int):")
    my.n <- n # must update outside function!!
    V <- pareto_P_var(k, n)
    print("Computed Var!!!")

    dist.pareto <- table(simpar$n.pareto)/iters
    dist.pois <- dpois(simrange, simpar$p.pareto*n)
    dist.norm <- dnorm(simrange, mean = simpar$p.pareto*n, sd = sqrt(V$V))
    y.max <- max(max(dist.pareto), max(dist.pois), max(dist.norm))*1.01
    
        
#    png(paste0(figs_dir, 'p_', as.character(k), '_', as.character(n), '_dist.png'), 
#        res=300,  width = 300, height=300, units="mm")
#    par(mar=c(5,6,4,1)+.1) # increase left margin
    plot(as.numeric(names(dist.pareto)), as.numeric(dist.pareto), ylim=c(0, y.max), 
         pch=19, cex.lab=2, cex.axis=2, xlab="", ylab="")
    if(n == n.vec[3])
      title(xlab=TeX(paste0("$k=", as.character(k), "$")), cex.lab=2) 
    if(k == k.vec[1])
      title(ylab=TeX(paste0("$n=", as.character(n), "$")), cex.lab=2)
    grid()
    lines(simrange, dist.pois, col="red", lwd=1.5)
    lines(simrange, dist.norm, col="blue", lwd=1.5)# new: normal approximation    
    
    if((n == n.vec[3]) & (k == k.vec[3])) # add legend
      legend(min(simrange), 0.99*y.max, lwd=rep(1, 3), c("Sim.", "Poisson", "Gaussian"), cex=2, 
             col=c("black", "red", "blue"), lty=c(NA,1,1), pch=c(19, NA, NA), box.lwd = 0,box.col = "white",bg = "white")
  
#    legend(min(simrange), 0.99*max.y, lwd=rep(1, num.k),  paste0(rep("k=", num.k), as.character(1:num.k)), col=chr.col.vec, 
#           cex=2, box.lwd = 0,box.col = "white",bg = "white") 
    
      
#    dev.off()
  }
dev.off()

  

#####################################################
### Ended figures ### 
#####################################################

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







# Test function: 
k=3; n=3;  pareto_E_Z1Z2(k,n, FALSE);  pareto_E_Z1Z2(k,n, TRUE)$e_k_n; pareto_E_Z1Z2(k,n, TRUE, TRUE)$e_k_n;

k=2; n=10;  A <-  pareto_E_Z1Z2(k,n, TRUE); B <-  pareto_E_Z1Z2(k,n, FALSE); # C <- pareto_E_Z1Z2(k,n, TRUE, FALSE, TRUE); #  B <- pareto_E_Z1Z2(k,n, TRUE, TRUE)$e_k_n;
print(A$e_k_n)
print(B$e_k_n)
print(C$e_k_n)
print(C$run.time)

pareto_P2(n, k)
pareto_P2(n, k)^2


plot(sort(abs(A$e_k_n_vec)))

I <- order(-abs(A$e_k_n_vec))
ev <- A$e_k_n_vec[I]
NEG <- which(ev < 0)

plot(abs(ev))
points(NEG, abs(ev[NEG]), col="red")

cv <- cumsum(ev)
plot(cumsum(ev))

k=2; n=10;  A <- pareto_E_Z1Z2(k,n, TRUE); print(A$e_k_n); B <- pareto_E_Z1Z2_python(as.integer(k), as.integer(n)); print(as.numeric(as.character(B))) 
k=2; n=20;  A <- pareto_E_Z1Z2(k,n, TRUE); print(A$e_k_n)
k=2; n=30;  A <- pareto_E_Z1Z2(k,n, TRUE); print(A$e_k_n)
k=2; n=40;  A <- pareto_E_Z1Z2(k,n, TRUE); print(A$e_k_n)







