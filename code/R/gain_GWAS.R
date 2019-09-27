library(MASS)
library(nleqslv)
library(pracma)
source('gain_moments.R')


# Predict the mean gain expected if we perform a GWAS of a certain sample size 
# h2snp - total SNP heritability
# r2ps - variance explained by the polygenic score from GWAS 
# M - effective number of markers (default: -1 is computed inside function)
# N - GWAS sample size 
# N - sample size for which we want to predict
# n - # of embryos
#
mean.gain.from.GWAS <- function(h2snp, r2ps, N_discovery, M=-1,  N, n)
{
  # First compute the expected variance explained in a GWAS of size N: 
  if(M==-1)
    M <- N_discovery * h2snp * (h2snp - r2ps) / r2ps # evaluate effective number of SNPs
  r2 <- h2snp / (1 + M/(N*h2snp)) # variance explained for a given N
  
  # Next, evaluate the gain: 
  return( list(E.G=gain.moments(n, r2, 'exact', 1)$E.G, M=M) )
}

# Usage example (Supp. Fig. 3 from thepaper)
plot_gain <- function(IQ, WGS, label, n)
{
  maxN = 8
  M=-1
  if (IQ)
  {
    h2snp = 0.19
    sigma_z = 15
    known_r2 = 0.052
    known_N = 269867
    ylabel = "Mean gain (points)"
    tstr = sprintf("IQ (array)")
    fname = "IQ"
  } else {
    fname = "Height_chip"
    h2snp = 0.483
    sigma_z = 6
    known_r2 = 0.244
    known_N = 693529
    tstr = sprintf("Height (array)")
    if(WGS) # whole genome sequencing
    {
      M = 10 *  mean.gain.from.GWAS(h2snp, known_r2, known_N, -1, 1, n)$M # assume 10-fold increase in effective # of markers 
      known_N = 21620
      h2snp = 0.79 # heritability including rare variants 
      fname = "Height_wgs"
      tstr = sprintf("Height (WGS)")
    }
    ylabel = "Mean gain (cm)"
  }
  
  Ns <- logspace(4, maxN, 50)
  G <- numeric(length(Ns))
  for(i in seq_along(Ns))
    G[i] <- sigma_z * mean.gain.from.GWAS(h2snp, known_r2, known_N, M, Ns[i], n)$E.G # use function 

  xticks = 10^(4:maxN)
  par(mar=c(5,5.5,2.5,2))
  plot(Ns,G,type="l",col="blue",xlim=c(1e4,max(Ns)),lwd=2,xlab="N (Discovery GWAS)",ylab=ylabel,cex.lab=2,log='x',xaxt="n",ylim=c(0,8),yaxs='i',yaxt="n")
  axis(1,xticks,parse(text=sprintf('10^%d',as.integer(log10(xticks)))),las=2,cex.axis=2,las=0)
  axis(2,0:8,cex.axis=2,las=0)
  title(tstr,font.main=1,cex.main=2)
  ind = which.min(abs(Ns-known_N))
  segments(Ns[ind],0,Ns[ind],G[ind],col="red",lwd=1.5)
  segments(1,G[ind],Ns[ind],G[ind],col="red",lwd=1.5)
}

# Plot the 3 scenarios 
n = 10
par(mfrow=c(1,3))
plot_gain(0,0,'A', n)
plot_gain(1,0,'B', n)
plot_gain(0,1,'C', n)
