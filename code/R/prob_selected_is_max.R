# Probability that selected embryo has maximal phenotype 
# Y ~ N(0,1), epsilon ~ N(0, sqrt(h2ps/(2 - h2z - h2ps))) 
# n - # of embryos 
# h2z - heritability of trait
# h2ps - heritability explained by polygenic score PS
# method.str - how to compute (numeric or simulations)
#
library(MASS)
library(mvtnorm) # for multidimensional Gaussians
library(cubature) # for multidimensional integration 
source('truncated_sum_norm_cdf.R')


# Compute probability that the selected embryo yields the maximal trait value
prob.selected.is.max <- function(n, h2z, h2ps, method.str='numeric')
{	
  if(method.str == 'numeric')
  {
    sigma <- sqrt(h2ps/(2 - h2z - h2ps))
    if(h2ps < 0.0000000001) # avoid numerical problems 
      return(1/n)
    return(hcubature(prob.selected.integrand, lower=c(-20,-20*sigma), upper=c(20,20*sigma), n=n, h2z=h2z, h2ps=h2ps)$integral)
  }
  else # simulations 
  {
    iters <- 500000
    Sigma.PS <- 0.5 * h2ps * (matrix(1, n, n) + diag(n))
    Sigma.eps <- 0.5* ((h2z-h2ps) * matrix(1, n, n) + (2-h2ps-h2z) * diag(n))
    PS <- mvrnorm(n = iters, mu=rep(0,n), Sigma=Sigma.PS)
    eps <- mvrnorm(n = iters, mu=rep(0,n), Sigma=Sigma.eps)
    z <- PS + eps
    return( mean(max.col(PS) == max.col(z)) ) # find indices maximizing and compare
  }
}


# Internal function: Probability of correct selection conditioned on y_max, eps of y_max
prob.selected.integrand <- function(arg.vec, n, h2z, h2ps)
{
  eps <- arg.vec[1]
  y <- arg.vec[2]
  sigma <- sqrt(h2ps/(2 - h2z - h2ps))
  return(  (n * pnorm(y/sigma)^(n-1) * dnorm(y/sigma) * dnorm(eps) / sigma) * 
             truncated.sum.norm.cdf(y+eps, h2z, h2ps, y)^(n-1) )
} 

h2ps <- 0.27
h2z <- 0.8
n <- 7

prob.selected.is.max(n, h2z, h2ps) # problem when h2ps is small 
prob.selected.is.max(n, h2z, h2ps, 'simulation')

