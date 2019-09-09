library(mvtnorm) # for multidimensional Gaussians
library(MASS)
source('truncated_sum_norm_cdf.R')

# Probability that selected embryo has phenotype smaller than the average phenotype of all siblings
# n - # of embryos 
# h2z - heritability of trait
# h2ps - heritability explained by polygenic score PS
# method.str - how to compute (numeric or simulations)
#
prob.gain.negative <- function(n, h2z, h2ps, method.str='numeric')
{	
  if(method.str == 'numeric')
  {
    sigma <- sqrt(h2ps/(2 - h2z - h2ps))
    if(h2ps < 0.0000000001) # avoid numerical problems 
      return(1/n)
    return( integrate(prob.negative.integrand, -20*sigma, 20*sigma, n=n, h2z=h2z, h2ps=h2ps)$value)
  }
  else # simulations 
  {
    iters <- 500000
    Sigma.PS <- 0.5 * h2ps * (matrix(1, n, n) + diag(n))
    Sigma.eps <- 0.5 * ((h2z-h2ps) * matrix(1, n, n) + (2-h2ps-h2z) * diag(n))
    PS <- mvrnorm(n = iters, mu=rep(0,n), Sigma=Sigma.PS)
    eps <- mvrnorm(n = iters, mu=rep(0,n), Sigma=Sigma.eps)
    z <- PS + eps
    z.selected <- z[(max.col(PS)-1) * iters + c(1:iters)]
    return( mean(z.selected < rowMeans(z)) ) # find indices maximizing and compare
  }
}

# 
prob.negative.integrand <- function(y, n, h2z, h2ps)
{
  sigma <- sqrt(h2ps/(2 - h2z - h2ps))
  # Use truncated Gaussians approximation
  x <- sqrt(n-1)*(y + (sigma) * dnorm(y/sigma) / pnorm(y/sigma)) / 
    sqrt( n + sigma^2 * (1-  (y/sigma) * dnorm(y/sigma) / pnorm(y/sigma) - dnorm(y/sigma)^2 / pnorm(y/sigma)^2)  )
  return( (n/sigma) * pnorm(y/sigma)^(n-1) * dnorm(y/sigma) * (1 - pnorm(x) ) )
} 


# Example (real data parameters)
h2ps <- 0.27
h2z <- 0.8
n <- 7

prob.gain.negative(n, h2z, h2ps, 'simulation')
prob.gain.negative(n, h2z, h2ps) # numeric approximation 


