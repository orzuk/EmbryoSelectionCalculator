# Conditional cumulative probability of Y+epsilon < x given that Y < y
# Alternative: here eps ~ N(0, 1) ; y ~ N(0, sigma)
# h2z - heritability of trait
# h2ps - heritability explained by polygenic score PS
#
truncated.sum.norm.cdf <- function(x, h2z, h2ps, y)
{	
  sigma <- sqrt(h2ps/(2 - h2z - h2ps))
  return( integrate(dnorm.prod, -20*sigma, y, y=y, x=x, sigma=sigma )$value ) # integral starts from 20 instead of -infinity for numerical reasons 
}

# Internal function: Density of y+eps conditioned on <y
dnorm.prod <- function(t, y, x, sigma)
{
  return( dnorm(t/sigma) * pnorm(x-t) / (sigma*pnorm(y/sigma)) )
} 

# Usage example: 
truncated.sum.norm.cdf(1, 0.8, 0.3, 200)


