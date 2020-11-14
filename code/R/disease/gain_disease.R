# Compute relative increase in disease probability
# n - # of embryos
# K - disease prevalence
# h2 - variance explained by polygenic score (on liability scale)
gain.disease <- function(n, K, h2)
{
  return(  1 - (1/K) * integrate( gain.d, -10, 10, K=K, h2=h2, n=n )$value )
}

# integrand
gain.d <- function(t, K, h2, n)
{  
  return( dnorm(t) * pnorm( (qnorm(K)-t*sqrt(1-h2/2))*sqrt(2/h2)  )^n )
}

# Example
gain.disease(5, 0.01, 0.1)

