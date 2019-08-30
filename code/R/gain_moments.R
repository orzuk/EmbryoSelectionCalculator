# Compute first two moments of gain (for quantitative traits)
# n - # of embryos
# subtract.mean - subtract mean of x_i
# method.str - approx or exact
# sib.flag - include correlations between siblings 
#
gain.moments <- function(n, subtract.mean=1, method.str, sib.flag=1)
{
  switch (method.str,
          'approx'={
            E.G = 0.5^sib.flag * (qnorm(1-1/n) - digamma(1) / (n * dnorm( qnorm(1/n))))        
            a.n = n * dnorm( qnorm(1/n)) # Approximation from paper (not good for small n) 
            Var.G = 0.5^sib.flag * pi^2/6 / a.n^2
          },
          'exact'={
            E.G <- Var.G <- numeric(length(n))
            for(i in 1:length(n)) # vectorize
            {
              E.G[i] = (n[i]/ (sqrt(2)^sib.flag)) * integrate( gauss.poly, -10, 10, n=c(n[i]-1,1,1))$value
              Var.G[i] = 0.5^sib.flag * n[i] * integrate(gauss.poly, -10, 10, n=c(n[i]-1,1,2))$value - E.G[i]^2
            }
          }
  ) # end switch 
  
  if(subtract.mean) # correct variance 
    Var.G = Var.G - 1/(2^sib.flag*n)
  else
    Var.G = Var.G + 0.5*sib.flag # add constant contribution
  return(list(E.G=E.G, Var.G=Var.G))
}  

# Helper functios for integration:
# Phi(t)^n[1] * phi(t)^n[2] * t^n[3]  
gauss.poly <- function(t, n)
{
  return (pnorm(t)^n[1]*dnorm(t)^n[2]*t^n[3])
}
