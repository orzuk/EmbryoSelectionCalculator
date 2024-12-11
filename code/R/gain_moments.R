# Compute first two moments of gain (for quantitative traits)
# n - # of embryos
# h2ps - heritability explained by the polygenic score (We assume Var(Z)=1. If not, then this should be: h2ps*Var(Z))
# method.str - approx or exact
# subtract.mean - subtract mean of x_i. If zero, we compute moments of the top score
# sib.flag - include correlations between siblings 
#
gain_moments <- function(n, h2ps=1, method.str, subtract.mean=1, sib.flag=1)
{
  if(!exists("method.str"))
    method.str <- 'approx'
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
              E.G[i] = (n[i]/ (sqrt(2)^sib.flag)) * integrate( gauss_poly, -10, 10, n=c(n[i]-1,1,1))$value
              Var.G[i] = 0.5^sib.flag * n[i] * integrate(gauss_poly, -10, 10, n=c(n[i]-1,1,2))$value - E.G[i]^2
            }
          }
  ) # end switch 
  
  if(subtract.mean) # correct variance 
    Var.G = Var.G - 1/(2^sib.flag*n)
  else
    Var.G = Var.G + 0.5*sib.flag # add constant contribution
  return(list(E.G=E.G*sqrt(h2ps), Var.G=Var.G*h2ps)) # correct for variance explained
}  

# Helper functios for integration:
# Phi(t)^n[1] * phi(t)^n[2] * t^n[3]  
gauss_poly <- function(t, n)
{
  return (pnorm(t)^n[1]*dnorm(t)^n[2]*t^n[3])
}



# Mean gap between selected embryo (top PS) and best embryo (top phenotype z)
# n - # of embryos
# Sigma.T - variance-covariance matrix of traits
# w - weight vector for each trait
# method.str - approx or exact
#
multi_trait_gain_mean <- function(n, Sigma.T, w, method.str)
{	
  if(!exists("method.str"))
    method.str <- 'approx'
  return ( gain_moments(n, t(w) %*% Sigma.T %*% w / 2, method.str)$E.G )
}


# Mean gap between embryo from selected gametes (top PS) and random embryos (top phenotype z)
# n - # of embryos
# Sigma.T - variance-covariance matrix of traits
# w - weight vector for each trait
# method.str - approx or exact
#
multi_trait_gamete_gain_mean <- function(n, Sigma.T, w, method.str)
{	
  if(!exists("method.str"))
    method.str <- 'approx'
  return ( gain_moments(n, t(w) %*% Sigma.T %*% w / 2, method.str)$E.G )  # Change formula here 
}
