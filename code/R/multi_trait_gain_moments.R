# Mean gap between selected embryo (top PS) and best embryo (top phenotype z)
# n - # of embryos
# Sigma - variance-covariance matrix of traits
# w - weight vector for each trait
# method.str - approx or exact
#
multi.trait.gain.mean <- function(n, Sigma, w, method.str)
{	
  return ( gain.moments(n, t(w) %*% Sigma %*% w / 2, method.str)$E.G )
}

