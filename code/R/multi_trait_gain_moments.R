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

# Usage example 
T <- 4 # number of traits 
w <- c(0.1, 0.2, 0.3, 0.4)
Sigma <- matrix(0, T, T)
for(i in c(1:T))
  Sigma[i,i] <- 1
Sigma[1,2] <- Sigma[1,3] <- Sigma[2,1] <- Sigma[3,1] <- 0.5 

isposdef(Sigma) # check

G <- multi.trait.gain.mean(10, Sigma, w, 'exact')


