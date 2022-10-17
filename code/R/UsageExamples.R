# Run examples: 
source('gain_moments.R')
# Usage example 
T <- 4 # number of traits 
w <- c(0.1, 0.2, 0.3, 0.4)
Sigma <- matrix(0, T, T)
for(i in c(1:T))
  Sigma[i,i] <- 1
Sigma[1,2] <- Sigma[1,3] <- Sigma[2,1] <- Sigma[3,1] <- 0.5 

isposdef(Sigma) # check

G <- multi_trait_gain_mean(10, Sigma, w, 'exact')


select_vars <- function(x)
{
  return( x[which.min( (x - 1)**2 )])
}
  

