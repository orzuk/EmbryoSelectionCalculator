## Decile Exclusion No Limit
# Implemented in Equations 3.5 and 3.6 (here I'm using Shai's formulation)


# q is the cutoff percentage (the quantile) for the PS score distribution
# k is the cutoff percentage for exclusion of embryos
# qnorm(q) gives you the z score

# for the situation where a z score _less_ than z_k results in disease
# the PS cutoff progresses rightward as q increases
# the error integral represents all error values that, when added to the PS, result in a z score less than z_k; hence, D=1.

trunc_selection_function = function(x,q,k,r2)
{
  zk = qnorm(k) #kth quantile of the standard normal -- i.e., the liability score cutoff
  zq_ps = qnorm(q,mean=0,sd=sqrt(r2)) # the cutoff z score for the PGS
  y = dnorm(x,sd=sqrt(r2)) / pnorm(zq_ps,sd=sqrt(r2),lower.tail = F) # P(X=x)/P(X>zq_ps) - i.e., f_ps
  y = y * pnorm(zk-x,sd=sqrt(1-r2)) # multiply by the error term integral, which the P(error term < z_k-x); i.e., the prob that the individuals total z < z_k 
  return(y)
}


### Loop to Run to check function is working correctly

# # Sample Parameters
# r2 = 0.1 # h^2_{pgs}
# qs = seq(0, 0.3, 0.01) # the percentile of PS to exclude
# k = 0.05 # disease prevalence
# risk = numeric(length(qs)) # the average risk of non-excluded embryos
# gain = numeric(length(qs))


# for (i in seq_along(qs))
# {
#   q = qs[i]
#   zq_ps = qnorm(q,mean=0,sd=sqrt(r2)) 
#   risk[i] = integrate(trunc_selection_function,zq_ps,Inf,q,k,r2)$value
#   gain[i] <- 1 - risk[i]/k
# }