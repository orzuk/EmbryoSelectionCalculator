gain_conditional_parent_pgs <- function(mean_parental_pgs, N, h2_pgs, k, ..) {
  return(1 - (integrate(selected_conditional_prevalence, -Inf, Inf, mean_parental_pgs = mean_parental_pgs, 
                        N = N, h2_pgs = h2_pgs, k = k)$value /
                integrate(base_prevalence, -Inf, Inf, mean_parental_pgs = mean_parental_pgs, h2_pgs = h2_pgs, k = k)$value))
}

selected_conditional_prevalence <- function(mean_parental_pgs, N, h2_pgs, K, t) {
#  z_K <- qnorm(1-K)
  z_K <- qnorm(K)
  return((N*pnorm(t, mean_parental_pgs, sqrt(h2_pgs/2))^(N-1))*
           dnorm(t, mean_parental_pgs, sqrt(h2_pgs/2))*pnorm((z_K-t)/sqrt(1-h2_pgs)))
}

base_prevalence <- function(mean_parental_pgs, h2_pgs, k, t) {
  z_k <- qnorm(1-k)
  return(dnorm(t, mean_parental_pgs, sqrt(h2_pgs/2))*pnorm((z_k-t)/sqrt(1-h2_pgs)))
}

# numerator --> P(disease_{selected embryo} | parental pgs)
integrate(selected_conditional_prevalence, -Inf, Inf, mean_parental_pgs = 0, N = 10, h2_pgs = 0.131, K = 0.1)$value

