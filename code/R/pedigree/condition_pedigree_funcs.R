# Compute conditional gain given the parents scores
gain_conditional_pedigreet_pgs <- function(pedigree_G, pedigree_ps, N, h2_pgs, k, ..) {

    return(1 - (integrate(selected_conditional_prevalence, -Inf, Inf, pedigree_ps = mean_parental_pgs, 
                        N = N, h2_pgs = h2_pgs, k = k)$value /
                integrate(base_prevalence, -Inf, Inf, mean_parental_pgs = mean_parental_pgs, h2_pgs = h2_pgs, k = k)$value))
}


gain_conditional_pedigreet_pgs_and_trait <- function(pedigree_G, pedigree_ps, pedigree_y, N, h2_pgs, k, ..) {
  # T.B.D. implementation
  return(0)
  }


# Usage example: set pedigree with values, and then compute the embryos scores and traits distribution: 


