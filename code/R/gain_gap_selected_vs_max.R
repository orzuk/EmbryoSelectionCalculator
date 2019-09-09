# Mean gap between selected embryo (top PS) and best embryo (top phenotype z)
# n - # of embryos
# sigma.z - variance of trait
# h2z - heritability of trait
# h2ps - heritability explained by polygenic score PS
# method.str - approx or exact
#
gain.gap.selected.vs.max.embryo <- function(n, sigma_z, h2z, h2ps, method.str)
{	
	return ( gain.moments(n, 1, method.str)$E.G * sigma_z^2 * ( sqrt(1-h2z/2) - sqrt(h2ps) ) )
}


# Usage example: 
gain.gap.selected.vs.max.embryo(1000, 1, 0.8, 0.3, 'exact')
  