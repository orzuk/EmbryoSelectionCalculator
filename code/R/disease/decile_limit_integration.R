sigma_z <- 1
h2_ps <- 0.02
h_ps <- sqrt(h2_ps)
k <- 0.15 # prevalence
N <- seq(1, 20, 1) # enough to show first few 
gain <- numeric(length = length(N))
q_exclude <- 0.5 # 0.1 is decile
gamma <- qnorm(q_exclude, 0, h_ps) # threshold for exclusion. Doesn't have to be the same k

# Used inside integral of f_x_max
f_x_function <- function(s, t, n, h_ps, sigma_z) {
  return( pnorm(sqrt(2)*s/(sigma_z*h_ps))^(n-1) * 
            dnorm((2*s-t)/(sigma_z*h_ps)) )
}

# Density of max(Y_1,..y_n) + Z
f_x_max <- function(t, n, h_ps, sigma_z) {
  sigma_z_2 <- sigma_z^2
  a <- (2*n/(sigma_z_2*h_ps)^2)*dnorm(t/(1*sigma_z*h_ps))
  for(i in seq_along(t)) # cannot vectorize because t is inside the intergral 
  {
    a[i] <- a[i] * integrate(f_x_function, -Inf, Inf,  ### This is where the integral diverges when limits of integration are -Inf, Inf
                             t[i], n, h_ps, sigma_z, rel.tol = 10e-4)$value
  }
  return(a)
}

F_max_integrand <- function(s, t, n, h_ps, sigma_z) {
  return((pnorm((sqrt(2)*(t-s))/(sigma_z*h_ps))^n)*dnorm(sqrt(2)*(s)/(sigma_z*h_ps)))
}

# Cumulative of max(Y_1,..y_n) + Z
F_x_max <- function(t, n, h_ps, sigma_z) {
  a <- sqrt(2) / (sigma_z*h_ps) # no multiplication by n 
  return(a*integrate(F_max_integrand, -Inf, Inf, t, n, h_ps, sigma_z)$value)
}

first_integrand<- function(t, n, h_ps, sigma_z, k) {
  y <- pnorm((qnorm(k)-t)/sqrt(1-h_ps^2)) 
  y <- f_x_max(t, n, h_ps, sigma_z) * y  # first argument is t? 
  return(y)
}

second_integrand <- function(t, h_ps, k) {
  dnorm(t)*pnorm((qnorm(k)-t*h_ps)/sqrt(1-h_ps^2))
}

# New function for disease risk 
disease_risk <- function(n, h_ps, sigma_z, k, q_exclude=0.1, limit=1)   {
  gamma <- qnorm(q_exclude, 0, h_ps) # not k
  lower_bound <- gamma/h_ps  
  risk <- limit * (pnorm(gamma/h_ps)^n*integrate(first_integrand, -Inf, gamma, n, h_ps, sigma_z, k)$value /
                     F_x_max(gamma, n, h_ps, sigma_z)) + 
    (1-pnorm(gamma/h_ps)^n)^(limit) * 
    integrate(second_integrand, lower_bound, Inf, h_ps, k)$value / (1-pnorm(lower_bound))  
  return(risk)
}

gain_no_limit <- gain  
for(i in seq_along(N)) {
  gain[i] <- 1-(disease_risk(N[i], h_ps, sigma_z, k, q_exclude)/k)
  gain_no_limit[i] <- 1-(disease_risk(N[i], h_ps, sigma_z, k, q_exclude, 0)/k)
}

plot(N, gain)
lines(N, gain_no_limit, col='red')

