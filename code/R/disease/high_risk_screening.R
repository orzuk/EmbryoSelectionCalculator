h2_ps <- 0.3
h_ps <- sqrt(h2_ps)
k <- 0.05
N <- seq(1, 15, 1)
gain <- numeric(length = length(N))
gamma <- qnorm(k, 0, h_ps)

f_x_function <- function(s, t, n, h_ps) {
  (pnorm(s)^(n-1))*dnorm(s)*dnorm((sqrt(2)*t/h_ps)-s)
}

f_x_max <- function(s, t, n, h_ps) {
  return((sqrt(2)*n/h_ps)*integrate(f_x_function, -Inf, Inf, t, n, h_ps, abs.tol = 1e-4)$value) 
}

first_integrand<- function(t, n, h_ps, k) {
  y <- pnorm((qnorm(k)-t)/sqrt(1-(h_ps)^2)) 
  y <- f_x_max(s, t, n, h_ps) * y
}
  
second_integrand <- function(t, n, h_ps, k) {
  dnorm(t)*pnorm((qnorm(k)-t*h_ps)/sqrt(1-(h_ps)^2))
}
  
  
for(i in seq_along(N)) {
  n <- N[i]
  lower_bound <- gamma/h_ps
  risk <- ((pnorm(gamma)^n*integrate(first_integrand, -Inf, gamma, n, h_ps, k)$value)/
          f_x_max(s, gamma, n, h_ps)) + # <-- this is where I have a question
          ((1-pnorm(gamma)^n)*integrate(second_integrand, lower_bound, Inf, n, h_ps, k)$value)/(1-pnorm(lower_bound))      
  gain[i] <- 1-(risk/k)
}

plot(N, gain)


