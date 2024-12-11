source("score_generation_functions.R")
library(MASS)
library(mvtnorm)
library(ggplot2)
library(lqmm) # for making matrices positive definite so mvrnorm works
library(tmvtnorm) # For the MVN integration approach

### RISK PREDICTION USING MONTE CARLO SIMULATION, for general multivariate case with general loss #####

# Use simulation to predict the risk of disease 1 through M (R_1....R_M) for
# an embryo with polygenic score vector, Z. 


# Input
# Z = {z_1, z_2, ..., z_M}, PGS risk scores
# E = E_E + E_G used for ~MVN(Z, E).
# Z, which functions as the mean, mu. 

# simulate n_sims times

# Output:
# the simulated risk (R_1,..., R_M) and probability of having no diseases at all R_0.

simulate_PGS_disease_risk <- function(Z, K, E, n_sims = 1000) {
  threshold <- qnorm(1-K)
  samples <- Z + mvrnorm(n_sims, mu = rep(0, length(K)), Sigma = E) # generate a ton of samples with mean Z
  D_scores <- return_D_scores(samples, threshold) # binary scoring of each embryo score vector
  R_0 <- sum(rowSums(D_scores) == 0)/n_sims # return the percentage of embryos which have no diseases
  R <- as.matrix(colSums(D_scores)/n_sims) # vector of probabilities for disease 1 thru M
  R_list <- list() # concatenate the two lists and return
  R_list$R_0 <- R_0
  R_list$R <- R
  return(R_list)
}


########### EXAMPLES ############

#### EXAMPLE 1: M = 5 CASE #####
# cor <- matrix(c(1.000, 0.3, 0.4, 0.5, 0.6,
#                 0.3, 1.000, 0.5, 0.3, 0.2,
#                 0.4, 0.5, 1.000, 0.400, 0.60,
#                 0.5, 0.3, 0.400, 1.000, 0.40,
#                 0.6, 0.2, 0.600, 0.400, 1.00), nrow = 5, ncol = 5)

cor <- matrix(0.1, nrow = 5, ncol = 5)
diag(cor) <- 1
M <- nrow(cor)


### EXAMPLE: 2 traits
cor <- matrix(0.1, nrow = 2, ncol = 2)
diag(cor) <- 1
M <- nrow(cor)

### SIMULATION WITH MULTIPLE PARAMETERS ###

K <- c(1:10 %o% 10^(-2:-3)); K <- K[1:(length(K)-1)] # vector of K values to simulate with

Z_frame <- expand.grid(z = seq(0, 3, 0.1), K = K, PGS = seq(0.05, 0.9, 0.05), R_0 = 0, R_1 = 0, R_2 = 0, R_3 = 0, R_4 = 0, R_5 = 0, MV_CDF_R0 = 0)
for(i in 1:nrow(Z_frame)) {
  K <- rep(Z_frame$K[i], M)
  threshold <- qnorm(1-K)
  Z <- rep(Z_frame$z[i], M)
  pgs <- rep(Z_frame$PGS[i], M)
  h2 <- pgs + runif(M, 0.01, 0.09)
  matrices <- generate_matrices_for_sim(cor, pgs, h2, 10)  
  E <- matrices$Sigma_G + matrices$E_Matrix
  test_sim_risk <- simulate_PGS_disease_risk(Z, K, E, n_sims = 1000)
  Z_frame$R_0[i] <- test_sim_risk$R_0
  R <- test_sim_risk$R
  Z_frame$R_1[i] <- R[1]
  Z_frame$R_2[i] <- R[2]
  Z_frame$R_3[i] <- R[3]
  Z_frame$R_4[i] <- R[4]
  Z_frame$R_5[i] <- R[5]
  Z_frame$MV_CDF_R0[i] <- pmvnorm(lower = -Inf, upper = threshold, mean = Z, sigma = E)
    
}

# Plots where R_0 is calculated via integration using pmvnorm

ggplot(data = Z_frame, aes(x = z, y = R_1)) + geom_point(aes(color = PGS)) + facet_wrap(~K, nrow =3) +
  ggtitle("P(Have Disease #1) as a Function of Z, K, and h^2_{pgs}") + xlab("Z-score of Embryo") + 
  ylab("P(Have Disease #1") + scale_color_gradient(name = "h^2_{pgs}")


# Simulation-based plots

ggplot(data = Z_frame, aes(x = z, y = R_1)) + geom_point(aes(color = PGS)) + facet_wrap(~K, nrow =3) +
  ggtitle("P(Have Disease #1) as a Function of Z, K, and h^2_{pgs}") + xlab("Z-score of Embryo") + 
  ylab("P(Have Disease #1") + scale_color_gradient(name = "h^2_{pgs}")


ggplot(data = Z_frame, aes(x = z, y = MV_CDF_R0)) + geom_point(aes(color = PGS)) + facet_wrap(~K, nrow =3) + 
  ggtitle("P(Disease Free) as a Function of Z, K, and h^2_{pgs} \n Number of Traits: 2") + xlab("Z-score of Embryo") + 
  ylab("P(Disease Free)") + scale_color_gradient(name = "h^2_{pgs}")

# Take a subset of the data:
subset_0.1 <- subset(Z_frame, K == 0.02)
ggplot(data = subset_0.1, aes(x = z, y = R_0, colour = PGS)) + geom_point() +  ggtitle("P(Disease Free) as a Function of Z and h^2_{pgs}, K = 0.1") + xlab("Z-score of Embryo") + 
  ylab("P(Disease Free)") + scale_color_gradient(name = "h^2_{pgs}")


### Determining the Minimum Effective Simulation Size
## Using Point Prediction Examples

nsims <- seq(100, 10000, 100)
sim_size <- matrix(ncol = 2, nrow = length(nsims))
sim_size[,1] <- nsims

Z <- rnorm(5)
K <- runif(5, 0.01, 0.1)
for ( i in 1:nrow(sim_size)) {
  nsims <- sim_size[i,1]
  risk_variance <- replicate(100, simulate_PGS_disease_risk(Z, K, E, n_sims = nsims)[[1]])
  sim_size[i, 2] <- var(risk_variance)
}

sim_size <- as.data.frame(sim_size)
ggplot(data = sim_size, aes(x = sim_size[,1], y = sim_size[,2])) + geom_point() + xlab("Simulation Size") + ylab("Variance") + ggtitle("Variance of Risk Prediction Estimate as Function of Simulation Size")


### We hit diminishing returns to sample size at around 2500. 

