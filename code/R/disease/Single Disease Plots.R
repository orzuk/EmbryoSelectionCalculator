setwd("C:/Users/Adam/OneDrive/Documents/Embryo Sims/functions_and_simulations")
source("score_generation_functions.R")
source("multiple_selection_function.R")
source("gain_disease.R")
source("decile_exclusion_selection.R")
source("trunc_selection_function.R")
library("MASS")
library("ggplot2")
library("lqmm")
library("dbplyr")
library("data.table")
library("mvtnorm")
library("beepr")
library("cowplot")


## SCZ: 
k_SCZ = 0.01
h2_pgs_SCZ = 0.1
## MDD: 
k_MDD = 0.15
h2_pgs_MDD = 0.02
## IBD: 
k_IBD = 0.005
h2_pgs_IBD = 0.131


### 1. LOWEST RISK: simulation and numerical integration


max_N <- 50
low_risk_data <- expand.grid(N_embryos = seq(1, max_N, 1), pgs = c(0.02, 0.1, 0.131), random_average_disease = 0, selected_average_disease = 0)
low_risk_data$K <- c(rep(0.15, max_N), rep(0.01, max_N), rep(0.005, max_N))
low_risk_data$name <- c(rep("MDD", max_N), rep("SCZ", max_N), rep("IBD", max_N))
cor <- matrix(c(1))
nsims1 <- 50000
N_hi <- max(low_risk_data$N_embryos)
M <- 1 # number of diseases

## Running the simulation
for(a in seq(N_hi, nrow(low_risk_data), by = N_hi)) {
  
  selection_stats <- as.data.frame(matrix(0, nrow = N_hi, ncol = 3))
  
  # the total D_sums that are eventually used to calculate disease rates
  colnames(selection_stats) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")
  pgs <- low_risk_data$pgs[a] #can set this to be a vector of different values
  h2 <- pgs + 0.4
  N_embryos <- low_risk_data$N_embryos[a]
  matrices <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)  
  Sigma <- matrices$Sigma # the PGS covariance matrix
  Sigma_G <- matrices$Sigma_G # The cov matrix for the genetic variance not captured by the PGS
  E_Matrix <- matrices$E_Matrix # The remaining variance covariance matrix
  Kin_Matrix <- matrices$Kin_Matrix # The N x N relatedness matrix; 1's on the diag, 0.5 on the off-diag
  E <- Sigma_G + E_Matrix # Error matrix is the sum of E_G and E_E
  threshold <- qnorm(1-low_risk_data$K[a])
  system.time(
    for(i in 1:nsims1) {
      score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, N_hi)
      PGS_scores <- score_matrices$PGS_scores
      D_scores <- score_matrices$D_scores
      pass_to_data_frame <- speedy_selection(PGS_scores, D_scores, E, threshold, N_hi)
      selection_stats[1:N_hi, 1:3] <- selection_stats[1:N_hi, 1:3] + pass_to_data_frame[1:N_hi, 1:3]
    }
  )
  colnames(selection_stats) <- c("Embryo Number", "Selected_D_sum", "Random_D_sum")
  index <- 1
  for (p in seq(a+1-N_hi, a, 1)) {
    low_risk_data$selected_average_disease[p] <- selection_stats$Selected_D_sum[index]/nsims1
    low_risk_data$random_average_disease[p] <- selection_stats$Random_D_sum[index]/nsims1
    index <- index + 1
  }
}
beep("coin")
low_risk_data$average_absolute_disease_reduction <- low_risk_data$random_average_disease - low_risk_data$selected_average_disease
low_risk_data$average_relative_disease_reduction <- low_risk_data$average_absolute_disease_reduction/low_risk_data$random_average_disease

low_risk_data$ARR_expected <- low_risk_data$K - low_risk_data$selected_average_disease
low_risk_data$RRR_expected <- low_risk_data$ARR_expected/low_risk_data$K

low_risk_data$gain_d <- 0

#### Add the numerical integration data
for(i in 1:nrow(low_risk_data)) {
  low_risk_data$integration_RRR[i] <- gain.disease(low_risk_data$N_embryos[i], low_risk_data$K[i], low_risk_data$pgs[i]) #plots expected gain from numerical integration approach
  low_risk_data$integration_ARR[i] <- low_risk_data$K[i]*low_risk_data$integration_RRR[i]
}



#### 2. DECILE EXCLUSION (limited by N): simulation and numerical integration

nsims2 <- 100000

max_N <- 50
decile_exclusion_data <- expand.grid(N_embryos = seq(1, max_N, 1), pgs = c(0.02, 0.1, 0.131), random_average_disease = 0, selected_average_disease = 0)
decile_exclusion_data$K <- c(rep(0.15, max_N), rep(0.01, max_N), rep(0.005, max_N))
decile_exclusion_data$name <- c(rep("MDD", max_N), rep("SCZ", max_N), rep("IBD", max_N))
cor <- matrix(c(1))
N_hi <- max(decile_exclusion_data$N_embryos)
M <- 1

## Running the simulation

for(a in seq(N_hi, nrow(decile_exclusion_data), by = N_hi)) {
  
  selection_stats <- as.data.frame(matrix(0, nrow = N_hi, ncol = 3))
  
  # the total D_sums that are eventually used to calculate disease rates
  colnames(selection_stats) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")
  pgs <- decile_exclusion_data$pgs[a] #can set this to be a vector of different values
  h2 <- pgs + 0.4
  N_embryos <- decile_exclusion_data$N_embryos[a]
  matrices <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)  
  Sigma <- matrices$Sigma # the PGS covariance matrix
  Sigma_G <- matrices$Sigma_G # The cov matrix for the genetic variance not captured by the PGS
  E_Matrix <- matrices$E_Matrix # The remaining variance covariance matrix
  Kin_Matrix <- matrices$Kin_Matrix # The N x N relatedness matrix; 1's on the diag, 0.5 on the off-diag
  E <- Sigma_G + E_Matrix # Error matrix is the sum of E_G and E_E
  threshold <- qnorm(1-decile_exclusion_data$K[a])
  for(i in 1:nsims2) {
    score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, N_hi)
    PGS_scores <- score_matrices$PGS_scores
    Z_scores <- score_matrices$Z_scores
    D_scores <- score_matrices$D_scores
    pass_to_data_frame <- decile_exclusion_selection(pgs, PGS_scores, D_scores, E, threshold, N_hi)
    selection_stats[1:N_hi, 1:3] <- selection_stats[1:N_hi, 1:3] + pass_to_data_frame[1:N_hi, 1:3]
  }
  colnames(selection_stats) <- c("Embryo Number", "Selected_D_sum", "Random_D_sum")
  index <- 1
  for (p in seq(a+1-N_hi, a, 1)) {
    decile_exclusion_data$selected_average_disease[p] <- selection_stats$Selected_D_sum[index]/nsims2
    decile_exclusion_data$random_average_disease[p] <- selection_stats$Random_D_sum[index]/nsims2
    index <- index + 1
  }
}

decile_exclusion_data$average_absolute_disease_reduction <- decile_exclusion_data$random_average_disease - decile_exclusion_data$selected_average_disease
decile_exclusion_data$average_relative_disease_reduction <- decile_exclusion_data$average_absolute_disease_reduction/decile_exclusion_data$random_average_disease

decile_exclusion_data$ARR_expected <- decile_exclusion_data$K - decile_exclusion_data$selected_average_disease
decile_exclusion_data$RRR_expected <- decile_exclusion_data$ARR_expected/decile_exclusion_data$K
beep("coin")

## Numerical Integration (Eqs. 3.9-3.10)
### still need to figure this out


## 3. DECILE EXCLUSION NO EMBRYO LIMIT: 

#### numerical integration data (relative risk reduction)

## MDD
k <- 0.15; r2 <- 0.02; q <- 0.1; zq_ps <- qnorm(q, 0, sqrt(r2))
MDD_expected_gain <- 1-((integrate(trunc_selection_function, zq_ps, Inf, q, k, r2)$value)/k) # expected relative gain


## SCZ
k <- 0.01; r2 <- 0.1; q <- 0.1; zq_ps <- qnorm(q, 0, sqrt(r2))
SCZ_expected_gain <- 1-((integrate(trunc_selection_function, zq_ps, Inf, q, k, r2)$value)/k) # expected relative gain


## IBD
k <- 0.005; r2 <- 0.131; q <- 0.1; zq_ps <- qnorm(q, 0, sqrt(r2))
IBD_expected_gain <- 1-((integrate(trunc_selection_function, zq_ps, Inf, q, k, r2)$value)/k) # expected relative gain




### PLOTS ###

## Some of the plots will throw warnings due to text parsing in the annotations. 
## Because of this, to save a graph, you'll have to plot it and then save it.


### 1a. Lowest Risk, Absolute Gain
subset_low_risk_MDD = subset(low_risk_data, name == "MDD")
subset_low_risk_SCZ = subset(low_risk_data, name == "SCZ")
subset_low_risk_IBD = subset(low_risk_data, name == "IBD")

anno_MDD_1 <- bquote("Disease Prevalence:" ~ .(k_MDD))
anno_MDD_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_MDD))
anno_MDD_3 <- bquote("Number of Simulations: " ~ .(nsims1))
MDD_ARR_1 <- ggplot(data = subset_low_risk_MDD, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "red", "Numerical Integration" = "blue")) +
  geom_line(aes(y = integration_ARR,  color = 'Numerical Integration')) +
  geom_point(aes(y = ARR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Absolute Disease Gain") +
  labs(title = "Major Depression: Absolute Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.021, label = anno_MDD_1) +
  annotate("text", x = 35, y = 0.018, label = anno_MDD_2) +
  annotate("text", x = 35, y = 0.015, label = anno_MDD_3)
                          

anno_SCZ_1 <- bquote("Disease Prevalence:" ~ .(k_SCZ))
anno_SCZ_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_SCZ))
anno_SCZ_3 <- bquote("Number of Simulations: " ~ .(nsims1))
SCZ_ARR_1 <- ggplot(data = subset_low_risk_SCZ, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "red", "Numerical Integration" = "blue")) +
  geom_line(aes(y = integration_ARR,  color = 'Numerical Integration')) +
  geom_point(aes(y = ARR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Absolute Disease Gain") +
  labs(title = "Schizophrenia: Absolute Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.004, label = anno_SCZ_1) +
  annotate("text", x = 35, y = 0.0035, label = anno_SCZ_2) +
  annotate("text", x = 35, y = 0.003, label = anno_SCZ_3)


anno_IBD_1 <- bquote("Disease Prevalence:" ~ .(k_IBD))
anno_IBD_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_IBD))
anno_IBD_3 <- bquote("Number of Simulations: " ~ .(nsims1))
IBD_ARR_1 <- ggplot(data = subset_low_risk_IBD, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "red", "Numerical Integration" = "blue")) +
  geom_line(aes(y = integration_ARR,  color = 'Numerical Integration')) +
  geom_point(aes(y = ARR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Absolute Disease Gain") +
  labs(title = "Inflammatory Bowel Disease: Absolute Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.00225, label = anno_IBD_1) +
  annotate("text", x = 35, y = 0.002, label = anno_IBD_2) +
  annotate("text", x = 35, y = 0.00175, label = anno_IBD_3)




### 1b. Lowest Risk, Relative Gain 

anno_MDD_1 <- bquote("Disease Prevalence:" ~ .(k_MDD))
anno_MDD_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_MDD))
anno_MDD_3 <- bquote("Number of Simulations: " ~ .(nsims1))
MDD_RRR_1 <- ggplot(data = subset_low_risk_MDD, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "purple", "Numerical Integration" = "orange")) +
  geom_line(aes(y = gain_d,  color = 'Numerical Integration')) +
  geom_point(aes(y = RRR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Relative Disease Gain") +
  labs(title = "Major Depression: Relative Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.16, label = anno_MDD_1) +
  annotate("text", x = 35, y = 0.14, label = anno_MDD_2) +
  annotate("text", x = 35, y = 0.12, label = anno_MDD_3)


anno_SCZ_1 <- bquote("Disease Prevalence:" ~ .(k_SCZ))
anno_SCZ_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_SCZ))
anno_SCZ_3 <- bquote("Number of Simulations: " ~ .(nsims1))
SCZ_RRR_1 <- ggplot(data = subset_low_risk_SCZ, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "purple", "Numerical Integration" = "orange")) +
  geom_line(aes(y = gain_d,  color = 'Numerical Integration')) +
  geom_point(aes(y = RRR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Relative Disease Gain") +
  labs(title = "Schizophrenia: Relative Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.32, label = anno_SCZ_1) +
  annotate("text", x = 35, y = 0.28, label = anno_SCZ_2) +
  annotate("text", x = 35, y = 0.24, label = anno_SCZ_3)


anno_IBD_1 <- bquote("Disease Prevalence:" ~ .(k_IBD))
anno_IBD_2 <- bquote(h[pgs]^2 ~ "=" ~ .(h2_pgs_IBD))
anno_IBD_3 <- bquote("Number of Simulations: " ~ .(nsims1))
IBD_RRR_1 <- ggplot(data = subset_low_risk_IBD, aes(x = N_embryos)) + 
  scale_color_manual(name = "Gain Calculation Method", values = c("Simulation" = "purple", "Numerical Integration" = "orange")) +
  geom_line(aes(y = gain_d,  color = 'Numerical Integration')) +
  geom_point(aes(y = RRR_expected, color = 'Simulation')) +
  xlab("Number of Embryos") + ylab("Relative Disease Gain") +
  labs(title = "Inflammatory Bowel Disease: Relative Disease Gain as Function of Number of Embryos", subtitle = "Selection Method: Lowest Risk") +
  annotate("text", x = 35, y = 0.42, label = anno_IBD_1) +
  annotate("text", x = 35, y = 0.37, label = anno_IBD_2) +
  annotate("text", x = 35, y = 0.32, label = anno_IBD_3)



### 2a. Decile exclusion (limited by N), absolute gain

subset_decile_exclusion_MDD = subset(decile_exclusion_data, name == "MDD")
subset_decile_exclusion_SCZ = subset(decile_exclusion_data, name == "SCZ")
subset_decile_exclusion_IBD = subset(decile_exclusion_data, name == "IBD")


#### 2b. Decile exclusion (limited by N), relative gain



#### 3a. Decile exclusion no limit, absolute gain



#### 3b. Decile exclusion no limit, relative gain




