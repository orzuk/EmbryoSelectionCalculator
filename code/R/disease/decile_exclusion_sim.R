library(MASS)
library(lqmm)
library(ggplot2)
source("limit_exclusion_selection.R")
source("score_generation_functions.R")
source("decile_limit_integration.R")
source('condition_parental_pgs_integration.R')

# Sample disease parameters
pgs <- 0.131
h_ps <- sqrt(pgs)
h2 <- 0.3 
k <- 0.005 # quite low prevalence. numerical errors? 
max_N <- 5 # for debug ! 20
nsims <- 25000 # 0

# Set score generation parameters
cor <- matrix(1); M <- 1
matrices <- generate_matrices_for_sim(cor, pgs, h2, max_N)
Sigma <- matrices$Sigma # the PGS covariance matrix
Sigma_G <- matrices$Sigma_G # The cov matrix for the genetic variance not captured by the PGS
E_Matrix <- matrices$E_Matrix # The remaining variance covariance matrix
Kin_Matrix <- matrices$Kin_Matrix # The N x N relatedness matrix; 1's on the diag, 0.5 on the off-diag
E <- Sigma_G + E_Matrix # Error matrix is the sum of E_G and E_E
threshold <- qnorm(1-k)

# Run simulation
decile_limit_selection_stats <- as.data.frame(matrix(0, nrow = max_N, ncol = 3))
colnames(decile_limit_selection_stats) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")

for(i in 1:nsims) {
  if(i%%1000==0){
    print(c("i=", i))
  }
  score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, max_N)
  PGS_scores <- score_matrices$PGS_scores
  D_scores <- score_matrices$D_scores
  pass_to_data_frame <- decile_exclusion_selection(pgs, PGS_scores, D_scores, max_N)
  decile_limit_selection_stats[1:max_N, 1:3] <- decile_limit_selection_stats[1:max_N, 1:3] + pass_to_data_frame[1:max_N, 1:3]
}

for(p in seq(1, max_N, 1)) {
  decile_limit_selection_stats$N_embryos[p] <- p
  decile_limit_selection_stats$selected_average_disease[p] <- decile_limit_selection_stats$Selected_D_sum[p]/nsims
  decile_limit_selection_stats$random_average_disease[p] <- decile_limit_selection_stats$Random_D_sum[p]/nsims
  decile_limit_selection_stats$ARR_expected[p] <- k - decile_limit_selection_stats$selected_average_disease[p]
  decile_limit_selection_stats$RRR_expected[p] <- 1 - decile_limit_selection_stats$selected_average_disease[p]/k
  decile_limit_selection_stats$integration_RRR[p] <- 1 - (disease_risk(p, h_ps, 1, k, 0.1, 1)/k)
  decile_limit_selection_stats$integration_ARR[p] <- decile_limit_selection_stats$integration_RRR[p]*k
}

ggplot(data = decile_limit_selection_stats, aes(x = N_embryos)) + ggtitle(label = "Decile Exclusion Gain: Simulation vs. Integration") + 
  geom_point(aes(y = RRR_expected, color = "Simulation")) + 
  geom_point(aes(y = integration_RRR, color = "Integration")) +
  scale_color_manual(name = "Gain Calculation Method", values = c("Integration" = "blue", "Simulation" = "red"))
