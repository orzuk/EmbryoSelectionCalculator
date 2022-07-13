source("score_generation_functions.R")
source("decile_exclusion_selection.R")
library("MASS")
library("ggplot2")
library("lqmm")

nsims2 <- 1000000
IBD_data <- expand.grid(N_embryos = seq(1, 2, 1), pgs = 0.131, K = 0.005, random_average_disease = 0, selected_average_disease = 0, fraction_0 = 0, fraction_1 = 0, fraction_2 = 0)
cor <- matrix(c(1))
N_hi <- max(IBD_data$N_embryos) # redundant, yes
M <- 1
fraction_tracker <- as.data.frame(matrix(0, nrow = 3, ncol = 5))
rownames(fraction_tracker) <- c("0_Embryos", "1_Embryos", "2_Embryos")
colnames(fraction_tracker) <- c("Embryos_in_top_decile", "Selected_diseased", "Random_diseased", "ARR", "RRR")
K <- 0.005

## Running the simulation

for(a in seq(N_hi, nrow(IBD_data), by = N_hi)) {
  
  selection_stats <- as.data.frame(matrix(0, nrow = N_hi, ncol = 3))
  # the total D_sums that are eventually used to calculate disease rates
  colnames(selection_stats) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")
  pgs <- IBD_data$pgs[a] #can set this to be a vector of different values
  h2 <- pgs + 0.4 # value doesn't matter here as long as h2 < 1
  N_embryos <- IBD_data$N_embryos[a]
  matrices <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)  
  Sigma <- matrices$Sigma # the PGS covariance matrix
  Sigma_G <- matrices$Sigma_G # The cov matrix for the genetic variance not captured by the PGS
  E_Matrix <- matrices$E_Matrix # The remaining variance covariance matrix
  Kin_Matrix <- matrices$Kin_Matrix # The N x N relatedness matrix; 1's on the diag, 0.5 on the off-diag
  E <- Sigma_G + E_Matrix # Error matrix is the sum of E_G and E_E
  threshold <- qnorm(1-IBD_data$K[a])
  decile_risk_score <- qnorm(0.9, 0, sqrt(pgs))
  
  for(i in 1:nsims2) {
    score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, N_hi)
    PGS_scores <- score_matrices$PGS_scores
    Z_scores <- score_matrices$Z_scores
    D_scores <- score_matrices$D_scores
    number_embryo_in_decile <- sum(PGS_scores > decile_risk_score)
    fraction_tracker$Embryos_in_top_decile[number_embryo_in_decile + 1] <-  fraction_tracker$Embryos_in_top_decile[number_embryo_in_decile + 1] + 1
    pass_to_data_frame <- decile_exclusion_selection(pgs, PGS_scores, D_scores, N_hi)
    fraction_tracker$Selected_diseased[number_embryo_in_decile + 1] <- fraction_tracker$Selected_diseased[number_embryo_in_decile + 1] + pass_to_data_frame$Selected_D_sum[N_hi] # selected D sum for N = 2
    fraction_tracker$Random_diseased[number_embryo_in_decile + 1] <- fraction_tracker$Random_diseased[number_embryo_in_decile + 1] + pass_to_data_frame$Random_D_sum[N_hi] # random D sum for N = 2
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

fraction_tracker$prop_selected <- fraction_tracker$Selected_diseased/fraction_tracker$Embryos_in_top_decile
fraction_tracker$prop_random <- fraction_tracker$Random_diseased/fraction_tracker$Embryos_in_top_decile
fraction_tracker$percentage_total <- fraction_tracker$Embryos_in_top_decile/nsims2


fraction_tracker$ARR <- (fraction_tracker$Random_diseased - fraction_tracker$Selected_diseased)/fraction_tracker$Embryos_in_top_decile # absolute reduction 
average_ARR <- as.matrix(t(fraction_tracker$ARR)) %*% as.matrix(fraction_tracker$percentage_total) # average ARR
average_RRR <- average_ARR/k # ARR/k

