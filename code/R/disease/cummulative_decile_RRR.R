source("score_generation_functions.R")
source("decile_exclusion_selection.R")
source("trunc_selection_function.R")
library("MASS")
library("ggplot2")
library("lqmm")


### IBD
h2_pgs_IBD <- 0.131
pgs <- h2_pgs_IBD
k_IBD <- 0.005
N_embryos <- 1
cor <- matrix(c(1))
M <- 1 # the number of diseases
q <- 0.1; zq_ps <- qnorm(q, 0, sqrt(pgs))
IBD_expected_gain <- 1-((integrate(trunc_selection_function, zq_ps, Inf, q, k_IBD, pgs)$value)/k_IBD) # expected relative gain


nsims3 <- 1000000
families_IBD <- as.data.frame(matrix(nrow = nsims3, ncol = 3))
colnames(families_IBD) <- c("Number", "Disease_Status", "Cummulative_RRR")

threshold <- qnorm(1-k_IBD)
decile_threshold <- qnorm(0.9, 0, sqrt(h2_pgs_IBD))
h2 <- pgs + 0.4
matrices <- generate_matrices_for_sim(cor, h2_pgs_IBD, h2, N_embryos)  
Sigma <- matrices$Sigma # the PGS covariance matrix
Sigma_G <- matrices$Sigma_G # The cov matrix for the genetic variance not captured by the PGS
E_Matrix <- matrices$E_Matrix # The remaining variance covariance matrix
Kin_Matrix <- matrices$Kin_Matrix # The N x N relatedness matrix; 1's on the diag, 0.5 on the off-diag
E <- Sigma_G + E_Matrix # Error matrix is the sum of E_G and E_E
cummulative_sick_families <- 0


for(i in 1:nrow(families_IBD)) {
  boolean <- 0
  selected_D_score <- 0
  while(boolean == 0) {
    score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, N_embryos)
    PGS_scores <- score_matrices$PGS_scores
    D_scores <- score_matrices$D_scores
    if(PGS_scores[1] < decile_threshold) {
      selected_D_score <- D_scores[1]
      boolean <- 1
    }
  }
  cummulative_sick_families <- cummulative_sick_families + selected_D_score # update X_s
  families_IBD$Number[i] <- i # update index of X
  families_IBD$Disease_Status[i] <- selected_D_score # The disease status of family X
  families_IBD$Cummulative_RRR[i] <- 1 - ((cummulative_sick_families/families_IBD$Number[i])/k_IBD) # 1 - ((X_s/X)/K)
}

subset <- families_IBD[seq(1, nsims3, 100),] # make plotting the data a bit easier; might want to start the sequence at 100,000 or so to zoom in on where the graph converges.

ggplot(data = subset, aes(x = Number, y = Cummulative_RRR)) + geom_point() + geom_hline(yintercept = IBD_expected_gain)


