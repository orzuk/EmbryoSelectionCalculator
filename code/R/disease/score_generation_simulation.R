#### README#######

# This file generates risk scores using 'multiple_disease_functions.R' to test 
# that the embryo score generation works correctly.


source("score_generation_functions.R")
library(MASS)
library(lqmm)
set.seed(1234)


### SETTING SAMPLE PARAMETERS ###

# Get sample disease genetic correlation matrix:
cor <- matrix(c(1.000, 0.3, 0.4, 0.5, 0.6,
                0.3, 1.000, 0.5, 0.3, 0.2,
                0.4, 0.5, 1.000, 0.400, 0.60,
                0.5, 0.3, 0.400, 1.000, 0.40,
                0.6, 0.2, 0.600, 0.400, 1.00), nrow = 5, ncol = 5)


M <- nrow(cor)
## Generate a vector of sample h^2_pgs
pgs <- runif(M, 0.1, 0.5)

## generate the total heritabilities for the diseases.
h2 <- pgs + runif(M, 0.1, 0.3)

# Choose a number of embryos
N_embryos <- 10

# generate disease prevalences
K <- runif(M, 0.01, 0.1)
threshold <- qnorm(1-K)



#### RUNNING THE SIMULATION ITSELF ####

### Generate the matrices for the simulation:
test <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)
Sigma <- test[[1]]; Sigma_G <- test[[2]]; E_Matrix <- test[[3]]; Kin_Matrix <- test[[4]];

### The simulation
iters <- 10000
disease_tracker <- matrix(nrow = iters, ncol = M, 0) # track the number of embryos with scores over the threshold, per disease, per round of embryo generation.
Z_list <- list(n = iters) ## store each round of Z scores
score_list <- list(n = iters) ## store each round of total scores

for(i in 1:iters) {
  disease_scores <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold)
  Z_list[[i]]<- disease_scores[[1]]
  score_list[[i]] <- disease_scores[[2]]
  D_scores <- disease_scores[[3]]
  disease_tracker[i, ] <- colSums(D_scores)
}


########### TESTS ########

# Ensure prevalences of diseases among embryos equal the expected prevalences
# Can increase the number of iterations to make the difference smaller
empirical_prevalences <- colSums(disease_tracker)/(N_embryos*iters)
empirical_prevalences; K
empirical_prevalences - K

# Make sure that the empirical correlations among the diseases equal the expected correlations
embryo_df <- as.data.frame(matrix(nrow = 0, ncol = M)) # a df to store all the embryo scores

for(i in 1:length(score_list)) {
  embryo_df <- rbind(embryo_df, score_list[[i]])
}

empirical_cor <- cor(embryo_df)
final_matrix <- E_Matrix + Sigma + Sigma_G
empirical_cor; final_matrix
cor(embryo_df) - final_matrix


















