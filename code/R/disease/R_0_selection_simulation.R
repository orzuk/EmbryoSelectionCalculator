setwd('C:\\Users\\Or Zuk\\Dropbox\\EmbryoSelection\\Code\\R\\Adam') # change to your path with the code
source("simulate_risk_prediction.R")
source("score_generation_functions.R")
library(MASS)
library(ggplot2)
library(dbplyr)

### SELECTION FUNCTION 1: MAXIMIZE PREDICTED R_0 #### 

# embryos = list of 2 matrices, the Z scores and the total liability scores
# K = prevalence
# E = covariance matrix of error

select_embryos_R_0 <- function(Z_scores, D_scores, K, E) {
  
  best_embryo_index <- 0 # index of selected embryo
  predicted_max_R_0 <- 0 # predicted R_0 of selected embryo
  selected_D_sum <- 0 # sum of disease score for selected embryo
  
  ran_num <-floor(runif(1, 1, nrow(Z_scores)+1))
  random_D_sum <- 0
  random_R_0 <- 0
  for(i in 1:nrow(Z_scores)) {
    Z_predicted_R_0 <- simulate_PGS_disease_risk(Z_scores[i,], K, E, n_sims = 50)[[1]] # The PGS predicted risk, which we use for selection; we want the R_0 element of this
    actual_D_sum <- sum(D_scores[i,]) # Sum of D scores, the actual predicted number of diseases
    if(Z_predicted_R_0 > predicted_max_R_0) { # If Z_predicted R_0 is better than current R_0
      best_embryo_index <- i 
      predicted_max_R_0 <- Z_predicted_R_0
      selected_D_sum <- actual_D_sum # save the ACTUAL predicted R_0
    }
    if(i == ran_num) {
      random_D_sum <- selected_D_sum
      random_R_0 <- Z_predicted_R_0
    }
  }
  embryo_list <- list(best_embryo_index, predicted_max_R_0, selected_D_sum, ran_num, random_R_0, random_D_sum)
  return(embryo_list)
}


#### EMBRYO SELECTION SIMULATION
#### Sample Parameters ####
cor <- matrix(c(1.000, 0.3, 0.4, 0.5, 0.6,
                0.3, 1.000, 0.5, 0.3, 0.2,
                0.4, 0.5, 1.000, 0.400, 0.60,
                0.5, 0.3, 0.400, 1.000, 0.40,
                0.6, 0.2, 0.600, 0.400, 1.00), nrow = 5, ncol = 5)
M <- nrow(cor) 
K <- runif(M, 0.01, 0.1)

### CRITERION 1: Maximize R_0 ######

## SIMULATION FRAME TO STORE VALUES
sim_frame <- expand.grid(pgs = seq(0.00, 0.6, 0.03), N_embryos = seq(3, 10, 1))

## Seeing how PGS and N_embryos affects the increase in R_0 with selected vs. random embryos
nsims <- 20 # 2000

for(j in 1:nrow(sim_frame)) {
  N_embryos <- sim_frame$N_embryos[j]
  pgs <- rep(sim_frame$pgs[j], 5)
  h2 <- pgs + 0.2 # for the sake of argument, but this could be specified
  matrices <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)  
  Sigma <- matrices[[1]]
  Sigma_G <- matrices[[2]]
  E_Matrix <- matrices[[3]]
  Kin_Matrix <- matrices[[4]]
  E <- matrices[[2]] + matrices[[3]] 
  D_sum_selected <- 0
  D_sum_random <- 0
  threshold <- qnorm(1-K)
  for (i in 1:nsims) {
    score_matrices <- generate_total_scores(Sigma, Sigma_G, E_Matrix, Kin_Matrix, threshold, N_embryos)
    selection_stats <- select_embryos_R_0(score_matrices[[1]], score_matrices[[3]], K, E)
    D_sum_selected<- D_sum_selected + selection_stats[[3]]
    D_sum_random <- D_sum_random + selection_stats[[6]]
  }
  sim_frame$selected_average_disease[j] <- D_sum_selected/nsims # The average # of diseases per embryo
  sim_frame$random_average_disease[j] <- D_sum_random/nsims
}
 
sim_frame$average_absolute_disease_reduction <- sim_frame$random_average_disease - sim_frame$selected_average_disease
sim_frame$average_relative_disease_reduction <- sim_frame$average_absolute_disease_reduction/sim_frame$random_average_disease



### Visualization of Risk Reduction
# Absolute decrease in average disease
ggplot(data = sim_frame, aes(x = pgs, y = average_absolute_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) +
  stat_smooth(method="lm", formula=y ~ poly(x, 2, raw=TRUE)) + xlab("PGS Score") + ylab("Absolute Reduction in Average Disease Incidence Per Embryo") +
  ggtitle("Absolute Reduction in Disease Incidence per Embryo as Function of PGS and N")

# Relative decrease
ggplot(data = sim_frame, aes(x = pgs, y = average_relative_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) + 
  stat_smooth(method="lm") + xlab("PGS Score") + ylab("Relative Reduction in Average Disease Incidence Per Embryo") +
  ggtitle("Relative Reduction in Disease Incidence per Embryo as Function of PGS and N")

# look at just N_embryos = 5
subset_data <- subset.data.frame(sim_frame, N_embryos == 5)
ggplot(data = subset_data, aes(x = pgs, y = average_relative_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) + geom_smooth(method = "lm")


# function to predict how much absolute increase in R_0 you'd get from selecting using PGS vs. at random
risk_reduction_model <- lm(data = sim_frame, average_relative_disease_reduction ~ pgs + N_embryos + pgs:N_embryos)
summary(risk_reduction_model) # significant N_embryo x pgs interaction on absolute increase in R_0

# same function, but relative increase
risk_reduction_model <- lm(data = sim_frame, proportional_R_0_increase ~ pgs + N_embryos + pgs:N_embryos)
summary(risk_reduction_model) # significant N_embryo x pgs interaction for proportional increase in R_0
