setwd("Embryo Sims")
source("score_generation_functions.R") 
library(MASS)
library(ggplot2)
library(lqmm)
library(dbplyr)
library(mvtnorm)
library(beepr)
## INPUT ##

# an N-max by M matrix of PGS_scores
# an N-max by M matrix of D scores
# K: a vector of disease prevalences
# E: E_Matrix + E_G, used to calculate the truncated MVN probability
# N_hi: the highest number of embryos to evaluate
  
  
### SELECTION FUNCTION ####
#### YOU MUST START WITH 1 EMBRYO; FUNCTION WILL THROW ERROR IF YOU DON'T
speedy_selection <- function(PGS_scores, D_scores, E, threshold, N_hi) {
  # the temporary dataframe to store values for embryos N_lo through N_hi
  pass_to_sim_frame <- as.data.frame(matrix(nrow = N_hi, ncol = 4))
  colnames(pass_to_sim_frame) <- c("N_embryos", "Selected_D_sum", "Random_D_sum", "Predicted_R_0_selected")
  store_embryo_stats <- as.data.frame(matrix(nrow = N_hi, ncol = 3)); colnames(store_embryo_stats) <- c("Embryo_num", "Predicted_R_0", "D_sum")
  for(j in 1:N_hi) {
    pass_to_sim_frame$N_embryos[j] <- j # this is the number of embryos being selected from
    store_embryo_stats$Embryo_num[j] <- j # this is the index of the newest embryo
    # start with completely empty frame
    ran_num <- sample(1:j, 1)
    # If the randomly selected embryo hasn't yet been filled out, fill it out
    store_embryo_stats$Predicted_R_0[j] <- pmvnorm(lower = -Inf, upper = threshold, mean = PGS_scores[j,], sigma = E)[[1]]
    store_embryo_stats$D_sum[j] <- sum(D_scores[j,]) #store D_sum of newest embryo
    pass_to_sim_frame$Random_D_sum[j] <- store_embryo_stats$D_sum[ran_num] # store D_sum of randomly selected embryo
    if(j == 1) {
      pass_to_sim_frame$Predicted_R_0_selected[j] <- store_embryo_stats$Predicted_R_0[j]
      pass_to_sim_frame$Selected_D_sum[j] <- store_embryo_stats$D_sum[j]
    } else {
      if(store_embryo_stats$Predicted_R_0[j] > pass_to_sim_frame$Predicted_R_0_selected[j-1]) {
        pass_to_sim_frame$Predicted_R_0_selected[j] <- store_embryo_stats$Predicted_R_0[j]
        pass_to_sim_frame$Selected_D_sum[j] <- store_embryo_stats$D_sum[j]
      } else {
        pass_to_sim_frame$Predicted_R_0_selected[j] <- pass_to_sim_frame$Predicted_R_0_selected[j-1]
        pass_to_sim_frame$Selected_D_sum[j] <- pass_to_sim_frame$Selected_D_sum[j-1]
      }
    }
  }
  return(pass_to_sim_frame)
}



#### THE SIMULATION ####
sim_frame <- expand.grid(N_embryos = seq(1, 20, 1), pgs = seq(0.01, 0.61, 0.05), random_average_disease = 0, selected_average_disease = 0)

### Specify parameters
cor <- matrix(0.1, nrow = 5, ncol = 5)
diag(cor) <- 1
M <- nrow(cor)
K <- runif(M, 0.09, 0.1)
threshold <- qnorm(1-K)
nsims <- 2000
N_hi <- max(sim_frame$N_embryos)

## Running the simulation
system.time(
for(a in seq(N_hi, nrow(sim_frame-N_hi), by = N_hi)) {
  
  selection_stats <- as.data.frame(matrix(0, nrow = N_hi, ncol = 3))
  
  # the total D_sums that are eventually used to calculate disease rates
  colnames(selection_stats) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")
  pgs <- rep(sim_frame$pgs[a], M) #can set this to be a vector of different values
  h2 <- pgs + 0.1
  N_embryos <- sim_frame$N_embryos[a]
  matrices <- generate_matrices_for_sim(cor, pgs, h2, N_embryos)  
  Sigma <- matrices$Sigma
  Sigma_G <- matrices$Sigma_G
  E_Matrix <- matrices$E_Matrix
  Kin_Matrix <- matrices$Kin_Matrix
  E <- Sigma_G + E_Matrix 
  system.time(
    for(i in 1:nsims) {
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
    sim_frame$selected_average_disease[p] <- selection_stats$Selected_D_sum[index]/nsims
    sim_frame$random_average_disease[p] <- selection_stats$Random_D_sum[index]/nsims
    index <- index + 1
  }
}
)
beep("coin")
sim_frame$average_absolute_disease_reduction <- sim_frame$random_average_disease - sim_frame$selected_average_disease
sim_frame$average_relative_disease_reduction <- sim_frame$average_absolute_disease_reduction/sim_frame$random_average_disease


## Runtime stats
## N = [1:3], PGS = seq(0.01, 0.71, 0.05), nsims = 10000 --> 22 minutes
## N = [1:15], PGS = seq(0.01, 0.6, 0.05), nsims = 2000 --> 27 minutes
## N = [1:20], "", "" --> 34 minutes


### PLOTS ####

## Confirm that the disease prevalence is equal for N=1 thru N_hi and stays constant as a function of h^2_{pgs}
ggplot(data = sim_frame, aes(x = pgs, y = random_average_disease)) + geom_point() + stat_smooth(method = "lm") + facet_wrap(~N_embryos)


### Visualization of Risk Reduction
# Absolute decrease in average disease
ggplot(data = sim_frame, aes(x = pgs, y = average_absolute_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) +
  stat_smooth(method = "lm") + xlab("h^2_{pgs}") + ylab("Absolute Reduction in Average Disease Incidence Per Embryo") +
  ggtitle("Absolute Reduction in Disease Incidence per Embryo as Function of PGS and N")


### polynomial fits on a single plot, 
ggplot(data = sim_frame, aes(x = pgs, y = average_absolute_disease_reduction, colour = as.factor(N_embryos))) + geom_point() + stat_smooth(method= "lm", formula = y ~ poly(x, 2)) +
   xlab("h^2_{pgs}") + ylab("Absolute Reduction in Average Disease Incidence Per Embryo") +
  ggtitle("Absolute Reduction in Disease Incidence per Embryo as Function of PGS and N")


# Relative decrease in disease
ggplot(data = sim_frame, aes(x = pgs, y = average_relative_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) + 
  stat_smooth(method="lm") + xlab("h^2_{pgs}") + ylab("Relative Reduction in Average Disease Incidence Per Embryo") +
  ggtitle("Relative Reduction in Disease Incidence per Embryo as Function of PGS and N")

# look at just N_embryos = 15
subset_data <- subset.data.frame(sim_frame, N_embryos == 15)
ggplot(data = subset_data, aes(x = pgs, y = average_relative_disease_reduction)) + geom_point() + facet_wrap(~N_embryos) + geom_smooth(method = "lm")


# function to predict how much absolute increase in R_0 you'd get from selecting using PGS vs. at random
risk_reduction_model <- lm(data = sim_frame, average_relative_disease_reduction ~ pgs + N_embryos + N_embryos:pgs)
summary(risk_reduction_model) # significant N_embryo x pgs interaction on absolute increase in R_0

# same function, but relative increase
risk_reduction_model <- lm(data = sim_frame, proportional_R_0_increase ~ pgs + N_embryos + N_embryos:pgs)
summary(risk_reduction_model) # significant N_embryo x pgs interaction for proportional increase in R_0





