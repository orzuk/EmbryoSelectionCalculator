library("MASS")

## INPUT ##

# an N-max by M matrix of PGS_scores
# an N-max by M matrix of D scores
# K: a vector of disease prevalences
# E: E_Matrix + E_G, used to calculate the truncated MVN probability
# N_hi: the highest number of embryos to evaluate
multiple_selection_function <- function(PGS_scores, D_scores, E, threshold, N_hi) {
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
