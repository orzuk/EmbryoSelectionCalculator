
decile_exclusion_selection <- function(pgs, PGS_scores, D_scores, N_hi) {
  # the temporary dataframe to store values for embryos N_lo through N_hi
  pass_to_sim_frame <- as.data.frame(matrix(nrow = N_hi, ncol = 3))
  colnames(pass_to_sim_frame) <- c("N_embryos", "Selected_D_sum", "Random_D_sum")
  store_embryo_stats <- as.data.frame(matrix(nrow = N_hi, ncol = 3)); colnames(store_embryo_stats) <- c("Embryo_num", "Top_decile", "D_sum")
  decile_risk_score <- qnorm(0.9, 0, sqrt(pgs)) # 90th percentile of genetic risk score distribution is our cutoff
  decile_boolean <- 0 # indicates if we've found an embryo with a PGS in the top decile; if TRUE, means that we'll be running the first conditional from now on
  
  for(j in 1:N_hi) {
    # Add each new embryo for 1 through N_hi to the data frame
    pass_to_sim_frame$N_embryos[j] <- j # this is the number of embryos being selected from
    store_embryo_stats$Embryo_num[j] <- j # take the newest embryo's index
    store_embryo_stats$D_sum[j] <- D_scores[j] # store the D_sum of this embryo
    store_embryo_stats$Top_decile[j] <- ifelse(PGS_scores[j] >= decile_risk_score, 1, 0) # boolean for if more than one embryo is in the top decile
    if(decile_boolean == 0) {
      if(sum(store_embryo_stats$Top_decile[1:j]) >= 1 & sum(store_embryo_stats$Top_decile[1:j]) != j) {
        decile_boolean <- 1
      }
    }
    if(decile_boolean == 1 ) { # if at least one, but not all, of the embryos is in the top decile
      below_top_decile <- subset(store_embryo_stats, Top_decile == 0) # take a subset of the data where top_decile = F
      ran_num <- sample(1:nrow(below_top_decile), 1) # randomly sample from this subset
      pass_to_sim_frame$Selected_D_sum[j] <- below_top_decile$D_sum[ran_num] # and store this selection as selected_D_sum
      ran_num <- sample(1:j, 1) # then take a random sample from _all_ the embryos (not just the subset)
      pass_to_sim_frame$Random_D_sum[j] <- store_embryo_stats$D_sum[ran_num] # and store the selected value as Random_D_sum
      
    } else if(sum(store_embryo_stats$Top_decile[1:j]) == 0) { # if none of the embryos are in the top decile 
      ran_num_1 <- sample(1:j, 1) # pick a random embryo
      ran_num_2 <- sample(1:j, 1) # pick another random embryo
      pass_to_sim_frame$Random_D_sum[j]<- store_embryo_stats$D_sum[ran_num_1] # and store it for both random and selected
      pass_to_sim_frame$Selected_D_sum[j] <- store_embryo_stats$D_sum[ran_num_2]
    
    } else { # if all PGS are in the top decile
      ran_num_1 <- sample(1:j, 1) # pick a random embryo
      lowest_score_index <- which.min(PGS_scores[1:j]) # find the lowest PGS score
      pass_to_sim_frame$Random_D_sum[j]<- store_embryo_stats$D_sum[ran_num_1] # store the random D score
      pass_to_sim_frame$Selected_D_sum[j] <- store_embryo_stats$D_sum[lowest_score_index] # pass the D score of the embryo with the lowest PGS
    }
  }
  return(pass_to_sim_frame)
}