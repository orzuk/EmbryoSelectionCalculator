library(ggplot2)
# Set sample parameters
k <- 0.005; M <- 1
pgs <- 0.131; h2 <- 0.3; threshold <- qnorm(1-k)
max_embryos <- 5

# generate decile mean parental PGS values
deciles <- seq(0.05, 0.95, 0.1)
memoize_parental_pgs_decile <- qnorm(deciles, 0, sqrt(pgs))

# for each decile score, generate a shit ton of sims and blah

nsims <- 500000
decile_bin_selection_results <- as.data.frame(matrix(0, nrow = 10, ncol = 8))
colnames(decile_bin_selection_results) <- c("decile", "mid_parental_pgs", "random_disease", "selected_disease", "base_prevalence", "selected_prevalence", "RRR", "ARR")
decile_bin_selection_results$mid_parental_pgs <- memoize_parental_pgs_decile
decile_bin_selection_results$decile <- seq(1, 10, 1)
embryo_list <- list(length(nsims*10)) # store embryos generated

# Simulate
for(i in 1:nrow(decile_bin_selection_results)) {
  mean_parental_pgs <- decile_bin_selection_results$mid_parental_pgs[i]
  for(j in 1:nsims) {
    x_scores <- rnorm(max_embryos, mean_parental_pgs, sqrt(pgs/2)) 
    e_scores <- rnorm(max_embryos, 0, sqrt(1-pgs))
    z_scores <- x_scores + e_scores
    d_scores <- z_scores > threshold
    decile_bin_selection_results$selected_disease[i] <- decile_bin_selection_results$selected_disease[i] + d_scores[which.min(x_scores)] # looking only at N = 5
    decile_bin_selection_results$random_disease[i] <- decile_bin_selection_results$random_disease[i] + d_scores[5]
  }
  decile_bin_selection_results$base_prevalence[i] <- decile_bin_selection_results$random_disease[i]/nsims
  decile_bin_selection_results$selected_prevalence[i] <- decile_bin_selection_results$selected_disease[i]/nsims
  decile_bin_selection_results$ARR[i] <- decile_bin_selection_results$base_prevalence[i] - decile_bin_selection_results$selected_prevalence[i]
  decile_bin_selection_results$RRR[i] <- 1 - decile_bin_selection_results$selected_prevalence[i]/decile_bin_selection_results$base_prevalence[i]
}


ggplot(data = decile_bin_selection_results) + geom_point(aes(x = as.factor(decile), y = RRR), color = "red") +
  ggtitle(paste("Relative Risk Reduction as Function of Parental Mean PGS \n h^2_pgs:", pgs, " k:", k, "N:", 5, "Simulations:", nsims)) +
  xlab("Mean Parental PGS Decile") + ylab("Relative Risk Reduction")

ggplot(data = decile_bin_selection_results) + 
  geom_point(aes(x = as.factor(decile), y = base_prevalence, color = "Unselected")) +
  geom_point(aes(x = as.factor(decile), y = selected_prevalence, color = "Selected")) + 
  ggtitle(paste("Absolute Risk Reduction as Function of Parental Mean PGS \n h^2_pgs:", pgs, "  k:", k, "  N:", 5, "  Simulations:", nsims)) +
  xlab("Mean Parental PGS Decile") + ylab("Absolute Risk Reduction") + scale_color_manual(name = "Prevalence", values = c("Unselected" = "red", "Selected" = "blue")) + 
  geom_segment(aes(x = as.factor(decile), xend = as.factor(decile), y = selected_prevalence, yend = base_prevalence))

