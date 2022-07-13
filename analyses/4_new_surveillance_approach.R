# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(zoo)

# Loading in custom functions
source(here("functions", "time_series_characterisation_functions.R"))

# Loading in Metadata and Location Data
ts_metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE)))
cluster_membership <- readRDS(here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_membership[, c("id", "cluster")]
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))
overall <- overall %>%
  left_join(cluster_membership, by = "id")

# Extracting Mean Realisations
id <- overall$id
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(id), ncol = (12 * interpolating_points + 1))
overdisp <- vector(mode = "numeric",length = length(id))
for (i in 1:length(id)) {
  
  index <- id[i]
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
  }  
  
  # Extracting the mean fitted time-series
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  ordered_timepoints <- all_timepoints[order(all_timepoints)]
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  negbinom_intensity_mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  mean_realisation[i, ] <- negbinom_intensity_mean
  overdisp[i] <- mean(MCMC_output[, "overdispersion"])
}

# Interpolate to Get Daily Biting Rate and then Normalise By Total Density
approx_min_output <- function(x) {
  temp <- approx(x, n = 365)
  return(temp$y)
}
conv_daily_to_monthly <- function(x) {
  month_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  temp <- c()
  for (i in 1:12) {
    if (i == 1) {
      temp <- c(temp, mean(x[1:month_length[i]]))
    } else {
      temp <- c(temp, mean(x[sum(month_length[1:(i-1)]):sum(month_length[1:i])]))
    }
  }
  return(temp)
}
daily_vector_density <- t(apply(mean_realisation, 1, approx_min_output))
normalised_daily_vector_density <- t(apply(daily_vector_density, 1, normalise_total))
monthly_vector_density <- t(apply(normalised_daily_vector_density, 1, conv_daily_to_monthly))
normalised_monthly_vector_density <- t(apply(monthly_vector_density, 1, normalise_total))

# Surv Stuff - see https://en.wikipedia.org/wiki/Poisson_binomial_distribution

# want to vary:
#   -> effort (number of nights sampled and how many months you sample for)
#   -> timing (when you sample)
#   -> metric for ID (number of mosquitoes expected to be caught)

# surv strategies
#   -> random contiguous (i.e. sequential months)
#   -> random non-contiguous (i.e. random months over the course of the year)
#   -> ento targeted contiguous (sequential months, knowledge of vector peak)
#   -> ento targeted non-contiguous (random months, knowledge of vector peak)
#   -> rainfall targeted contiguous (sequential months, 1 month before/after rainfall peak)

EIR <- 5
sporozoite_prev <- 0.05
ABR <- EIR/sporozoite_prev
num_months_sampled <- c(1, 2, 3, 4, 5, 6) # varying effort 
num_nights_per_month <- c(1, 2, 3, 4, 5, 6) # varying effort

# For single time series and contiguous sampling (so 12 possible combinations of sampled months all differing
# by starting point)
prob_not_sampled <- array(data = NA, dim = c(6, 6, 12))
monthly_density <- normalised_monthly_vector_density[1, ]
monthly_prob_sampled <- test_month_dens * ABR/c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) # not sure if this is quite right - need to double check
for (i in 1:length(num_months_sampled)) {

  # Extract num months to be sampled and Generate matrix of all possible contiguous months to sample at 
  months_sampled <- num_months_sampled[i]
  sampling_mat <- matrix(0, 12, months_sampled)
  sampling_mat <- col(sampling_mat) + row(sampling_mat) - 1
  sampling_mat[sampling_mat > 12] <- sampling_mat[sampling_mat > 12] - 12
  
  # Extract num nights sampling per month for each run and duplicate entries in matrix as required
  for (j in 1:length(num_nights_per_month)) {
    
    nights_per_month <- num_nights_per_month[j]
    if (months_sampled == 1 & nights_per_month == 1) {
      overall_sampling_mat <- t(t(apply(sampling_mat, 1, rep, each = nights_per_month)))
    } else {
      overall_sampling_mat <- t(apply(sampling_mat, 1, rep, each = nights_per_month))
    }
    
    for (k in 1:12) {
      temp_probs <-  monthly_prob_sampled[overall_sampling_mat[k, ]]
      prob_not_sampled[i, j, k] <- prod(1 - temp_probs)
    }
  }
}


# Reordering mean fitted time-series to start at max
reordered_mean_realisation <- matrix(nrow = length(id), ncol = (12 * interpolating_points + 1))
start_index <- apply(normalised_output, 1, function(x) which(x == max(x)))
end_index <- dim(reordered_mean_realisation)[2]
for (i in 1:length(id)) {
  reordered_mean_realisation[i, ] <- normalised_output[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
}
one_output <- reordered_mean_realisation[cluster_membership$cluster == 1, ]
mean_one_output <- apply(one_output, 2, mean)

two_output <- reordered_mean_realisation[cluster_membership$cluster == 2, ]
mean_two_output <- apply(two_output, 2, mean)

aug_norm_output <- rbind(normalised_output, mean_one_output, mean_two_output)

summary_probs <- array(dim = c(length(id) + 2, 12, 3))
p_overall <- 1
for (k in 1:67) {
  
  # Generating Monthly Relative Probabilities of Finding Stephensi and Multiplying By Defined Overall Prob Detection Given Present
  time_series <- aug_norm_output[k, ]
  avg_month_prob <- vector(mode = "numeric", length = 12L)
  counter <- 1
  for (i in 1:length(avg_month_prob)) {
    if (i < 12) {
      avg_month_prob[i] <- mean(time_series[c(counter, counter + 1)])
      counter <- counter + 2
    } else {
      avg_month_prob[i] <- mean(time_series[c(counter, counter + 1, counter + 2)])
      counter <- counter + 2
    }
  }
  p_zero_caught <- 1 - (avg_month_prob/max(avg_month_prob) * p_overall) # option1 <- Alt = 1 - (test * 0.5)
  
  # Iterating Over All Possible Start Months and All Possible Number of Consecutive Months Sampled
  #   rows are the month we start at
  #   columns are the number of months we consider
  prob_matrix <- matrix(data = NA, nrow = 12, ncol = 12)
  for (i in 1:12) {
    for (j in 1:12) {
      if ((i+j-1) > 12) {
        indices <- c(i:12, 1:(j - length(i:12)))
      } else {
        indices <- i:(i+j-1)
      }
      temp_p <- prod(p_zero_caught[indices])
      prob_matrix[i, j] <- temp_p
    }
  }
  temp_median_prob <- apply(prob_matrix, 2, quantile, prob = c(0.25, 0.5, 0.75))
  summary_probs[k, , ] <- t(temp_median_prob)
  print(k)
}

lower_probs <- summary_probs[1:65, , 1]
median_probs <- summary_probs[1:65, , 2]
upper_probs <- summary_probs[1:65, , 3]

cluster_one_median <- data.frame(time = seq(1, 12), median = summary_probs[66, , 2], lower = summary_probs[66, , 1], upper = summary_probs[66, , 3])
cluster_two_median <- data.frame(time = seq(1, 12), median = summary_probs[67, , 2], lower = summary_probs[67, , 1], upper = summary_probs[67, , 3])

median_summary <- data.frame(id = rep(id, 3), stat = c(rep("med", 65), rep("low", 65), rep("high", 65)),
                             cluster = rep(cluster_membership$cluster, each = 3), 
                             rbind(median_probs, lower_probs, upper_probs)) %>%
  pivot_longer(cols = X1:X12, names_to = "time", values_to = "prob") %>%
  mutate(time = as.numeric(gsub("X", "", time))) %>%
  group_by(time, stat, cluster) %>%
  summarise(median = median(prob),
            lower = min(prob),
            upper = max(prob)) %>%
  pivot_wider(names_from = "stat", 
              values_from = median:upper)#
medians <- data.frame(id = id, median_probs) %>%
  pivot_longer(cols = X1:X12, names_to = "time", values_to = "prob") %>%
  mutate(time = as.numeric(gsub("X", "", time)))

fig4c <- ggplot(median_summary) +
  geom_ribbon(aes(x = time, ymin = lower_med, ymax = upper_med, fill = factor(cluster)), alpha = 0.1) +
  geom_path(aes(x = time, y = median_med, colour = factor(cluster)), size = 1) +
  scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12)) +
  scale_fill_manual(values = palette()[2:1]) + 
  scale_colour_manual(values = palette()[2:1]) + 
  labs(y = "Probability of Missing Anopheles Stephensi", x = "Number of Consecutive Months Sampled") +
  theme_bw() +
  theme(legend.position = "none")
