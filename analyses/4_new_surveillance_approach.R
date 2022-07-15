# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(zoo); library(RColorBrewer)

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
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2")) %>%
  left_join(cluster_membership, by = "id")

# Extracting the Mean Fitted Profile for Each Time Series
id <- overall$id
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(id), ncol = (12 * interpolating_points + 1))
overdisp <- vector(mode = "numeric",length = length(id))
for (i in 1:length(id)) {
  
  # Loading in and processing the fitted time-series
  index <- id[i]
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
  }  
  
  # Extracting the mean fitted temporal profile
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  negbinom_intensity_mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  mean_realisation[i, ] <- negbinom_intensity_mean
}

# Interpolate to Get Daily Biting Rate and then Normalise By Total Density
daily_vector_density <- t(apply(mean_realisation, 1, approx_min_output)) # interpolate to get daily biting rate
monthly_vector_density <- t(apply(daily_vector_density, 1, conv_daily_to_monthly)) # average over month
normalised_monthly_vector_density <- t(apply(monthly_vector_density, 1, normalise_total)) # normalise within each time-series
peak_density_month <- apply(normalised_monthly_vector_density, 1, function(x) which(x == max(x))) # get month of peak density

# Extracting Rainfall Data (also sort out leap year stuff)
months_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
rainfall_storage <- matrix(nrow = dim(overall)[1], ncol = length(months_length)) 
for (i in 1:length(overall$id)) {
  index <- overall$id[i]
  temp <- c()
  rf <- read.csv(here(paste0("data/location_specific_rainfall/rainfall_ts", index, ".csv")))
  rf <- rf %>%
    group_by(daymonth_id) %>%
    summarise(rainfall = mean(rainfall))
  counter <- 1
  count_vec <- counter
  for (j in 1:length(months_length)) {
    indices <- counter:(counter + months_length[j] - 1)
    temp <- c(temp, sum(rf$rainfall[indices]))
    counter <- counter + months_length[j]
  }
  rainfall_storage[i, ] <- temp
}
peak_rainfall_month <- apply(rainfall_storage, 1, function(x) which(x == max(x)))

# Exploring the Impact of Different Surveillance Strategies:
#   Premise behind this is that we want to vary:
#     -> effort (number of nights sampled and how many months you sample for)
#     -> timing (when you sample, whether months are contiguous etc)
#     -> metric for ID (number of mosquitoes expected to be caught)

# Factors to vary:
#  Number of months sampled
#  Number of nights within each month sampled
#  Whether sampling is done sequentially monthly or months are picked randomly
#  For sequentially, whether timing is random, rainfall driven or entomological peak driven

# surv strategies
#   -> random contiguous (i.e. sequential months)
#   -> random non-contiguous (i.e. random months over the course of the year)
#   -> ento targeted contiguous (sequential months, knowledge of vector peak)
#   -> ento targeted non-contiguous (random months, knowledge of vector peak)
#   -> rainfall targeted contiguous (sequential months, 1 month before/after rainfall peak)


# Simulating Strategies With Sequential Monthly Sampling 
EIR <- 5
sporozoite_prev <- 0.05
ABR <- EIR/sporozoite_prev
num_months_sampled <- 1:12 # varying effort 
num_nights_per_month <- 1:10 # varying effort

# For single time series and contiguous sampling (so 12 possible combinations of sampled months all differing
# by starting point)
prob_not_sampled_basic <- array(data = NA, dim = c(65, 12, 10, 12))
prob_not_sampled_poisson <- array(data = NA, dim = c(65, 12, 10, 12))
for (t in 1:65) {
  
  monthly_density <- normalised_monthly_vector_density[t, ]
  monthly_prob_sampled <- monthly_density * ABR/c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) # not sure if this is quite right - need to double check
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
        temp_probs[temp_probs > 1] <- 1
        prob_not_sampled_basic[t, i, j, k] <- prod(1 - temp_probs)
        
        temp_lambda_indiv <- monthly_prob_sampled[overall_sampling_mat[k, ]]
        temp_lambda_combined <- sum(temp_lambda_indiv) 
        prob_not_sampled_poisson[t, i, j, k] <- exp(-temp_lambda_combined) # p(k=0) from Poisson dist; equivalent to dpois(x = 0, lambda = temp_lambda_combined)
      }
    }
  }
  print(t)
}

# 1st dim is the time series
# 2nd dim is the number of months sampled
# 3rd dim is the number of nights per month sampled
# 4th dim is which month you start at 
not_sampled_vector_peak_absolute_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_rainfall_peak_absolute_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_rainfall_peak_plus_absolute_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_annual_avg_absolute_summary <- array(data = NA, dim = c(65, 6, 6))

not_sampled_vector_annual_diff_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_vector_rainfall_diff_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_vector_rainfall_peak_plus_diff_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_rainfall_annual_diff_summary <- array(data = NA, dim = c(65, 6, 6))
not_sampled_rainfall_peak_plus_annual_diff_summary <- array(data = NA, dim = c(65, 6, 6))

for (i in 1:65) {
  ts_peak_vector_month <- peak_density_month[i] 
  ts_peak_rainfall_month <- peak_rainfall_month[i] 
  ts_peak_rainfall_plus_month <- (ts_peak_rainfall_month + 1) %% 12
  ts_peak_rainfall_plus_month <- ifelse(ts_peak_rainfall_plus_month == 0, 12, ts_peak_rainfall_plus_month)
  
  not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6] <- 1- prob_not_sampled_poisson[i, 1:6, 1:6, ts_peak_vector_month]
  not_sampled_rainfall_peak_absolute_summary[i, 1:6, 1:6] <- 1 - prob_not_sampled_poisson[i, 1:6, 1:6, ts_peak_rainfall_month]
  not_sampled_rainfall_peak_plus_absolute_summary[i, 1:6, 1:6] <- 1 - prob_not_sampled_poisson[i, 1:6, 1:6, ts_peak_rainfall_plus_month]
  not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6] <- 1 - apply(prob_not_sampled_poisson[i, 1:6, 1:6, 1:12], c(1, 2), mean)
  
  not_sampled_vector_annual_diff_summary[i, 1:6, 1:6] <- not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6] - not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6]
  not_sampled_vector_rainfall_diff_summary[i, 1:6, 1:6] <- not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6] - not_sampled_rainfall_peak_absolute_summary[i, 1:6, 1:6]
  not_sampled_vector_rainfall_peak_plus_diff_summary[i, 1:6, 1:6] <- not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6] - not_sampled_rainfall_peak_plus_absolute_summary[i, 1:6, 1:6]
  not_sampled_rainfall_annual_diff_summary[i, 1:6, 1:6] <- not_sampled_rainfall_peak_absolute_summary[i, 1:6, 1:6] - not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6]
  not_sampled_rainfall_peak_plus_annual_diff_summary[i, 1:6, 1:6] <- not_sampled_rainfall_peak_plus_absolute_summary[i, 1:6, 1:6] - not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6]
}

apply(not_sampled_annual_avg_peak_diff_summary, c(2, 3), mean)
apply(not_sampled_rainfall_peak_diff_summary, c(2, 3), mean)
apply(not_sampled_vector_rainfall_peak_plus_diff_summary, c(2, 3), mean)

apply(not_sampled_vector_peak_absolute_summary, c(2, 3), mean)/
  apply(not_sampled_annual_avg_absolute_summary, c(2, 3), mean)

apply(not_sampled_vector_peak_absolute_summary, c(2, 3), mean)/
  apply(not_sampled_rainfall_peak_absolute_summary, c(2, 3), mean)



apply(not_sampled_annual_avg_peak_diff_summary, c(2, 3), mean)
apply(not_sampled_annual_avg_peak_diff_summary, c(2, 3), mean)/apply(not_sampled_annual_avg_absolute_summary, c(2, 3), mean)


which(not_sampled_rainfall_peak_ratio_summary[, 1, 1] == max(not_sampled_rainfall_peak_ratio_summary[, 1, 1]))

hist(not_sampled_rainfall_peak_ratio_summary[, 1, 1], xlim = c(0, ))


not_sampled_vector_peak_absolute_summary[1, , ]
not_sampled_annual_avg_absolute_summary[1, , ]

# peak vector
temp <- data.frame(1 - prob_not_sampled_poisson[6, 1:6, 1:6, ts_peak_vector_month])
temp$Y <- paste0("Y", 1:6)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:6))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
    position = "right",
    labels = 1:10) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")

# annual average
temp2 <- data.frame(1 - apply(prob_not_sampled_poisson[6, 1:6, 1:6, 1:12], c(1, 2), mean)) # average for the year
temp2$Y <- paste0("Y", 1:6)
temp2$Y <- factor(temp2$Y, levels = paste0("Y", 1:6))
temp_long2 <- temp2 %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long2$X  <- factor(temp_long2$X, levels = colnames(temp2)[-length(colnames(temp2))])
ggplot(temp_long2) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long2$X)),
    position = "right", 
    labels = 6:1) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")

# rainfall peak
temp3 <- data.frame(1 - prob_not_sampled_poisson[6, 1:6, 1:6, ts_peak_rainfall_month])
temp3$Y <- paste0("Y", 1:6)
temp3$Y <- factor(temp3$Y, levels = paste0("Y", 1:6))
temp_long3 <- temp3 %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long3$X  <- factor(temp_long3$X, levels = colnames(temp3)[-length(colnames(temp3))])
ggplot(temp_long3) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long3$X)),
    position = "right",
    labels = 6:1) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")

temp[, -7]/temp2[, -7]
temp[, -7]/temp3[, -7]

temp2[, -7]/temp3[, -7]

# Reordering rainfall time-series to start at max of vector density (note NOT max rainfall, which is what I calculate above)
# reordered_rainfall <- matrix(nrow = length(retain_index), ncol = length(months_length))
# end_rain_index <- dim(reordered_rainfall)[2]
# for (i in 1:length(retain_index)) {
#   reordered_rainfall[i, ] <- norm_rainfall_storage[i, c(start_index[i]:end_rain_index, 1:(start_index[i]-1))]
# }

# Reordering
# reordered_mean_realisation <- matrix(nrow = 65, ncol = (12 * interpolating_points + 1))
# start_index <- apply(mean_realisation, 1, function(x) which(x == max(x)))
# start_index_month <- round(start_index/25 * 12)
# end_index <- dim(reordered_mean_realisation)[2]
# for (i in 1:65) {
#   reordered_mean_realisation[i, ] <- mean_realisation[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
# }

# not_sampled_vector_annual_ratio_summary <- array(data = NA, dim = c(65, 6, 6))
# not_sampled_vector_rainfall_ratio_summary <- array(data = NA, dim = c(65, 6, 6))
# not_sampled_rainfall_peak_plus_ratio_summary <- array(data = NA, dim = c(65, 6, 6))
# not_sampled_rainfall_annual_ratio_summary <- array(data = NA, dim = c(65, 6, 6))

# not_sampled_vector_annual_ratio_summary[i, 1:6, 1:6] <- not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6]/not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6]
# not_sampled_vector_rainfall_ratio_summary[i, 1:6, 1:6] <- not_sampled_vector_peak_absolute_summary[i, 1:6, 1:6]/not_sampled_rainfall_peak_absolute_summary[i, 1:6, 1:6]
# not_sampled_rainfall_annual_ratio_summary[i, 1:6, 1:6] <- not_sampled_rainfall_peak_absolute_summary[i, 1:6, 1:6]/not_sampled_annual_avg_absolute_summary[i, 1:6, 1:6]

