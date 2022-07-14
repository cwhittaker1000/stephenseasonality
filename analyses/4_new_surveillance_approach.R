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

# Reordering
reordered_mean_realisation <- matrix(nrow = 65, ncol = (12 * interpolating_points + 1))
start_index <- apply(mean_realisation, 1, function(x) which(x == max(x)))
start_index_month <- round(start_index/25 * 12)
end_index <- dim(reordered_mean_realisation)[2]
for (i in 1:65) {
  reordered_mean_realisation[i, ] <- mean_realisation[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
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
num_months_sampled <- 1:12 # varying effort 
num_nights_per_month <- 1:10 # varying effort

# For single time series and contiguous sampling (so 12 possible combinations of sampled months all differing
# by starting point)
prob_not_sampled <- array(data = NA, dim = c(65, 12, 10, 12))
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
        prob_not_sampled[t, i, j, k] <- prod(1 - temp_probs)
      }
    }
  }
  print(t)
}

# 1st dim is the time series
# 2nd dim is the number of months sampled
# 3rd dim is the number of nights per month sampled
# 4th dim is which month you start at 
temp <- data.frame(1 - prob_not_sampled[1, 1:6, 1:6, 7])
temp$Y <- paste0("Y", 1:6)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:6))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(y$col_nums)),
                   position = "right",
                   labels = 1:10) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")

plot(reordered_mean_realisation[6, ])
lines(mean_realisation[1, ])

peak_month <- start_index_month[6] # for time series 1, peak month is month 6

temp <- data.frame(1 - prob_not_sampled[6, 1:6, 1:6, peak_month])
temp$Y <- paste0("Y", 1:6)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:6))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(y$col_nums)),
    position = "right",
    labels = 1:10) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")


temp2 <- data.frame(1 - apply(prob_not_sampled[6, 1:6, 1:6, 1:12], c(1, 2), mean)) # average for the year
temp2$Y <- paste0("Y", 1:6)
temp2$Y <- factor(temp2$Y, levels = paste0("Y", 1:6))
temp_long2 <- temp2 %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long2$X  <- factor(temp_long2$X, levels = colnames(temp2)[-length(colnames(temp2))])
ggplot(temp_long2) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(y$col_nums)),
    position = "right",
    labels = 1:10) +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels =  1:12) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled")

temp[, -7]/temp2[, -7]
