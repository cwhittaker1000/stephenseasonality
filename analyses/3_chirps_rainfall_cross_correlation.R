#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools)

#######################################################################################################
##                                                                                                   ##
##                           Raw Monthly Time Series and Rainfall Data                               ##
##                                                                                                   ##
#######################################################################################################
# Load metadata and processed (but non-smoothed) monthly mosquito catch data
metadata <- readRDS(here("data", "processed", "metadata_and_processed_counts.rds"))
processed_metadata <- metadata %>%
  select(id, city, Jan:Dec)  %>%
  tidyr::pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "catch")
cluster <- readRDS(here("data", "processed", "cluster_membership.rds")) %>%
  select(id, cluster)

# Generating monthly rainfall (for comparison to raw, unsmoothed catch data)
monthly_sum_rainfall_storage <- matrix(nrow = nrow(metadata), ncol = 12)
rainfall_files <- list.files(here("data", "processed", "location_specific_rainfall"))
id <- metadata$id
for (i in 1:nrow(metadata)) {
  temp_rainfall <- read.csv(paste0(here("data", "processed", "location_specific_rainfall"), "/rainfall_ts", id[i], ".csv"), stringsAsFactors = FALSE)
  temp_rainfall$month <- format(as.Date(temp_rainfall$day), "%m")
  temp_rainfall <- temp_rainfall %>%
    group_by(month) %>%
    summarise(rainfall = sum(rainfall, na.rm = TRUE))
  monthly_sum_rainfall_storage[i, ] <- temp_rainfall$rainfall
  print(i)
}
ccf_monthly_sum_rainfall_storage <- cbind(metadata$id, monthly_sum_rainfall_storage)
colnames(ccf_monthly_sum_rainfall_storage) <- c("id", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
ccf_monthly_sum_rainfall_storage <- as.data.frame(ccf_monthly_sum_rainfall_storage) %>%
  tidyr::pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "rainfall")

# Combining dataframes together
overall <- processed_metadata %>%
  left_join(ccf_monthly_sum_rainfall_storage, by = c("id", "month")) %>%
  left_join(cluster, by = "id")

cross_cor <- overall %>%
  drop_na(catch) %>%
  group_by(id) %>%
  summarise(cross = ccf(catch, rainfall, lag.max = 0, plot = FALSE)$acf)

cross <- cross_cor %>%
  left_join(cluster, by = "id") %>%
  group_by(cluster) %>%
  summarise(n = n(), mean = mean(cross))

#######################################################################################################
##                                                                                                   ##
##                       Smoothed Time Series and Rainfall Data Comparison                           ##
##                                                                                                   ##
#######################################################################################################

# Loading in smoothed time-series
example <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", 1, ".rds"))
length <- length(example$all_timepoints)
smoothed_mean <- matrix(nrow = length(metadata$id), ncol = length)
smoothed_lower <- matrix(nrow = length(metadata$id), ncol = length)
smoothed_upper <- matrix(nrow = length(metadata$id), ncol = length)
counter <- 1
prior <- "informative"
for (i in metadata$id) {
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", i, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", i, ".rds"))
  }  
  
  # Extracting the mean fitted time-series
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  ordered_timepoints <- all_timepoints[order(all_timepoints)]
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  f_lower <- apply(f, 2, quantile, 0.25)
  f_upper <- apply(f, 2, quantile, 0.75)
  
  mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  lower <- as.numeric(exp(f_lower)[order(all_timepoints)])
  upper <- as.numeric(exp(f_upper)[order(all_timepoints)])
  
  smoothed_mean[counter, ] <- mean
  smoothed_lower[counter, ] <- lower
  smoothed_upper[counter, ] <- upper
  
  counter <- counter + 1
  print(i)
  
}

ccf_smoothed <- smoothed_mean[, -1]
ccf_smoothed <- data.frame(id = cluster$id, cluster = cluster$cluster, ccf_smoothed) %>%
  pivot_longer(X1:X24, names_to = "point", values_to = "catch")

# Generating summed rainfall over specific time-periods (for comparison to smoothed catch data)
sum_rainfall_storage <- matrix(nrow = length(metadata$id), ncol = length-1)
rainfall_files <- list.files(here("data", "processed", "location_specific_rainfall"))
if (length == 25) {
  interpolation_points <- 2
  months_length <- c(15, 16, 14, 14, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16)
  leap_months_length <- c(15, 16, 14, 15, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16)
} else if (length == 37) {
  interpolation_points <- 3
} else {
  stop("how many interpolating points???")
}
id <- metadata$id
for (i in 1:length(metadata$id)) {
  temp_rainfall <- read.csv(paste0(here("data", "processed", "location_specific_rainfall"), "/rainfall_ts", id[i], ".csv"), stringsAsFactors = FALSE)
  temp_summed <- c()
  rainfall <- temp_rainfall$rainfall
  leap <- ifelse(length(rainfall == 365), FALSE, TRUE)
  if(leap) {
    subset <- cumsum(leap_months_length)
  } else {
    subset <- cumsum(months_length)
  }
  for (j in 1:length(months_length)) {
    if (j == 1) {
      temp_summed <- c(temp_summed, sum(rainfall[1:subset[j]]))
    } else {
      temp_summed <- c(temp_summed, sum(rainfall[(subset[j-1]+1):subset[j]]))
    }
  }
  sum_rainfall_storage[i, ] <- temp_summed
  print(i)
}

ccf_sum_rainfall_storage <- data.frame(id = cluster$id, cluster = cluster$cluster, sum_rainfall_storage) %>%
  pivot_longer(X1:X24, names_to = "point", values_to = "rainfall")

# Combining dataframes together
overall_smoothed <- ccf_smoothed %>%
  left_join(ccf_sum_rainfall_storage, by = c("id", "cluster", "point")) 

cross_cor <- overall_smoothed %>%
  drop_na(catch) %>%
  group_by(id) %>%
  summarise(cross = ccf(catch, rainfall, lag.max = 0, plot = FALSE)$acf)

cross <- cross_cor %>%
  left_join(cluster, by = "id") %>%
  group_by(cluster) %>%
  summarise(n = n(), mean = mean(cross))

# Plotting the results
rain_plot <- data.frame(sum_rainfall_storage)
colnames(rain_plot) <- paste0("X", seq(0.25, 11.75, length.out = length(colnames(rain_plot))))
rain_plot <- data.frame(id = cluster$id, cluster = cluster$cluster, rain_plot) %>%
  pivot_longer(X0.25:X11.75, names_to = "point", values_to = "rainfall")

max_rain <- rain_plot %>%
  group_by(id) %>%
  summarise(max_rain = max(rainfall, na.rm = TRUE),
            total_rain = sum(rainfall, na.rm = TRUE))

smooth_mean_plot <- data.frame(smoothed_mean)
colnames(smooth_mean_plot) <- paste0("X", seq(0, 12, length.out = length(colnames(smooth_mean_plot))))
smooth_mean_plot <- data.frame(id = cluster$id, cluster = cluster$cluster, smooth_mean_plot) %>%
  pivot_longer(X0:X12, names_to = "point", values_to = "mean_fitted_catch")

smooth_lower_plot <- data.frame(smoothed_lower)
colnames(smooth_lower_plot) <- paste0("X", seq(0, 12, length.out = length(colnames(smooth_lower_plot))))
smooth_lower_plot <- data.frame(id = cluster$id, cluster = cluster$cluster, smooth_lower_plot) %>%
  pivot_longer(X0:X12, names_to = "point", values_to = "lower_fitted_catch")

smooth_upper_plot <- data.frame(smoothed_upper)
colnames(smooth_upper_plot) <- paste0("X", seq(0, 12, length.out = length(colnames(smooth_upper_plot))))
smooth_upper_plot <- data.frame(id = cluster$id, cluster = cluster$cluster, smooth_upper_plot) %>%
  pivot_longer(X0:X12, names_to = "point", values_to = "upper_fitted_catch")

raw_plot <- metadata %>% select(Jan:Dec)
colnames(raw_plot) <- paste0("X", seq(0.5, 11.5, length.out = length(colnames(raw_plot))))
raw_plot <- data.frame(id = cluster$id, cluster = cluster$cluster, raw_plot) %>%
  pivot_longer(X0.5:X11.5, names_to = "point", values_to = "raw_catch")

overall_plot <- smooth_mean_plot %>%
  left_join(smooth_lower_plot, by = c("id", "cluster", "point")) %>%
  left_join(smooth_upper_plot, by = c("id", "cluster", "point")) %>%
  left_join(raw_plot, by = c("id", "cluster", "point")) %>%
  full_join(rain_plot, by = c("id", "cluster", "point")) 
overall_plot$point <- as.numeric(gsub("X", "", overall_plot$point))

total_raw_catch <- overall_plot %>%
  group_by(id) %>%
  summarise(total_raw_catch = sum(raw_catch, na.rm = TRUE),
            max_raw_catch = max(raw_catch/total_raw_catch, na.rm = TRUE))

overall_plot <- overall_plot %>%
  left_join(total_raw_catch, by = "id") %>%
  mutate(mean_norm = mean_fitted_catch/total_raw_catch,
         lower_norm = lower_fitted_catch/total_raw_catch,
         upper_norm = upper_fitted_catch/total_raw_catch,
         raw_norm = raw_catch/total_raw_catch)

max_record <- overall_plot %>%
  group_by(id) %>%
  summarise(max_val = max(upper_norm, raw_norm, na.rm = TRUE),
            )

overall_plot <- overall_plot %>%
  left_join(max_rain, by = "id") %>%
  left_join(max_record, by = "id") %>%
  mutate(norm_rain = rainfall/max_rain,
         scaled_rain = norm_rain * max_val)

city <- data.frame(id = metadata$id, city = metadata$city)

overall_plot <- overall_plot %>%
  left_join(city, by = "id")

# Example Plot
x <- overall_plot[overall_plot$id %in% c(1:100), ]
y <- overall_plot[overall_plot$id %in% c(1:100) & !is.na(overall_plot$mean_fitted_catch), ]

ggplot() +
  geom_bar(data = x, aes(x = point, y = scaled_rain), stat = "identity", col = adjustcolor("#D1E6F0", alpha.f = 1), fill = "#D1E6F0", alpha = 1, width = 0.49) +
  geom_point(data = x, aes(x = point, y = raw_norm), colour = "black", size = 2) +
  geom_ribbon(data = y, aes(x = point, ymin = lower_norm, 
                            ymax = upper_norm), fill = "#E0521A", alpha = 0.2) +
  geom_line(data = y, aes(x = point, y = mean_norm), col = "#E0521A", size = 2) +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                     breaks = seq(0, 11, length.out = 12)) +
  facet_wrap(~id, ncol = 12, scales = "free") +
  labs(y = "Normalised Value", x= "") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggplot() +
  geom_bar(data = x, aes(x = point, y = scaled_rain), stat = "identity", col = adjustcolor("#D1E6F0", alpha.f = 1), fill = "#D1E6F0", alpha = 1, width = 0.49) +
  geom_point(data = x, aes(x = point, y = raw_norm), colour = "black", size = 2) +
  geom_ribbon(data = y, aes(x = point, ymin = lower_norm, 
                            ymax = upper_norm, fill = city), alpha = 0.2) +
  geom_line(data = y, aes(x = point, y = mean_norm, col = city), size = 2) +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                     breaks = seq(0, 11, length.out = 12)) +
  facet_wrap(~id, ncol = 12, scales = "free") +
  labs(y = "Normalised Value", x= "") +
  theme_bw() +
  theme(panel.grid = element_blank())

# #6EB4D1 another alternative blue
# geom_bar(data = x, aes(x = point, y = 2 * rainfall/sum(rainfall, na.rm = TRUE)), stat = "identity", 
#          col = adjustcolor("#D1E6F0", alpha.f = 1), fill = "#D1E6F0", alpha = 1, width = 0.5) 