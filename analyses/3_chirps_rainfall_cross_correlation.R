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
metadata <- readRDS(here("data", "processed", "metadata_and_processed_counts.rds")) %>%
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
monthly_sum_rainfall_storage <- cbind(metadata$id, monthly_sum_rainfall_storage)
colnames(monthly_sum_rainfall_storage) <- c("id", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
monthly_sum_rainfall_storage <- as.data.frame(monthly_sum_rainfall_storage) %>%
  tidyr::pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "rainfall")

# Combining dataframes together
overall <- metadata %>%
  left_join(monthly_sum_rainfall_storage, by = c("id", "month")) %>%
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
smoothed <- matrix(nrow = length(unique(metadata$id)), ncol = length)
counter <- 1
prior <- "informative"
for (i in unique(metadata$id)) {
  
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
  negbinom_intensity_mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  smoothed[counter, ] <- negbinom_intensity_mean
  counter <- counter + 1
  print(i)
  
}

smoothed <- smoothed[, -1]
smoothed <- data.frame(id = cluster$id, cluster = cluster$cluster, smoothed) %>%
  pivot_longer(X1:X24, names_to = "point", values_to = "catch")

# Generating summed rainfall over specific time-periods (for comparison to smoothed catch data)
sum_rainfall_storage <- matrix(nrow = length(unique(metadata$id)), ncol = length-1)
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
id <- unique(metadata$id)
for (i in 1:length(unique(metadata$id))) {
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

sum_rainfall_storage <- data.frame(id = cluster$id, cluster = cluster$cluster, sum_rainfall_storage) %>%
  pivot_longer(X1:X24, names_to = "point", values_to = "rainfall")

# Combining dataframes together
overall_smoothed <- smoothed %>%
  left_join(sum_rainfall_storage, by = c("id", "cluster", "point")) 

cross_cor <- overall_smoothed %>%
  drop_na(catch) %>%
  group_by(id) %>%
  summarise(cross = ccf(catch, rainfall, lag.max = 0, plot = FALSE)$acf)

cross <- cross_cor %>%
  left_join(cluster, by = "id") %>%
  group_by(cluster) %>%
  summarise(n = n(), mean = mean(cross))


monthly_sum_rainfall_storage <- cbind(metadata$id, monthly_sum_rainfall_storage)
colnames(monthly_sum_rainfall_storage) <- c("id", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
monthly_sum_rainfall_storage <- as.data.frame(monthly_sum_rainfall_storage) %>%
  tidyr::pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "rainfall")

ento_data <- metadata %>%
  dplyr::select(id, city, Jan:Dec) %>%
  pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "catch")

overall <- ento_data %>%
  left_join(monthly_sum_rainfall_storage, by = c("id", "month")) 

months <- overall$month

overall <- overall %>%
  group_by(id) %>%
  summarise(catch = catch/sum(catch),
            rainfall = rainfall/sum(rainfall))
overall$month <- months
overall$month <- factor(overall$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

par(mfrow = c(12, 6))
ggplot(data = overall[overall$id == 1, ]) +
  geom_point(aes(x = month, y = catch), colour = "#E0521A", size = 2) +
  geom_line(aes(x = month, y = catch, group = 1)) +
  geom_bar(aes(x = month, y = rainfall), stat = "identity", fill = "#6EB4D1", alpha = 0.5) +
  facet_wrap(~id, scales = "free") +
  scale_x_discrete(name = unique(new_overall$name),
                   labels = c("J", "", "", "F", "", "", "M", "", "", 
                              "A", "", "", "M", "", "", "J", "", "",
                              "J", "", "", "A", "", "", "S", "", "",
                              "O", "", "", "N", "", "", "D", "", ""))
