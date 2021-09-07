#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools)

# Load functions
source(here("functions", "functions/time_series_characterisation_functions.R"))
source(here("functions", "gp_fitting_functions.R"))

# Loading in extracted stephensi data and renaming variables
df <- read.csv(file = here("data", "raw", "extracted_data.csv"), stringsAsFactors = FALSE)
df <- df[1:85, ] %>%
  select(Time.Series.ID, Multiply.By., Month_Start, `Jan`:`Dec.4`) %>%
  rename(id = Time.Series.ID, multiply = Multiply.By., start = Month_Start) %>%
  mutate(multiply = as.numeric(multiply))

# Processing entomological data in the following ways:
#   1) Averaging months if a time-series spans >12 months and so some months have multiple values
#   2) Reordering data so that all months start in January
new_df <- matrix(nrow = 85, ncol = 12)
for (i in 1:85) {
  
  # Load in individual time series
  initial_ts <- as.numeric(df[i, 4:63])
  initial_ts <- na.trim(initial_ts)
  if (i == 16) { # use only Year 1988 for Ansari 1990
    initial_ts <- tail(initial_ts, 12)
  } else if (i == 70) {
    initial_ts <- tail(initial_ts, 24)
  }
  length <- length(initial_ts) 
  
  # Month the data starts at 
  start <- df[i, "start"]
  
  # Averaging same months from multiple years 
  if(length > 12) {
    temp <- vector(mode = "logical", length = 12)
    for (j in 1:12) {
      inds <- seq(j, length, 12)
      temp_vals <- initial_ts[inds]
      temp[j] <- mean(temp_vals)
    }
  } else {
    temp <- initial_ts
  }
  
  # Rearranging time series that don't start at January to make 
  #  them begin with the January observation
  if(start != 1) {
    temp <- temp[c((12-start+2):12, 1:(12-start+1))]
  } else{
    padding <- 12 - length(temp)
    temp <- c(temp, rep(NA, padding))
  }
  
  # Multiplying by count conversion factor
  multiply <- df[i, "multiply"]
  new_df[i, ] <- temp * multiply
  
}

#######################################################################################################
##                                                                                                   ##
##                            Negative Binomial Gaussian Process Fitting                             ##
##                                                                                                   ##
##    This next section fits a Negative Binomial Gaussian Process with Periodic Kernel to each of    ##
##    the time series, implemented in the probabilistic programming language STAN.                   ##
##                                                                                                   ##
#######################################################################################################
#options(mc.cores = parallel::detectCores() - 4)
GP_model <- stan_model(here("models", "Neg_Binom_Periodic_Kernel.stan"))
interpolating_points <- 2
prior <- "informative"
counter <- 1
for (i in counter:85) {
  input_time_series <- new_df[i, ]
  fitting_output <- periodic_gp_fit(input_time_series, interpolating_points, i, GP_model, 1000, 1, prior, TRUE)
}

# Visualising all the different time series and removing those with <11 months
sums <- apply(new_df, 1, sum, na.rm = TRUE)
length <- apply(new_df, 1, function(x) sum(!is.na(x)))
remove <- which(sums < 25 | length < 10)
retain_index <- seq(1, 85)[-remove]
for (i in retain_index) {
  mean_realisation_extract(i, new_df, prior, TRUE)
  browser()
}

#######################################################################################################
##                                                                                                   ##
##                               Time-Series Property Characterisation                               ##
##                                                                                                   ##
##    This next section generates a reduced representation of the time-series through                ##
##    statistically characterising their temporal properties using various summary statistics.       ##
##                                                                                                   ##
#######################################################################################################

# Generating the Time Series Features
mean_realisation <- matrix(nrow = length(retain_index), ncol = (12 * interpolating_points + 1))
features <- matrix(nrow = length(retain_index), ncol = 7)
colnames(features) <- c("entropy", "period", "prop_points", "jan_dist", "peaks", "mean", "weight")
for (i in 1:length(retain_index)) {
  
  index <- retain_index[i]
  
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
  
  # Fitting 1 and 2 component Von Mises distributions to the smoothed data
  normalised_output <- normalise_total(negbinom_intensity_mean)
  norm_output <- normalised_output/AUC(2*pi*ordered_timepoints/length(ordered_timepoints), normalised_output) # normalising so AUC sums to 1
  von_mises_params <- von_mises_fitting(norm_output, TRUE)
  diff <- as.numeric(abs(von_mises_params["2_1_Mean"] - von_mises_params["2_2_Mean"]))
  weight_temp <- as.numeric(von_mises_params["2_W"])
  if ((diff > 2.094395 & diff < 4.18879) & (weight_temp > 0.3 & weight_temp < 0.7)) {
    peaks <- 2
  } else {
    peaks <- 1
  }
  if (peaks == 1) {
    mean <- as.numeric(min(von_mises_params["1_Mean"], abs(2*pi - von_mises_params["1_Mean"])))
    weight <- 1
  } else if (peaks == 2) {
    mean <- -5
    weight <- as.numeric(max(von_mises_params["2_W"], 1 - von_mises_params["2_W"]))
  }
  
  # Extracting the features from the time-series
  # -> Note in contrast to Whittaker et al, these are ALL calculated on the fitted curves, rather than raw data
  # -> Need to figure out whether doing all this (including Von Mises fitting) to fitted curves is correct (rather than raw data)
  entropy <- entropic_measure(negbinom_intensity_mean)
  median_period <- median(MCMC_output[, "period"])
  prop_points_mean <- points_greater_than_mean(1.6, negbinom_intensity_mean)
  distance_from_jan <- calculate_peak_distance_from_jan(ordered_timepoints, negbinom_intensity_mean)
  
  # Adding single time-series features to overall features matrix
  features[i, ] <- c(entropy, median_period, prop_points_mean, distance_from_jan, peaks, mean, weight)
  print(i)
}

#######################################################################################################
##                                                                                                   ##
##                                   Clustering and Visualisation                                    ##
##                                                                                                   ##
##    This next section applies a PCA to the time-series feature data, then clusters similar         ##
##    time-series and plots them in the dimensionally reduced space.                                 ##
##                                                                                                   ##
#######################################################################################################

# Normalising features
normalised_features <- scale(features)
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Running PCA on normalised features
PCA <- prcomp(normalised_features)
summary <- summary(PCA)
loadings <- PCA$rotation[, 1:7]
PCA_output <- as.matrix(normalised_features) %*% loadings

# Clustering the data 
num_clust <- 4
clustering_results <- kmeans(PCA_output[, 1:4], num_clust, nstart = 20)

# Visualising the time-series belonging to each cluster
colours <- palette()[1:num_clust]
timepoints <- seq(0, 12, length = dim(normalised_output)[2])
cluster_membership <- clustering_results$cluster
table(cluster_membership)
par(mfrow = c(2, 2))
for (i in 1:num_clust) {
  cluster <- normalised_output[clustering_results$cluster == i, ]
  max <- max(cluster)
  plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", ylim = c(0, 15), lwd = 2, col = colours[i], 
       las = 1, xaxt = "n", xlab = "", ylab = "")
  for (j in 1:length(cluster[, 1])) {
    lines(timepoints, cluster[j, ] * 100, col = adjustcolor(colours[i], alpha.f = 0.2))
  }
  number_time_series <- length(cluster[, 1])
  text(1.5, 9.5, paste0("n = ", number_time_series), cex = 1.5, col = "grey20")
}

# Visualising the clusters of points
plot(PCA_output[, 2], PCA_output[, 1], col = clustering_results$cluster, pch = 20, xlab = "PCA Comp 2", ylab = "PCA Comp 1", cex = 2, las = 1)
clusters <- as.character(seq(1:num_clust))
legend("bottomright", as.character(seq(1:num_clust)), cex = 1.5, col = palette()[1:num_clust], pch = 20)


plot(normalised_output[4, ])




prop_exp <- PCA$sdev^2/sum(PCA$sdev^2)
plot(seq(1, 7, 1), prop_exp, ylab = "% Variance Explained", xlab = "PCA Component", las = 1)
lines(seq(1, 7, 1), prop_exp)
points(seq(1, 7, 1), prop_exp, pch = 20)

par(mfrow = c(1, 1), mar = c(3, 3, 3, 3))
plot(PCA_output[, 2], PCA_output[, 1], pch = 20, xlab = "PCA Comp 2 (20%)", ylab = "PCA Comp 1 (55%)", cex = 2)

