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
source(here("functions", "time_series_characterisation_functions.R"))
source(here("functions", "gp_fitting_functions.R"))

# Loading in extracted stephensi data and renaming variables
raw_df <- read.csv(file = here("data", "systematic_review_results", "extracted_entomological_data.csv"), stringsAsFactors = FALSE)
df <- raw_df[1:85, ] %>%
  dplyr::select(Time.Series.ID, Multiply.By., Month_Start, `Jan`:`Dec.4`) %>%
  rename(id = Time.Series.ID, multiply = Multiply.By., start = Month_Start) %>%
  mutate(multiply = as.numeric(multiply))

# Processing entomological data in the following ways:
#   1) Averaging months if a time-series spans >12 months and so some months have multiple values
#   2) Reordering data so that all months start in January
new_df <- matrix(nrow = 85, ncol = 12)
for (i in 1:85) {
  
  # Load in individual time series
  initial_ts <- as.numeric(df[i, 4:(dim(df)[2])])
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
set.seed(10)
fresh_run <- FALSE
prior <- "informative"
interpolating_points <- 2
if (fresh_run) {
  GP_model <- stan_model(here("models", "Neg_Binom_Periodic_Kernel.stan"))
  counter <- 1
  for (i in counter:85) {
    input_time_series <- new_df[i, ]
    fitting_output <- periodic_gp_fit(input_time_series, interpolating_points, i, GP_model, 1000, 1, prior, TRUE)
  }
}

# Visualising all the different time series and removing those with <11 months
sums <- apply(new_df, 1, sum, na.rm = TRUE)
length <- apply(new_df, 1, function(x) sum(!is.na(x)))
remove <- which(sums < 25 | length < 10)
retain_index <- seq(1, 85)[-remove]
for (i in retain_index) {
  mean_realisation_extract(i, new_df, prior, TRUE)
  #browser()
}

# Saving the processed (to make 12 months) but unsmoothed results + relevant metadata
metadata <- raw_df[retain_index, ] %>%
  dplyr::select(Time.Series.ID, Country, Admin.1, Admin.2, City., Year.Start, Year.End) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2, city = City., 
         start = Year.Start, end = Year.End)
overall <- cbind(metadata, new_df[retain_index, ])
colnames(overall)[8:19] <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
#saveRDS(overall, file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
table(metadata$country)  
table(metadata$city)  

#######################################################################################################
##                                                                                                   ##
##                               Time-Series Property Characterisation                               ##
##                                                                                                   ##
##    This next section generates a reduced representation of the time-series through                ##
##    statistically characterising their temporal properties using various summary statistics.       ##
##                                                                                                   ##
#######################################################################################################

# Extracting Mean Realisations
set.seed(10)
mean_realisation <- matrix(nrow = length(retain_index), ncol = (12 * interpolating_points + 1))
period <- vector(mode = "numeric", length = length(retain_index))
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
  period[i] <- median(MCMC_output[, "period"])

}
  
# Reordering mean fitted time-series to start at max
reordered_mean_realisation <- matrix(nrow = length(retain_index), ncol = (12 * interpolating_points + 1))
start_index <- apply(mean_realisation, 1, function(x) which(x == max(x)))
end_index <- dim(reordered_mean_realisation)[2]
for (i in 1:length(retain_index)) {
  reordered_mean_realisation[i, ] <- mean_realisation[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
}


# Extracting Rainfall Data (also sort out leap year stuff)
leap_years <- c(1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2016,	2020)
months_length <- c(15, 16, 14, 14, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16)
rainfall_storage <- matrix(nrow = dim(overall)[1], ncol = length(months_length)) # this is 24 long rather than 25 like ento fits - why is this?? Need to check
for (i in 1:length(overall$id)) {
  index <- overall$id[i]
  temp <- c()
  rf <- read.csv(paste0("data/location_specific_rainfall/rainfall_ts", index, ".csv"))
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
norm_rainfall_storage <- normalise_total(rainfall_storage)
smoothed_rainfall <- t(apply(norm_rainfall_storage, 1, raster::movingFun, n = 4, circular = TRUE))
rainfall_start_index <- apply(smoothed_rainfall, 1, function(x) which(x == max(x)))
rainfall_seas_3 <- apply(norm_rainfall_storage, 1, percent_incidence, 3, 2)
rainfall_seas_4 <- apply(norm_rainfall_storage, 1, percent_incidence, 4, 2)

# Generating Features
set.seed(10)
features <- matrix(nrow = length(retain_index), ncol = 10)
per_ind_3_months_vector <- c()
colnames(features) <- c("dens_peak_timing", "rainfall_peak_timing", "peak_diff", "entropy", "period", "prop_points", "peaks", "mean", "weight", "per_ind_4_months")
for (i in 1:length(retain_index)) {

  # Loading in Specific Time Series
  negbinom_intensity_mean <- reordered_mean_realisation[i, ]
  
  # Fitting 1 and 2 component Von Mises distributions to the smoothed data
  normed_output <- normalise_total(negbinom_intensity_mean)
  norm_output <- normed_output/AUC(2*pi*ordered_timepoints/length(ordered_timepoints), normed_output) # normalising so AUC sums to 1
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
  dens_peak <- start_index[i]
  rain_peak <- rainfall_start_index[i]
  peak_diff <- dens_peak - rain_peak
  entropy <- entropic_measure(negbinom_intensity_mean)
  prop_points_mean <- points_greater_than_mean(1.6, negbinom_intensity_mean)
  percent_incidence_out <- percent_incidence(negbinom_intensity_mean, 4)
  per_ind_3_months_vector <- c(per_ind_3_months_vector, percent_incidence(negbinom_intensity_mean, 3))
  
  # Adding single time-series features to overall features matrix
  features[i, ] <- c(dens_peak, rain_peak, peak_diff, entropy, period[i], prop_points_mean, peaks, mean, weight, percent_incidence_out)
  print(i)
}

features_df <- data.frame(id = metadata$id, country = metadata$country, admin1 = metadata$admin1, admin2 = metadata$admin2,
                          city = metadata$city, features, per_ind_3_months = per_ind_3_months_vector,
                          rainfall_seas_3 = rainfall_seas_3, rainfall_seas_4 = rainfall_seas_4)
saveRDS(features_df, file = here("data", "systematic_review_results", "metadata_and_rearranged_time_series_features.rds"))

hist(features_df$per_ind_4_months, breaks = 8)

which(features_df$per_ind_4_months == min(features_df$per_ind_4_months))
which(features_df$per_ind_4_months == max(features_df$per_ind_4_months))

plot(reordered_mean_realisation[11, ], type = "l", ylim = c(0, max(reordered_mean_realisation[11, ])))
plot(reordered_mean_realisation[35, ], type = "l", ylim = c(0, max(reordered_mean_realisation[35, ])))

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
normalised_features <- normalised_features[, -c(1, 2)]
normalised_output <- t(apply(reordered_mean_realisation, 1, normalise_total))

# Running PCA on normalised features
PCA <- prcomp(normalised_features)
summary <- summary(PCA)
loadings <- PCA$rotation[, 1:8]
PCA_output <- as.matrix(normalised_features) %*% loadings

# Clustering the data 
num_clust <- 2
clustering_results <- kmeans(PCA_output[, 1:4], num_clust, nstart = 20)

cluster_output <- data.frame(id = metadata$id, country = metadata$country, city = metadata$city, 
                             cluster = clustering_results$cluster, reordered_mean_realisation)
saveRDS(cluster_output, file = here("data", "systematic_review_results", "rearranged_cluster_membership.rds"))

# Visualising the time-series belonging to each cluster
colours <- palette()[1:num_clust]
timepoints <- seq(0, 12, length = dim(normalised_output)[2])
cluster_membership <- clustering_results$cluster
table(cluster_membership)
par(mfrow = c(2, 2))
for (i in 1:num_clust) {
  cluster <- normalised_output[clustering_results$cluster == i, ]
  max <- max(cluster)
  plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", ylim = c(0, 20), lwd = 2, col = colours[i], 
       las = 1, xaxt = "n", xlab = "", ylab = "")
  for (j in 1:length(cluster[, 1])) {
    lines(timepoints, cluster[j, ] * 100, col = adjustcolor(colours[i], alpha.f = 0.2))
  }
  number_time_series <- length(cluster[, 1])
  text(1.5, 9.5, paste0("n = ", number_time_series), cex = 1.5, col = "grey20")
}

# Visualising the clusters of points
one_cov <- cov(PCA_output[cluster_membership == 1, 1:2])
one_mean <- apply(PCA_output[cluster_membership == 1, 1:2], 2, mean)
one_ellipse <- ellipse::ellipse(x = one_cov, centre = one_mean, level = 0.95)
two_cov <- cov(PCA_output[cluster_membership == 2, 1:2])
two_mean <- apply(PCA_output[cluster_membership == 2, 1:2], 2, mean)
two_ellipse <- ellipse::ellipse(x = two_cov, centre = two_mean, level = 0.95)
ellipse_df <- data.frame(cluster = factor(c(rep(1, length(one_ellipse[, 1])),
                                            rep(2, length(two_ellipse[, 1])))),
                         rbind(one_ellipse, two_ellipse))
pca_df <- data.frame(cluster = factor(cluster_membership), PCA_output)

pca_plot <- ggplot() +
  geom_polygon(data = ellipse_df, aes(x = PC2, y = PC1, fill = cluster), alpha = 0.2) +
  geom_point(data = pca_df, aes(x = PC2, y = PC1, col = cluster), size = 2) +
  coord_cartesian(xlim = c(-2.25, 5), ylim = c(-2.5, 8)) + 
  scale_fill_manual(values = colours[1:2]) +
  scale_colour_manual(values = colours[1:2]) +
  labs(x = "Principal Component 2 (15% Total Variation)",
       y = "Principal Component 1 (54% Total Variation)") +
  theme_bw() +
  theme(legend.position = "none")

pca_plot +
  scale_fill_manual(values = palette()[3:4]) +
  scale_colour_manual(values = palette()[3:4]) +
  labs(x = "PC2", y = "PC1") 

cluster_catch_seasonality +
  scale_fill_manual(values = palette()[3:4]) +
  scale_colour_manual(values = palette()[3:4])

## PCA Cluster Results Plotting
cluster_df <- data.frame(cluster = factor(3 - cluster_membership), id = seq(1:dim(normalised_output)[1]), normalised_output)
cluster_df <- cluster_df %>%
  pivot_longer(cols= X1:X25, names_to = "timepoint", values_to = "density") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint))) %>%
  mutate(timepoint2 = ifelse(timepoint <= 12, timepoint + 13, timepoint - 12))
mean_df <- cluster_df %>%
  group_by(timepoint, timepoint2, cluster) %>%
  summarise(mean = mean(density))
cluster_time_series <- ggplot(cluster_df, aes(x = timepoint2, y = 100 * density)) +
  geom_line(aes(col = cluster, group = id), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 20)) +
  #scale_x_continuous(breaks = seq(1, 23, 2), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_x_continuous(breaks = seq(1, 23, 2), labels = 1:12) +
  facet_wrap(~cluster, nrow = 2) +
  theme_bw() +
  geom_line(data = mean_df, aes(x = timepoint2, y = 100 * mean, col = cluster), size = 2) +
  scale_color_manual(values = colours[2:1]) +
  labs(y = "Normalised Monthly Vector Density",
       x = "Peak-Standardised Time (Months)") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

set.seed(16)
z <- data.frame(id = c(rep("dens", 65), rep("rain", 65)),
                per_ind = c(features_df$per_ind_4_months, features_df$rainfall_seas_4),
                cluster = c(3-cluster_membership, 3-cluster_membership))
cluster_catch_seasonality <- ggplot(z, aes(x = factor(cluster), y = 100 * per_ind, col = factor(cluster)))  +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  facet_wrap(~id, nrow = 2,
             strip.position = "right",
             labeller = as_labeller(c(dens = "% Annual Density In 4 Months", rain = "% Annual Rainfall In 4 Months"))) +
  scale_color_manual(values = colours[2:1]) +
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2")) +
  geom_jitter(aes(x = factor(cluster), y = 100 * per_ind), size = 1, width = 0.25) +
  theme_bw() +
  scale_y_continuous(position = "right", limits = c(35, 100)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside")

al <- cowplot::plot_grid(cluster_time_series, cluster_catch_seasonality, ncol = 2, rel_widths = c(2, 1),
                   axis = "b", align = "h")

ab <- cowplot::plot_grid(pca_plot, al, ncol = 2, rel_widths = c(2, 2))
#ggsave width = 10, height = 5

test <- 

ggsave(cluster_catch_seasonality, file = "textboxplot.pdf", width = 2, height = 4)

mean(features_df$per_ind_3_months[cluster_membership == 2])
mean(features_df$per_ind_3_months[cluster_membership == 1])

table(features_df$city, cluster_membership)
chisq.test(x = features_df$city[!features_df$city == "Mixture/Unclear"],
           y = cluster_membership[!features_df$city == "Mixture/Unclear"],
           simulate.p.value = TRUE, B = 5000)



# Visualising the properties of each cluster
number_properties <- ncol(normalised_features)
property_names <- colnames(normalised_features)
par(mfrow = c(num_clust, number_properties + 1),
    mar = c(4.5, 2, 2, 2))
max <- 10 
mean_operation_values <- matrix(nrow = num_clust, ncol = number_properties)
colnames(mean_operation_values) <- property_names
for (i in 1:num_clust) {
  if (i < num_clust) {
    subsetter <- clustering_results$cluster == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2, ylab = "", xlab = "")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlab = "", ylab = "",
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  } else {
    subsetter <- clustering_results$cluster == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2, ylab = "", xlab = "")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlab = "", ylab = "",
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
      mtext(property_names[j], side = 1, outer = FALSE, cex = 1, font = 2, line = 3, col = "grey20")
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  }
}

# Visualising the time-series belonging to urban/rural
urban_rural <- overall$city
tab <- table(clustering_results$cluster, urban_rural)
tab <- tab[, 2:3]
chisq.test(tab)
x <- table(features[, "peaks"], urban_rural)
chisq.test(x[1:2, 2:3])
colours <- palette()[1:2]
timepoints <- seq(0, 12, length = dim(normalised_output)[2])
locs <- c("Urban", "Rural")
par(mfrow = c(1, 2))
for (i in 1:2) {
  set <- normalised_output[urban_rural == locs[i], ]
  max <- max(set)
  plot(timepoints, apply(set, 2, mean) * 100, type = "l", ylim = c(0, 10), lwd = 2, col = colours[i], 
       las = 1, xaxt = "n", xlab = "", ylab = "")
  for (j in 1:length(set[, 1])) {
    lines(timepoints, set[j, ] * 100, col = adjustcolor(colours[i], alpha.f = 0.2))
  }
  number_time_series <- length(set[, 1])
  text(1.5, 9.5, paste0("n = ", number_time_series), cex = 1.5, col = "grey20")
}

mean(features$per_ind_4_months[urban_rural == "Urban"]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Rural"]) # mixture of unimodal and bimodal time series for rural

mean(features$peaks[urban_rural == "Urban"]) # almost all urban time series are unimodal
mean(features$peaks[urban_rural == "Rural"]) # mixture of unimodal and bimodal time series for rural

mean(features$per_ind_4_months[clustering_results$cluster == 1])
mean(features$per_ind_4_months[clustering_results$cluster == 2])
mean(features$per_ind_4_months[clustering_results$cluster == 3])
mean(features$per_ind_4_months[clustering_results$cluster == 4])

mean(features$jan_dist[clustering_results$cluster == 1])
mean(features$jan_dist[clustering_results$cluster == 2])
mean(features$jan_dist[clustering_results$cluster == 3])
mean(features$jan_dist[clustering_results$cluster == 4])

# Exploring entropy in further detail - consider removing/check with Sam I'm doing it correct
# par(mfrow = c(2, 4))
# for (i in which(cluster_membership == 3)) {
#   index <- metadata$id[i]
#   mean_realisation_extract(index, new_df, prior, TRUE)
#   browser()
# }
# features[which(cluster_membership == 3), ]
# 
# interval <-  findInterval(features[, "entropy"], quantile(features[, "entropy"], probs=0:5/5))
# colours <- palette()[1:length(unique(interval))]
# timepoints <- seq(0, 12, length = dim(normalised_output)[2])
# cluster_membership <- interval
# table(cluster_membership)
# par(mfrow = c(2, 3))
# for (i in 1:(length(unique(interval)) - 1)) {
#   cluster <- normalised_output[cluster_membership == i, ]
#   max <- max(cluster)
#   plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", ylim = c(0, 10), lwd = 2, col = colours[i], 
#        las = 1, xaxt = "n", xlab = "", ylab = "")
#   for (j in 1:length(cluster[, 1])) {
#     lines(timepoints, cluster[j, ] * 100, col = adjustcolor(colours[i], alpha.f = 0.2))
#   }
#   number_time_series <- length(cluster[, 1])
#   entropy <- mean(features[cluster_membership == i, "entropy"])
#   text(1.5, 9.5, paste0("n = ", number_time_series), cex = 1.5, col = "grey20")
#   text(6.5, 9.5, paste0("Entropy = ", round(entropy, 0)), cex = 1.5, col = "grey20")
# }
