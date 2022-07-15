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
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

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
# for (i in retain_index) {
#   mean_realisation_extract(i, new_df, prior, TRUE)
  #browser()
#}

# Saving the processed (to make 12 months) but unsmoothed results + relevant metadata
metadata <- raw_df[retain_index, ] %>%
  dplyr::select(Time.Series.ID, Country, Admin.1, Admin.2, City., Year.Start, Year.End) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2, city = City., 
         start = Year.Start, end = Year.End)
overall <- cbind(metadata, new_df[retain_index, ])
colnames(overall)[8:19] <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
if (fresh_run) {
  saveRDS(overall, file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
}

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
lower_realisation <- matrix(nrow = dim(metadata)[1], ncol = (12 * interpolating_points + 1))
upper_realisation <- matrix(nrow = dim(metadata)[1], ncol = (12 * interpolating_points + 1))
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
  f_lower <- apply(f, 2, quantile, 0.125)
  f_upper <- apply(f, 2, quantile, 0.875)
  negbinom_intensity_mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  lower <- as.numeric(exp(f_lower)[order(all_timepoints)])
  upper <- as.numeric(exp(f_upper)[order(all_timepoints)])
  mean_realisation[i, ] <- negbinom_intensity_mean
  lower_realisation[i, ] <- lower
  upper_realisation[i, ] <- upper
  period[i] <- median(MCMC_output[, "period"])

}

# Plotting all the Fitted Time-Series
mean <- as.data.frame(mean_realisation)
mean$country <- metadata$country
mean$id <- retain_index
mean_pv <- mean %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "mean")

lower <- as.data.frame(lower_realisation)
lower$country <- metadata$country
lower$id <- retain_index
lower_pv <- lower %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "lower")

upper <- as.data.frame(upper_realisation)
upper$country <- metadata$country
upper$id <- retain_index
upper_pv <- upper %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "upper")

catch_data <- overall %>%
  filter(id %in% retain_index) %>%
  dplyr::select(id, country, Jan:Dec)
colnames(catch_data) <- c("id", "country", seq(2, 24, length.out = 12)) 
total_raw_catch <- catch_data %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "raw_catch") %>%
  group_by(id) %>%
  summarise(total_raw = sum(raw_catch, na.rm = TRUE))
raw_catch_data <- catch_data %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "raw_catch") %>%
  left_join(total_raw_catch, by = "id") %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  dplyr::select(country, id, timepoint, raw_catch)

overall_pv <- mean_pv %>%
  left_join(lower_pv, by = c("id", "country", "timepoint")) %>%
  left_join(upper_pv, by = c("id", "country", "timepoint")) %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  left_join(raw_catch_data, by =  c("id", "country", "timepoint"))

all_ts_plot <- ggplot(data = overall_pv) +
  geom_path(aes(x = timepoint, y = mean, col = factor(country)), size = 2) +
  geom_point(aes(x = timepoint, y = raw_catch, group = factor(id)), col = "black") +
  geom_ribbon(aes(x = timepoint, ymin = lower, ymax = upper, fill = country), alpha = 0.2) +
  facet_wrap(~id, scales = "free_y") +
  scale_y_continuous(limits=c(0, NA), position = "left") +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(2, 24, length.out = 12)) +
  scale_fill_manual(values = gg_color_hue(6)) +
  ylab("Monthly Catch") +
  theme_bw() +
  theme(#legend.position = "bottom",
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.justification = c(1, 0), legend.position = c(1, 0)) +
  guides(col = "none", fill = guide_legend(nrow = 1, ncol = 6, title = ""))
if (fresh_run) {
  ggsave2(file = here("figures", "Supp_Fig_AllTs.pdf"), plot = all_ts_plot, width = 16, height = 7, dpi = 500)
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
norm_rainfall_storage <- normalise_total(rainfall_storage)
smoothed_rainfall <- t(apply(norm_rainfall_storage, 1, raster::movingFun, n = 4, circular = TRUE))
rainfall_start_index <- apply(smoothed_rainfall, 1, function(x) which(x == max(x)))
rainfall_seas_3 <- apply(norm_rainfall_storage, 1, percent_incidence, 3, 2)
rainfall_seas_4 <- apply(norm_rainfall_storage, 1, percent_incidence, 4, 2)

# Reordering rainfall time-series to start at max of vector density (note NOT max rainfall, which is what I calculate above)
reordered_rainfall <- matrix(nrow = length(retain_index), ncol = length(months_length))
end_rain_index <- dim(reordered_rainfall)[2]
for (i in 1:length(retain_index)) {
  reordered_rainfall[i, ] <- norm_rainfall_storage[i, c(start_index[i]:end_rain_index, 1:(start_index[i]-1))]
}

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
if (fresh_run) {
  saveRDS(features_df, file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
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
normalised_features <- normalised_features[, -c(1, 2)]
normalised_output <- t(apply(reordered_mean_realisation, 1, normalise_total))

# Running PCA on normalised features
set.seed(21)
PCA <- prcomp(normalised_features)
summary <- summary(PCA)
loadings <- PCA$rotation[, 1:8]
PCA_output <- as.matrix(normalised_features) %*% loadings

# Clustering the data 
num_clust <- 2
clustering_results <- kmeans(PCA_output[, 1:4], num_clust, nstart = 20)
cluster_membership <- clustering_results$cluster
cluster_output <- data.frame(id = metadata$id, country = metadata$country, city = metadata$city, cluster = clustering_results$cluster, reordered_mean_realisation)
saveRDS(cluster_output, file = here("data", "systematic_review_results", "cluster_membership.rds"))

# Comparing Sample Sizes ## GO BACK AND CALCULATE MONTHLY CATCHES INSTEAD OF TOTAL AS SOME HAVE
## 11 MONTHS NOT 12
one <- cluster_membership == 1
one_total <- total_raw_catch$total_raw[one]
two <- cluster_membership == 2
two_total <- total_raw_catch$total_raw[two]
median(one_total)
IQR(one_total)
quantile(one_total, probs = c(0.25, 0.75))
median(two_total)
mean(one_total)
mean(two_total)
quantile(two_total, probs = c(0.25, 0.75))

df_catch_size <- data.frame(total_raw_catch, cluster = factor(cluster_membership))
catch_size_hist <- ggplot(df_catch_size, aes(x = log(total_raw))) +
  geom_histogram(aes(fill = cluster)) +
  labs(x = "Log Total Catch Size") +
  scale_fill_manual(values = c("#DF536B", "black")) +
  theme_bw() +
  theme(legend.position = "none")
df_catch_size$cluster <- factor(df_catch_size$cluster)
df_catch_size$cluster_membership <- ifelse(df_catch_size$cluster == 1, "Cluster 1", "Cluster 2")
catch_size_comparison <- ggplot(df_catch_size, aes(x = cluster_membership, y = log(total_raw), colour = factor(cluster_membership))) +
  geom_boxplot(fill = NA) +
  geom_jitter(aes(x = factor(cluster_membership), y = log(total_raw)), size = 1, width = 0.25) +
  labs(x = "", y = "Log Total Catch Size") +
  scale_colour_manual(values = c("#DF536B", "black")) +
  theme_bw() +
  theme(legend.position = "none")
t.test(total_raw ~ cluster_membership, data = df_catch_size)
mood.test(df_catch_size$total_raw[df_catch_size$cluster == 1], df_catch_size$total_raw[df_catch_size$cluster == 2])
catch_size_plots <- cowplot::plot_grid(catch_size_comparison, catch_size_hist, nrow = 1, ncol = 2)
if (fresh_run) {
  cowplot::ggsave2(file = here("figures", "Supp_Fig_CatchSizeComparison.pdf"), plot = catch_size_plots, width = 11, height = 5, dpi = 500)
}

# Plotting the First 2 Principle Components

## Generating Ellipe and Centre for PC1
one_cov <- cov(PCA_output[cluster_membership == 1, 1:2])
one_mean <- apply(PCA_output[cluster_membership == 1, 1:2], 2, mean)
one_ellipse <- ellipse::ellipse(x = one_cov, centre = one_mean, level = 0.95)

## Generating Ellipe and Centre for PC2
two_cov <- cov(PCA_output[cluster_membership == 2, 1:2])
two_mean <- apply(PCA_output[cluster_membership == 2, 1:2], 2, mean)
two_ellipse <- ellipse::ellipse(x = two_cov, centre = two_mean, level = 0.95)

## Plotting the PCA Loadings and Ellipses
# palette is "black"   "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "gray62"
colours <- palette()[1:2]
ellipse_df <- data.frame(cluster = factor(c(rep(1, length(one_ellipse[, 1])), rep(2, length(two_ellipse[, 1])))), rbind(one_ellipse, two_ellipse))
pca_df <- data.frame(cluster = factor(cluster_membership), PCA_output)
pca_plot <- ggplot() +
  geom_polygon(data = ellipse_df, aes(x = PC2, y = PC1, fill = cluster), alpha = 0.2) +
  geom_point(data = pca_df, aes(x = PC2, y = PC1, col = cluster), size = 2) +
  #coord_cartesian(xlim = c(-2.25, 5), ylim = c(-2.5, 8)) + 
  scale_fill_manual(values = colours[2:1]) +
  scale_colour_manual(values = colours[2:1]) +
  labs(x = "Principal Component 2 (15% Total Variation)",
       y = "Principal Component 1 (54% Total Variation)") +
  theme_bw() +
  theme(legend.position = "none")

# Clustered Time-Series Plotting 
cluster_df <- data.frame(cluster = factor(cluster_membership), id = seq(1:dim(normalised_output)[1]), normalised_output)
cluster_df <- cluster_df %>%
  pivot_longer(cols= X1:X25, names_to = "timepoint", values_to = "density") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint))) %>%
  mutate(timepoint2 = ifelse(timepoint <= 12, timepoint + 13, timepoint - 12))
mean_df <- cluster_df %>%
  group_by(timepoint, timepoint2, cluster) %>%
  summarise(mean = mean(density))
cluster_time_series <- ggplot(cluster_df, aes(x = timepoint2, y = 100 * density * 25/12)) +
  geom_line(aes(col = cluster, group = id), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 42)) +
  scale_x_continuous(breaks = seq(1, 23, 2), labels = 1:12) +
  facet_wrap(~cluster, nrow = 2) +
  theme_bw() +
  geom_line(data = mean_df, aes(x = timepoint2, y = 100 * mean * 25/12, col = cluster), size = 2) +
  scale_color_manual(values = colours[2:1]) +
  labs(y = "Normalised Monthly Vector Density", x = "Peak-Standardised Time (Months)") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank())

set.seed(16)
z <- data.frame(id = c(rep("dens", 65), rep("rain", 65)),
                per_ind = c(features_df$per_ind_4_months, features_df$rainfall_seas_4),
                cluster = c(cluster_membership, cluster_membership))
t.test(per_ind ~ factor(cluster), data = z[z$id == "dens", ])
t.test(per_ind ~ factor(cluster), data = z[z$id == "rain", ])

cluster_catch_seasonality <- ggplot(z, aes(x = factor(cluster), y = 100 * per_ind, col = factor(cluster)))  +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  facet_wrap(~id, nrow = 2, strip.position = "right",
             labeller = as_labeller(c(dens = "% Annual Density In 4 Months", rain = "% Annual Rainfall In 4 Months"))) +
  scale_color_manual(values = colours[2:1]) +
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2")) +
  geom_jitter(aes(x = factor(cluster), y = 100 * per_ind), size = 1, width = 0.25) +
  theme_bw() +
  scale_y_continuous(position = "right", limits = c(35, 100)) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")

cluster_boxplots <- cowplot::plot_grid(cluster_time_series, cluster_catch_seasonality, ncol = 2, rel_widths = c(2, 1), axis = "b", align = "h")
fig2_overall <- cowplot::plot_grid(pca_plot, cluster_boxplots, ncol = 2, rel_widths = c(2, 2.2))
if (fresh_run) {
  ggsave(fig2_overall, file = here("figures", "Fig2_Overall_Raw.pdf"), width = 10, height = 5)
}

# Visualising the properties of each cluster
pdf("Figures/Supp_Figure_Cluster_Properties.pdf", height = 5, width = 12, useDingbats = FALSE)
number_properties <- ncol(normalised_features)
property_names <- colnames(normalised_features)
par(mfrow = c(num_clust, number_properties + 1), mar = c(3.5, 2, 2, 1), oma = c(1, 2.5, 5, 1))
max <- 12
timepoints <- seq(0, 12, length = dim(normalised_output)[2])
mean_operation_values <- matrix(nrow = num_clust, ncol = number_properties)
colnames(mean_operation_values) <- c("Peak Timing", "Entropy", "Period", "Prop. Points\n1.6x Mean", 
                                     "Number of\nPeaks", "Von Mises\nMean", "Von Mises\nWeight", "% Incidence\nIn 4 Months")
for (i in 1:num_clust) {
  if (i == 1) {
    subsetter <- cluster_membership == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[2], lwd = 2, ylab = "", xlab = "", las = 1)
    mtext("Time (Months)", side = 1, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    mtext("Norm.Vector\nDensity", side = 2, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[2], main = "", xlab = "", ylab = "", las = 1,
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
      mtext(colnames(mean_operation_values)[j], side = 3, outer = FALSE, cex = 1, font = 2, line = 1, col = "grey20")
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  } else {
    subsetter <- cluster_membership == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[1], lwd = 2, ylab = "", xlab = "", las = 1)
    mtext("Time (Months)", side = 1, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    mtext("Norm.Vector\nDensity", side = 2, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[1], main = "", xlab = "", ylab = "", las = 1,
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  }
  print(i)
}
dev.off()

## Supplementary Figure for 4 Clusters
num_clust <- 4
clustering_results <- kmeans(PCA_output[, 1:4], num_clust, nstart = 20)
cluster_membership <- clustering_results$cluster
cluster_output <- data.frame(id = metadata$id, country = metadata$country, city = metadata$city, cluster = clustering_results$cluster, reordered_mean_realisation)

cluster_df <- data.frame(cluster = factor(cluster_membership), id = seq(1:dim(normalised_output)[1]), normalised_output)
cluster_df <- cluster_df %>%
  pivot_longer(cols= X1:X25, names_to = "timepoint", values_to = "vector_density") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint))) %>%
  mutate(timepoint2 = ifelse(timepoint <= 12, timepoint + 13, timepoint - 12))

rainfall_df <- data.frame(cluster = factor(cluster_membership), id = seq(1:dim(normalised_output)[1]), reordered_rainfall)
rainfall_df <- rainfall_df %>%
  pivot_longer(cols= X1:X24, names_to = "timepoint", values_to = "rainfall_density") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint))) %>%
  mutate(timepoint2 = ifelse(timepoint <= 11, timepoint + 12, timepoint - 11))

cluster_df <- cluster_df %>%
  left_join(rainfall_df, by = c("cluster", "id", "timepoint2"))

mean_df <- cluster_df %>%
  group_by(timepoint2, cluster) %>%
  summarise(mean = mean(vector_density))
mean_rainfall_df <- cluster_df %>%
  group_by(timepoint2, cluster) %>%
  summarise(mean = mean(rainfall_density))

four_cluster <- ggplot() +
  geom_line(data = cluster_df, aes(x = timepoint2, y = 100 * vector_density * 25/12, col = cluster, group = id), alpha = 0.2) +
  coord_cartesian(ylim = c(0, 20 * 25/12)) +
  scale_x_continuous(breaks = seq(1, 23, 2), labels = 1:12) +
  facet_wrap(~cluster, nrow = 2) +
  theme_bw() +
  geom_line(data = mean_df, aes(x = timepoint2, y = 100 * mean * 25/12, col = cluster), size = 2) +
  #geom_line(data = mean_rainfall_df, aes(x = timepoint2, y = 1000 * mean * 25/12, group = cluster), col = "black") +
  scale_color_manual(values = palette()[6:3]) +
  labs(y = "Normalised Monthly Vector Density", x = "Peak-Standardised Time (Months)") +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank())
if (fresh_run) {
  ggsave(four_cluster, file = here("figures", "Supp_Figure_Four_Cluster_Sensitivity.pdf"), width = 7, height = 5.8)
}

mean_timings <- data.frame(features_df, cluster = cluster_membership) %>%
  group_by(cluster) %>%
  summarise(mean_rainfall_timing = mean(rainfall_peak_timing) * 12/24,
            se_rainfall = mean_rainfall_timing/sqrt(n()),
            mean_vector_timing = mean(dens_peak_timing) * 12/25,
            se_vector = mean_vector_timing/sqrt(n()),
            mean_peak_diff = mean(peak_diff) * 12/25,
            se_peak = mean_peak_diff/sqrt(n()),
            n = n())
x <- data.frame(peak_diff = features_df$peak_diff, cluster = factor(cluster_membership))
aov <- aov(peak_diff ~ cluster, x)
TukeyHSD(aov, "cluster")
