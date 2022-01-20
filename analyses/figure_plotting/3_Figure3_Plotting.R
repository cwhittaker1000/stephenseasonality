#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools);
library(dismo); library(gbm); library(mltools); library(glmnet); library(caret); library(themis)
library(tidymodels); library(doParallel); library(vip); library(forcats); library(vip);
library(RColorBrewer); library(corrplot); library(DALEXtra)

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))
source(here("functions", "gp_fitting_functions.R"))

# Loading in extracted stephensi, processed but unsmoothed data and renaming variables
overall <- readRDS(file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city

features_df <- readRDS(file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
features <- features_df[, 7:dim(features_df)[2]]

cluster_output <- readRDS(file = here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_output$cluster

rf_ups_full <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
rf_no_ups_full <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData.rds"))
rf_ups_subset <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_SubsetData.rds"))
rf_no_ups_subset <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_NoUpsampling_SubsetData.rds"))

# Plotting Random Forest Results
a <- bind_rows(rf_ups_full$importance)
a$data <- "full_data"
b <- bind_rows(rf_no_ups_full$importance)
b$data <- "full_data"
c <- bind_rows(rf_ups_subset$importance)
c$data <- "subset_data"
d <- bind_rows(rf_no_ups_subset$importance)
d$data <- "subset_data"



%>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x[rev(order(x$mean)), ]

x <- bind_rows(iterations_ups$importance) %>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x[rev(order(x$mean)), ]


mean(rf_ups_full$test_roc_auc)
mean(rf_no_ups_full$test_roc_auc)
mean(rf_ups_subset$test_roc_auc)
mean(rf_no_ups_subset$test_roc_auc)



x <- bind_rows(iterations_ups$importance, iterations$importance) %>%
  pivot_wider(names_from = sampling, values_from = Importance) %>%
  group_by(Variable) %>%
  summarise(mean_ups = mean(upsampling), 
            mean_no = mean(no_upsampling))


x <- bind_rows(iterations_ups$importance, iterations$importance) 
ggplot(x, aes(x= Variable, fill = sampling, y = Importance)) +
  geom_boxplot()



mean(rf_ups_full$test_accuracy)
mean(rf_no_ups_full$test_accuracy)
mean(rf_ups_subset$test_accuracy)
mean(rf_no_ups_subset$test_accuracy)

mean(rf_ups_full$test_one_peak_accuracy)
mean(rf_no_ups_full$test_one_peak_accuracy)
mean(rf_ups_subset$test_one_peak_accuracy)
mean(rf_no_ups_subset$test_one_peak_accuracy)

mean(rf_ups_full$test_two_peak_accuracy)
mean(rf_no_ups_full$test_two_peak_accuracy)
mean(rf_ups_subset$test_two_peak_accuracy)
mean(rf_no_ups_subset$test_two_peak_accuracy)


# Extracting Mean Realisation for Each Time-Series
set.seed(10)
interpolating_points <- 2
mean_realisation <- matrix(nrow = dim(overall)[1], ncol = (12 * interpolating_points + 1))
prior <- "informative"
for (i in 1:length(overall$id)) {
  
  index <- overall$id[i]
  
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
  
  print(i)
}
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Extracting Rainfall Data (also sort out leap year stuff)
leap_years <- c(1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2016,	2020)
months_length <- c(15, 16, 14, 14, 15, 16, 15, 15, 15, 16, 15, 15, 
                   15, 16, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16)
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
norm_rainfall_storage <- t(apply(rainfall_storage, 1, normalise_total))

# Exploring and Assessing Urban Association With Unimodal Dynamics

## Rural/urban associations with cluster membership
tab <- table(cluster_membership, urban_rural)
tab <- tab[, 2:3]
chisq.test(tab)
chisq.test(tab, simulate.p.value = TRUE)

## Rural/urban association with number of peaks
tab <- table(features[, "peaks"], urban_rural)
chisq.test(tab[1:2, 2:3])
chisq.test(tab[1:2, 2:3], simulate.p.value = TRUE)

## Average Incidence Per 4 Months
table(features[, "peaks"], urban_rural)

mean(features$per_ind_4_months[urban_rural == "Urban"]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 1]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 2]) # almost all urban time series are unimodal

mean(features$per_ind_4_months[urban_rural == "Rural"]) # mixture of unimodal and bimodal time series for rural
mean(features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 1]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 2]) # almost all urban time series are unimodal

# Need to explore the rural 1 vs 2 peaks in more depth - what's up with that and why are we seeing such
# divergent dynamics across those 2 groups.

palette(c("#70AA4D", "#6F8D98"))
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 0.5, 0.5), oma = c(3, 1, 1, 5))
max <- 12
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
lty <- 3
timepoints <- seq(0, 12, length = dim(normalised_output)[2])

for (i in c("Urban", "Rural")) {
  
  if (i == "Urban") {
    cluster <- normalised_output[urban_rural == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[2], las = 1, xaxt = "n")
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[2], alpha.f = 0.2))
    }
    text(0, max-0.5, paste0("Urban Settings (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[urban_rural == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
  }
  else if (i == "Rural") {
    cluster <- normalised_output[urban_rural == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[1], las = 1, xaxt = "n")
    axis(4, at = seq(0, 12, 2), las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[1], alpha.f = 0.2))
      
    }
    text(0, max-0.5, paste0("Rural Setings (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[urban_rural == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
    
    one_cluster <- normalised_output[urban_rural == i & features[, "peaks"] == 1, ]
    number_time_series <- dim(one_cluster)[1]
    plot(timepoints, apply(one_cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[1], las = 1, xaxt = "n")
    axis(4, at = seq(0, 12, 2), las = 2)
    for (j in 1:length(one_cluster[, 1])) {
      lines(timepoints, one_cluster[j, ] * 100, col = adjustcolor(palette()[1], alpha.f = 0.2))
      
    }
    text(0, max-0.5, paste0("Cluster 2 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[urban_rural == i & features[, "peaks"] == 1, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
    
    two_cluster <- normalised_output[urban_rural == i & features[, "peaks"] == 2, ]
    number_time_series <- dim(two_cluster)[1]
    plot(timepoints, apply(two_cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[1], las = 1, xaxt = "n")
    axis(4, at = seq(0, 12, 2), las = 2)
    for (j in 1:length(two_cluster[, 1])) {
      lines(timepoints, two_cluster[j, ] * 100, col = adjustcolor(palette()[1], alpha.f = 0.2))
      
    }
    text(0, max-0.5, paste0("Cluster 2 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[urban_rural == i & features[, "peaks"] == 2, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
    
  }
  
  
  mtext("Normalised Catch (% of Annual Total)", side = 4, outer = TRUE, cex = 1, font = 2, line = 3, col = "grey20")
}


boxplot(features$per_ind_4_months ~ urban_rural * features[, "peaks"], ylim = c(0.4, 1))

two_rural <- which(urban_rural == "Rural" & features[, "peaks"] == 2)
x <- overall[two_rural, ]
dev.off()
plot(features$per_ind_4_months, per_rain_4_months[1, ])

per_rain_4_months <- t(apply(norm_rainfall_storage, 1, percent_incidence, number_of_months = 4, timepoints_per_month = 2))
  
per_rain_4_months <- per_rain_4_months[1, ]

mean(per_rain_4_months[urban_rural == "Rural" & features[, "peaks"] == 1])
mean(per_rain_4_months[urban_rural == "Rural" & features[, "peaks"] == 2])
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 1])
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 2])

mean(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 1])
mean(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 2])
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 1])
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 2])



plot(apply(norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 1, ], 2, mean), type = "l")
lines(apply(norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 2, ], 2, mean), type = "l", col = "red")

for (i in 1:dim(norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 2, ])[1]) {
  plot(norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 2, ][i, ], type = "l")
  browser()
}


par(mfrow = c(3, 5), mar = c(1, 2, 1, 1), oma = c(3, 3, 3, 3))
rainfall_subset <- norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 2, ]
mosquito_subset <- normalised_output[urban_rural == "Rural" & features[, "peaks"] == 2, ]
for (i in 1:dim(mosquito_subset)[1]) {
  
  scale_factor <- max(mosquito_subset[i, ] * 100)/max(rainfall_subset[i, ] * 100)
  
  plot(timepoints, mosquito_subset[i, ] * 100, type = "l", ylab = NA, xlab = NA,
       lwd = 5, col = palette()[2], las = 1, xaxt = "n", ylim = c(0, max(mosquito_subset[i, ] * 100)))
  lines(timepoints[-length(timepoints)], rainfall_subset[i, ] * 100 * scale_factor, type = "l", col = "black", 
        lwd = 1, lty = lty)
}

par(mfrow = c(3, 8), mar = c(1, 2, 1, 1), oma = c(3, 3, 3, 3))
rainfall_subset <- norm_rainfall_storage[urban_rural == "Rural" & features[, "peaks"] == 1, ]
mosquito_subset <- normalised_output[urban_rural == "Rural" & features[, "peaks"] == 1, ]
for (i in 1:dim(mosquito_subset)[1]) {
  
  scale_factor <- max(mosquito_subset[i, ] * 100)/max(rainfall_subset[i, ] * 100)
  
  plot(timepoints, mosquito_subset[i, ] * 100, type = "l", ylab = NA, xlab = NA,
       lwd = 5, col = palette()[2], las = 1, xaxt = "n", ylim = c(0, max(mosquito_subset[i, ] * 100)))
  lines(timepoints[-length(timepoints)], rainfall_subset[i, ] * 100 * scale_factor, type = "l", col = "black", 
        lwd = 1, lty = lty)
}

par(mfrow = c(3, 8), mar = c(1, 2, 1, 1), oma = c(3, 3, 3, 3))
rainfall_subset <- norm_rainfall_storage[urban_rural == "Urban" & features[, "peaks"] == 1, ]
mosquito_subset <- normalised_output[urban_rural == "Urban" & features[, "peaks"] == 1, ]
for (i in 1:dim(mosquito_subset)[1]) {
  
  scale_factor <- max(mosquito_subset[i, ] * 100)/max(rainfall_subset[i, ] * 100)
  
  plot(timepoints, mosquito_subset[i, ] * 100, type = "l", ylab = NA, xlab = NA,
       lwd = 5, col = palette()[2], las = 1, xaxt = "n", ylim = c(0, max(mosquito_subset[i, ] * 100)))
  lines(timepoints[-length(timepoints)], rainfall_subset[i, ] * 100 * scale_factor, type = "l", col = "black", 
        lwd = 1, lty = lty)
}

par(mfrow = c(3, 8), mar = c(1, 2, 1, 1), oma = c(3, 3, 3, 3))
rainfall_subset <- norm_rainfall_storage[urban_rural == "Urban" & features[, "peaks"] == 2, ]
mosquito_subset <- normalised_output[urban_rural == "Urban" & features[, "peaks"] == 2, ]
for (i in 1:dim(mosquito_subset)[1]) {
  
  scale_factor <- max(mosquito_subset[i, ] * 100)/max(rainfall_subset[i, ] * 100)
  
  plot(timepoints, mosquito_subset[i, ] * 100, type = "l", ylab = NA, xlab = NA,
       lwd = 5, col = palette()[2], las = 1, xaxt = "n", ylim = c(0, max(mosquito_subset[i, ] * 100)))
  lines(timepoints[-length(timepoints)], rainfall_subset[i, ] * 100 * scale_factor, type = "l", col = "black", 
        lwd = 1, lty = lty)
}

plot(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 1],
     features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 1], pch = 20, xlim = c(0, 1), ylim = c(0, 1))
points(per_rain_4_months[urban_rural == "Rural" & features[, "peaks"] == 2],
       features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 2], pch = 20, col = "blue")
points(per_rain_4_months[urban_rural == "Rural" & features[, "peaks"] == 1],
       features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 1], pch = 20, col = "green")

mean(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 2])

mean(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 1])
mean(per_rain_4_months[urban_rural == "Urban" & features[, "peaks"] == 2])


plot(norm_rainfall_storage[1, ])
percent_incidence(fitting_output = norm_rainfall_storage[1, ], number_of_months = 4, timepoints_per_month = 2)


