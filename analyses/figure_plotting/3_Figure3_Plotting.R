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
library(RColorBrewer); library(corrplot); library(DALEXtra); library(cowplot)

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))
source(here("functions", "gp_fitting_functions.R"))

# Loading in extracted stephensi, processed but unsmoothed data and renaming variables
overall <- readRDS(file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city

features_df <- readRDS(file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
features <- features_df[, 7:dim(features_df)[2]]
envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE)))
envt_variables <- features_df %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))


cluster_output <- readRDS(file = here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_output$cluster

rf_ups_full <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
rf_no_ups_full <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData.rds"))
rf_ups_subset <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_SubsetData.rds"))
rf_no_ups_subset <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_NoUpsampling_SubsetData.rds"))

# Plotting AUC Results
for (i in 1:length(rf_no_ups_full$test_roc_curve)) {
  temp <- rf_no_ups_full$test_roc_curve[[i]]
  temp_ups <- rf_ups_full$test_roc_curve[[i]]
  if (i == 1) {
    df <- temp
    df$iteration <- 1
    df_ups <- temp_ups
    df_ups$iteration <- 1
  } else {
    temp$iteration <- i
    df <- rbind(df, temp)
    temp_ups$iteration <- i
    df_ups <- rbind(df_ups, temp_ups)
  }
  print(i)
}

a_plot <- ggplot(df, aes(x = 1-specificity, y = sensitivity, col = factor(iteration))) +
  geom_path() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = "black", lty = 2) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  theme(legend.position = "none") +
  annotate("label", x = 0.20, y = 0.97,
           label = paste0("Mean AUC = ", round(mean(rf_no_ups_full$test_roc_auc), 2)),
           label.padding = unit(0.35, "lines"), label.r = unit(0, "lines"),
           label.size = unit(0.35, "lines"), size = 5)
b <- bind_rows(rf_no_ups_full$importance)
b$data <- "full_data"
imp <- b %>%
  group_by(Variable) %>%
  summarise(mean_Importance = mean(Importance),
            stdev_Importance = sd(Importance),
            stder_Importance = sd(Importance)/sqrt(n()))
imp$lower <- pmax(rep(0, length(imp$mean_Importance)), imp$mean_Importance - 1.96 * imp$stdev_Importance)
var_names <- imp$Variable[order(imp$mean_Importance)]
new_names <- c("Study\nfrom\nIndia", "LC 180", "LC 150", "LC 11", 
               "Study\nfrom\nIran", "Temperature\nSeasonality", "LC 130", "LC 110",
               "LC 122", "LC 120", "Rainfall\nColdest\nQuarter", "Rainfall\nSeasonality",
               "LC 20", "LC 30", "Population\nPer\nSquare Km", "LC 10")
b_plot <- ggplot(imp, aes(x = reorder(Variable, mean_Importance), y = mean_Importance, 
                fill = mean_Importance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = pmax(0, mean_Importance - 1.96 * stdev_Importance),
                    ymax = mean_Importance + 1.96 * stdev_Importance)) +
  scale_x_discrete(labels = new_names) +
  xlab("") + ylab("Variable Importance") +
  lims(y = c(0, 0.036)) +
  theme(legend.position = "none")

cowplot::plot_grid(a_plot, b_plot, rel_widths = c(1, 2), align = "h", axis = "b")

a_plot_ups <- ggplot(df_ups, aes(x = 1-specificity, y = sensitivity, col = factor(iteration))) +
  geom_path() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = "black", lty = 2) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  theme(legend.position = "none") +
  annotate("label", x = 0.20, y = 0.97,
           label = paste0("Mean AUC = ", round(mean(rf_ups_full$test_roc_auc), 2)),
           label.padding = unit(0.35, "lines"), label.r = unit(0, "lines"),
           label.size = unit(0.35, "lines"), size = 5)
b_ups <- bind_rows(rf_ups_full$importance)
b_ups$data <- "full_data"
imp_ups <- b_ups %>%
  group_by(Variable) %>%
  summarise(mean_Importance = mean(Importance),
            stdev_Importance = sd(Importance),
            stder_Importance = sd(Importance)/sqrt(n()))
imp_ups$lower <- pmax(rep(0, length(imp_ups$mean_Importance)), imp_ups$mean_Importance - 1.96 * imp_ups$stdev_Importance)
var_names_ups <- imp_ups$Variable[order(imp_ups$mean_Importance)]
new_names_ups <- c("Study\nfrom\nIndia", "LC 180", "LC 150", "LC 130", "LC 11", "Rainfall\nColdest\nQuarter", 
                   "LC 120", "LC 20", "LC 122", "LC 10", "Rainfall\nSeasonality", "LC 110",
                   "Temperature\nSeasonality", "LC 30", "Study\nfrom\nIran", "Population\nPer\nSquare Km")
b_plot_ups <- ggplot(imp_ups, aes(x = reorder(Variable, mean_Importance), y = mean_Importance, 
                     fill = mean_Importance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = pmax(0, mean_Importance - 1.96 * stdev_Importance),
                    ymax = mean_Importance + 1.96 * stdev_Importance)) +
  scale_x_discrete(labels = new_names_ups) +
  xlab("") + ylab("Variable Importance") +
  lims(y = c(0, 0.064)) +
  theme(legend.position = "none")


cowplot::plot_grid(a_plot, b_plot, 
                   a_plot_ups, b_plot_ups, nrow = 2, ncol =2,
                   rel_widths = c(1, 2), align = "h", axis = "b")



# Plotting Random Forest Results
a <- bind_rows(rf_ups_full$importance)
a$data <- "full_data"
b <- bind_rows(rf_no_ups_full$importance)
b$data <- "full_data"
c <- bind_rows(rf_ups_subset$importance)
c$data <- "subset_data"
d <- bind_rows(rf_no_ups_subset$importance)
d$data <- "subset_data"


mean(envt_variables$population_per_1km[urban_rural == "Urban"])
mean(envt_variables$population_per_1km[urban_rural == "Rural"])

mean(envt_variables$LC_10[urban_rural == "Urban"])
mean(envt_variables$LC_10[urban_rural == "Rural"])

mean(envt_variables$LC_210[urban_rural == "Urban"])
mean(envt_variables$LC_210[urban_rural == "Rural"])



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


