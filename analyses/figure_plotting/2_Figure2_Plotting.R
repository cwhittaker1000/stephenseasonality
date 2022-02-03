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

# Loading in extracted stephensi, processed but unsmoothed data and renaming variables
overall <- readRDS(file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
features_df <- readRDS(file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
cluster_output <- readRDS(file = here("data", "systematic_review_results", "cluster_membership.rds"))
num_clust <- max(cluster_output$cluster)

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

# Normalising features & fitted output
features <- features_df[, 7:dim(features_df)[2]]
normalised_features <- scale(features)
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Running PCA on normalised features
set.seed(10)
PCA <- prcomp(normalised_features)
summary <- summary(PCA)
loadings <- PCA$rotation[, 1:7]
PCA_output <- as.matrix(normalised_features) %*% loadings

# Visualising the time-series belonging to each cluster
colours <- palette()[1:num_clust]
timepoints <- seq(0, 12, length = dim(normalised_output)[2])
clusters <- cluster_output$cluster
table(clusters)

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
norm_rainfall_storage <- normalise_total(rainfall_storage)

# Plotting Figure 2 and the Clusters
pdf("Figures/Figure_2_Overall.pdf", height = 6.5, width = 10.5, useDingbats = FALSE)
palette(c("#E0521A", "#3F88C5", "#44BBA4", "#393E41"))
max <- 12
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
lty <- 3
layout.matrix <- matrix(c(1, 2, 5, 3, 4, 6), nrow = 2, ncol = 3, byrow = TRUE)
layout(mat = layout.matrix)
par(mar = c(0.5, 0.5, 0.5, 0.5), oma = c(3, 5, 1, 1))
max <- 12
lty <- 2
for (i in 1:num_clust) {
  if (i == 1) {
    cluster <- normalised_output[clusters == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    # cm <- which(clusters == i)
    # cross_cor <- mean(mapply(calc_cross_cor, rainfall = as.data.frame(t(norm_rainfall_storage[cm, ])), index = 1:length(cm), MoreArgs = list(mosquitoes = normalised_output[cm, -25])))
    # text(0, max-0.5, paste0("Cluster 1 (n = ", number_time_series, ")", "\nr = ", round(cross_cor, 2)), cex = 1.2, col = "grey20", adj = 0, font = 2)
    text(0, max-0.5, paste0("Cluster 1 (n = ", number_time_series, ")"), cex = 1.7, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1.5, lty = lty)
    axis(2, at = seq(0, 12, 2), las = 2)
  }
  else if (i == 2) {
    cluster <- normalised_output[clusters == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
      
    }
    # cm <- which(clusters == i)
    # cross_cor <- mean(mapply(calc_cross_cor, rainfall = as.data.frame(t(norm_rainfall_storage[cm, ])), index = 1:length(cm), MoreArgs = list(mosquitoes = normalised_output[cm, -25])))
    # text(0, max-0.5, paste0("Cluster 2 (n = ", number_time_series, ")", "\nr = ", round(cross_cor, 2)), cex = 1.2, col = "grey20", adj = 0, font = 2)
    text(0, max-0.5, paste0("Cluster 2 (n = ", number_time_series, ")"), cex = 1.7, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1.5, lty = lty)
  } else if (i == 3) {
    cluster <- normalised_output[clusters == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    axis(1, at = seq(0, 11, 1), labels = months, las = 2)
    axis(2, at = seq(0, 12, 2), las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    # cm <- which(clusters == i)
    # cross_cor <- mean(mapply(calc_cross_cor, rainfall = as.data.frame(t(norm_rainfall_storage[cm, ])), index = 1:length(cm), MoreArgs = list(mosquitoes = normalised_output[cm, -25])))
    # text(0, max-0.5, paste0("Cluster 3 (n = ", number_time_series, ")", "\nr = ", round(cross_cor, 2)), cex = 1.2, col = "grey20", adj = 0, font = 2)
    text(0, max-0.5, paste0("Cluster 3 (n = ", number_time_series, ")"), cex = 1.7, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1.5, lty = lty)
  } else {
    cluster <- normalised_output[clusters == i, ]
    number_time_series <- dim(cluster)[1]
    plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
    axis(1, at = seq(0, 11, 1), labels = months, las = 2)
    for (j in 1:length(cluster[, 1])) {
      lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
    }
    # cm <- which(clusters == i)
    # cross_cor <- mean(mapply(calc_cross_cor, rainfall = as.data.frame(t(norm_rainfall_storage[cm, ])), index = 1:length(cm), MoreArgs = list(mosquitoes = normalised_output[cm, -25])))
    # text(0, max-0.5, paste0("Cluster 4 (n = ", number_time_series, ")", "\nr = ", round(cross_cor, 2)), cex = 1.2, col = "grey20", adj = 0, font = 2)
    text(0, max-0.5, paste0("Cluster 4 (n = ", number_time_series, ")"), cex = 1.7, col = "grey20", adj = 0, font = 2)
    mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
    lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1.5, lty = lty)
  }
  mtext("Normalised Catch (% of Annual Total)", side = 2, outer = TRUE, cex = 1, font = 2, line = 2, col = "grey20")
}

par(mar = c(0.5, 4, 0.5, 5))
set.seed(11)
boxplot(features_df$per_ind_4_months ~ clusters, border = palette()[1:4], col = NA, outline = FALSE,
        xlab = "", lty = 1, horizontal = FALSE, las = 1, yaxt = "n", ylab = NA, xaxt = "n")
axis(4, at = seq(0, 1, 0.1), las = 2)
mtext("% of Total Annual Catch \nin 4 Months", side=4, line=4, cex.lab=1, las=0, col="black")
stripchart(features_df$per_ind_4_months ~ clusters, method = "jitter", jitter = 0.25,
           pch = 20, cex = 1.5, col = palette()[1:4], vertical = TRUE, add = TRUE)

percent_incidence_out <- apply(norm_rainfall_storage, 1, percent_incidence, 4, 2)
boxplot(percent_incidence_out ~ clusters, border = palette()[1:4], col = NA, outline = FALSE,
        xlab = "", lty = 1, horizontal = FALSE, las = 1, yaxt = "n", ylab = NA)
axis(4, at = seq(0, 1, 0.1), las = 2)
mtext("% of Total Annual Rainfall \nin 4 Months", side=4, line=4, cex.lab=1, las=0, col="black")
stripchart(percent_incidence_out ~ clusters, method = "jitter", jitter = 0.25,
           pch = 20, cex = 1.5, col = palette()[1:4], vertical = TRUE, add = TRUE)
dev.off()

# Visualising the properties of each cluster
pdf("Figures/Supp_Figure_Cluster_Properties.pdf", height = 8, width = 10, useDingbats = FALSE)
number_properties <- ncol(normalised_features)
property_names <- colnames(normalised_features)
par(mfrow = c(num_clust, number_properties + 1), mar = c(3.5, 2, 2, 1), oma = c(1, 2.5, 5, 1))
max <- 12
mean_operation_values <- matrix(nrow = num_clust, ncol = number_properties)
colnames(mean_operation_values) <- c("Period", "Prop. Points\n1.6x Mean", "Dist\nfrom Jan",
                                     "Number of\nPeaks", "Von Mises\nMean", "Von Mises\nWeight", "% Incidence\nIn 4 Months")
for (i in 1:num_clust) {
  if (i == 1) {
    subsetter <- clusters == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2, ylab = "", xlab = "", las = 1)
    mtext("Time (Months)", side = 1, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    mtext("Norm.Vector\nDensity", side = 2, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlab = "", ylab = "", las = 1,
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
      mtext(colnames(mean_operation_values)[j], side = 3, outer = FALSE, cex = 1, font = 2, line = 1, col = "grey20")
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  } else {
    subsetter <- clusters == i
    cluster_time_series <- normalised_output[subsetter, ]
    cluster_time_series_properties <- normalised_features[subsetter, ]
    plot(timepoints, apply(cluster_time_series, 2, mean) * 100, type = "l", ylim = c(0, max), col = palette()[i], lwd = 2, ylab = "", xlab = "", las = 1)
    mtext("Time (Months)", side = 1, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    mtext("Norm.Vector\nDensity", side = 2, outer = FALSE, cex = 0.75, font = 2, line = 2, col = "grey20")
    for (j in 1:number_properties) {
      hist(cluster_time_series_properties[, j], col = palette()[i], main = "", xlab = "", ylab = "", las = 1,
           xlim = c(min(cluster_time_series_properties[, j]), max(cluster_time_series_properties[, j])))
    }
    mean_operation_values[i, ] <- apply(cluster_time_series_properties, 2, mean)
  }
}

dev.off()

mean(features$peaks[urban_rural == "Urban"]) # almost all urban time series are unimodal
mean(features$peaks[urban_rural == "Rural"]) # mixture of unimodal and bimodal time series for rural

mean(features$per_ind_4_months[cluster_output$cluster == 1])
mean(features$per_ind_4_months[cluster_output$cluster == 2])
mean(features$per_ind_4_months[cluster_output$cluster == 3])
mean(features$per_ind_4_months[cluster_output$cluster == 4])

mean(features$jan_dist[cluster_output$cluster == 1])
mean(features$jan_dist[cluster_output$cluster == 2])
mean(features$jan_dist[cluster_output$cluster == 3])
mean(features$jan_dist[cluster_output$cluster == 4])

# Note: mapply works with either vectors or lists - dataframes are lists so this works (do need the t())
# calc_cross_cor <- function(rainfall, index, mosquitoes) {
#   return(cor(rainfall, mosquitoes[index, ]))
# }
# 
# lags <- vector(mode = "numeric", length = 65L)
# max <-  vector(mode = "numeric", length = 65L)
# ccf_1 <-  vector(mode = "numeric", length = 65L)
# for (i in 1:65) {
#   x <- stats::ccf(x = normalised_output[i, -25],
#                   y = norm_rainfall_storage[i, ], 
#                   type = "correlation", plot = FALSE)
#   max_ccf <- max(x$acf)
#   max_lag <- x$lag[x$acf == max_ccf]
#   lags[i] <- max_lag
#   max[i] <- max_ccf
#   ccf_1[i] <- x$acf[x$lag == -2]
# }
# 
# mean(lags[clusters == 1])
# mean(lags[clusters == 2])
# mean(lags[clusters == 3])
# mean(lags[clusters == 4])
# 
# mean(max[clusters == 1])
# mean(max[clusters == 2])
# mean(max[clusters == 3])
# mean(max[clusters == 4])
# 
# mean(ccf_1[clusters == 1])
# mean(ccf_1[clusters == 2])
# mean(ccf_1[clusters == 3])
# mean(ccf_1[clusters == 4])

# for (i in 1:num_clust) {
#   if (i == 1) {
#     cluster <- normalised_output[clusters == i, ]
#     number_time_series <- dim(cluster)[1]
#     plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
#     for (j in 1:length(cluster[, 1])) {
#       lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
#     }
#     text(0, max-0.5, paste0("Cluster 1 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
#     mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
#     lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
#   }
#   else if (i == 2) {
#     cluster <- normalised_output[clusters == i, ]
#     number_time_series <- dim(cluster)[1]
#     plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
#     axis(4, at = seq(0, 12, 2), las = 2)
#     for (j in 1:length(cluster[, 1])) {
#       lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
#       
#     }
#     text(0, max-0.5, paste0("Cluster 2 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
#     mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
#     lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
#   } else if (i == 3) {
#     cluster <- normalised_output[clusters == i, ]
#     number_time_series <- dim(cluster)[1]
#     plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
#     axis(1, at = seq(0, 11, 1), labels = months, las = 2)
#     for (j in 1:length(cluster[, 1])) {
#       lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
#     }
#     text(0, max-0.5, paste0("Cluster 3 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
#     mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
#     lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
#   } else {
#     cluster <- normalised_output[clusters == i, ]
#     number_time_series <- dim(cluster)[1]
#     plot(timepoints, apply(cluster, 2, mean) * 100, type = "l", yaxt = "n", ylim = c(0, max), lwd = 5, col = palette()[i], las = 1, xaxt = "n")
#     axis(4, at = seq(0, 12, 2), las = 2)
#     axis(1, at = seq(0, 11, 1), labels = months, las = 2)
#     for (j in 1:length(cluster[, 1])) {
#       lines(timepoints, cluster[j, ] * 100, col = adjustcolor(palette()[i], alpha.f = 0.2))
#     }
#     text(0, max-0.5, paste0("Cluster 4 (n = ", number_time_series, ")"), cex = 1.2, col = "grey20", adj = 0, font = 2)
#     mean_rainfall <- apply(norm_rainfall_storage[clusters == i, ], 2, mean)
#     lines(timepoints[-length(timepoints)], mean_rainfall * 100, type = "l", col = "black", lwd = 1, lty = lty)
#   }
#   mtext("Normalised Catch (% of Annual Total)", side = 4, outer = TRUE, cex = 1, font = 2, line = 3, col = "grey20")
# }