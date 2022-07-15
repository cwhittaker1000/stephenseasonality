# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels); library(cowplot)
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
dim <- 5
not_sampled_vector_peak_absolute_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_rainfall_peak_absolute_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_rainfall_peak_plus_absolute_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_annual_avg_absolute_summary <- array(data = NA, dim = c(65, dim, dim))

not_sampled_vector_annual_diff_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_vector_rainfall_diff_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_vector_rainfall_peak_plus_diff_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_rainfall_annual_diff_summary <- array(data = NA, dim = c(65, dim, dim))
not_sampled_rainfall_peak_plus_annual_diff_summary <- array(data = NA, dim = c(65, dim, dim))

for (i in 1:65) {
  ts_peak_vector_month <- peak_density_month[i] 
  ts_peak_rainfall_month <- peak_rainfall_month[i] 
  ts_peak_rainfall_plus_month <- (ts_peak_rainfall_month + 1) %% 12
  ts_peak_rainfall_plus_month <- ifelse(ts_peak_rainfall_plus_month == 0, 12, ts_peak_rainfall_plus_month)
  
  not_sampled_vector_peak_absolute_summary[i, 1:dim, 1:dim] <- 1- prob_not_sampled_poisson[i, 1:dim, 1:dim, ts_peak_vector_month]
  not_sampled_rainfall_peak_absolute_summary[i, 1:dim, 1:dim] <- 1 - prob_not_sampled_poisson[i, 1:dim, 1:dim, ts_peak_rainfall_month]
  not_sampled_rainfall_peak_plus_absolute_summary[i, 1:dim, 1:dim] <- 1 - prob_not_sampled_poisson[i, 1:dim, 1:dim, ts_peak_rainfall_plus_month]
  not_sampled_annual_avg_absolute_summary[i, 1:dim, 1:dim] <- 1 - apply(prob_not_sampled_poisson[i, 1:dim, 1:dim, 1:12], c(1, 2), mean)
  
  not_sampled_vector_annual_diff_summary[i, 1:dim, 1:dim] <- not_sampled_vector_peak_absolute_summary[i, 1:dim, 1:dim] - not_sampled_annual_avg_absolute_summary[i, 1:dim, 1:dim]
  not_sampled_vector_rainfall_diff_summary[i, 1:dim, 1:dim] <- not_sampled_vector_peak_absolute_summary[i, 1:dim, 1:dim] - not_sampled_rainfall_peak_absolute_summary[i, 1:dim, 1:dim]
  not_sampled_vector_rainfall_peak_plus_diff_summary[i, 1:dim, 1:dim] <- not_sampled_vector_peak_absolute_summary[i, 1:dim, 1:dim] - not_sampled_rainfall_peak_plus_absolute_summary[i, 1:dim, 1:dim]
  not_sampled_rainfall_annual_diff_summary[i, 1:dim, 1:dim] <- not_sampled_rainfall_peak_absolute_summary[i, 1:dim, 1:dim] - not_sampled_annual_avg_absolute_summary[i, 1:dim, 1:dim]
  not_sampled_rainfall_peak_plus_annual_diff_summary[i, 1:dim, 1:dim] <- not_sampled_rainfall_peak_plus_absolute_summary[i, 1:dim, 1:dim] - not_sampled_annual_avg_absolute_summary[i, 1:dim, 1:dim]
}

size <- 12
temp <- apply(not_sampled_annual_avg_absolute_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
a <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "magma", name = "p(detect)",
                       limits = c(0, 1)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Probability of detecting Anopheles\nstephensi when starting sampling\nat random month") 

temp <- apply(not_sampled_rainfall_peak_absolute_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
d <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(temp_long$X)),
    position = "left",
    labels = 1:dim) +
  scale_fill_viridis_c(option = "magma", name = "p(detect)",
                       limits = c(0, 1)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Probability of detecting Anopheles\nstephensi when starting sampling\nat peak rainfall month") 

temp <- apply(not_sampled_vector_peak_absolute_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
b <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(#limits = rev(unique(temp_long$X)),
    position = "left",
    labels = 1:dim) +
  scale_fill_viridis_c(option = "magma", name = "p(detect)",
                       limits = c(0, 1)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Probability of detecting Anopheles\nstephensi when starting sampling\nat peak vector density month") 

temp <- apply(not_sampled_vector_annual_diff_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
max_val <- max(temp_long$value)
c <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "viridis", name = "Increase in\np(detect)",
                       limits = c(0, max_val)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Increase in probability of detecting Anopheles\nstephensi by starting sampling at vector peak\nrather than a random month") 

temp <- apply(not_sampled_rainfall_annual_diff_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
e <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "viridis", name = "Increase in\np(detect)",
                       limits = c(0, max_val)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Increase in probability of detecting Anopheles\nstephensi by starting sampling at rainfall peak\nrather than a random month") 

temp <- apply(not_sampled_rainfall_annual_diff_summary, c(2, 3), mean)
temp <- data.frame(temp)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
e <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "viridis", name = "Increase in\np(detect)",
                       limits = c(0, max_val)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Increase in probability of detecting Anopheles\nstephensi by starting sampling at rainfall peak\nrather than a random month") 

fold_change_vec_ann <- apply(not_sampled_vector_peak_absolute_summary/not_sampled_annual_avg_absolute_summary, c(2, 3), mean)
temp <- data.frame(fold_change_vec_ann)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
max_val2 <- max(temp_long$value)
f <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "mako", name = "Fold Change in\np(detect)",
                       limits = c(1, max_val2)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Fold change in probability of detecting Anopheles\nstephensi by starting sampling at vector peak\nrather than a random month") 

fold_change_rain_ann <- apply(not_sampled_rainfall_peak_absolute_summary/not_sampled_annual_avg_absolute_summary, c(2, 3), mean)
temp <- data.frame(fold_change_rain_ann)
temp$Y <- paste0("Y", 1:dim)
temp$Y <- factor(temp$Y, levels = paste0("Y", 1:dim))
temp_long <- temp %>%
  pivot_longer(cols = -Y, names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = colnames(temp)[-length(colnames(temp))])
g <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(limits = rev(unique(temp_long$X)),
                   position = "left",
                   labels = 1:dim) +
  scale_fill_viridis_c(option = "mako", name = "Fold Change in\np(detect)",
                       limits = c(1, max_val2)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) + #month.abb
  labs(y = "Sampling Days Per Month",
       x = "Number of Months Sampled",
       title = "Fold change in probability of detecting Anopheles\nstephensi by starting sampling at rainfall peak\nrather than a random month") 

legend <- get_legend(a + theme(legend.box.margin = margin(0, 0, 0, 2)))
temp_plot <- plot_grid(a + theme(legend.position="none"), 
                       d + theme(legend.position = "none"),
                       b + theme(legend.position="none"), 
                       nrow = 1) 
temp_plot2 <- plot_grid(temp_plot, legend, rel_widths = c(3, 0.4))
temp_plot3 <- plot_grid(NA,
                        e + theme(legend.position="none"), 
                        c + theme(legend.position="none"),
                        ncol = 3)
legend2 <- get_legend(c + theme(legend.box.margin = margin(0, 0, 0, 2)))
temp_plot4 <- plot_grid(temp_plot3, legend2, rel_widths = c(3, 0.4))
temp_plot5 <- plot_grid(NA,
                        g + theme(legend.position="none"), 
                        f + theme(legend.position="none"),
                        ncol = 3)
legend3 <- get_legend(f + theme(legend.box.margin = margin(0, 0, 0, 2)))
temp_plot6 <- plot_grid(temp_plot5, legend3, rel_widths = c(3, 0.4))
plot_grid(temp_plot2, temp_plot4, temp_plot6, nrow = 3) 

# Single day per month sampling, seeing the effect of 
# 1st dim is the time series
# 2nd dim is the number of months sampled
# 3rd dim is the number of nights per month sampled
# 4th dim is which month you start at 

# rearrange so that time-series peaks are all aligned
prob_not_sampled_poisson_rearranged <- array(data = NA, dim = c(65, 12, 10, 12))
for (i in 1:65) {
  index <- peak_density_month[i]
  if (index == 1) {
    prob_not_sampled_poisson_rearranged[i, , , ] <- prob_not_sampled_poisson[i, , , ]
  } else {
    prob_not_sampled_poisson_rearranged[i, , , ] <- prob_not_sampled_poisson[i, , , c(index:12, 1:(index - 1))]
  }
}

x <- 1 - prob_not_sampled_poisson_rearranged[, 3, 1, ] # 3 months sampled, 1 night per month
x <- data.frame(id = 1:65, cluster = overall$cluster, x) %>%
  pivot_longer(-c(id, cluster), names_to = "month_start", values_to = "prob_detect")
x$month_start <- factor(x$month_start, levels = paste0("X", 1:12))
x$id <- factor(x$id)
x$cluster <- factor(x$cluster)
y <- x %>%
  group_by(cluster, month_start) %>%
  summarise(mean = mean(prob_detect))

vec <- ggplot(x, aes(x = month_start, y = prob_detect, group = id, col = cluster)) +
  geom_line(alpha = 0.2) +
  geom_line(data = y, aes(x = month_start, y = mean, group = cluster), size = 1.5) +
  scale_x_discrete(labels =  0:11) +
  scale_colour_manual(values = c("#DF536B", "black")) +
  labs(x = "Month Sampling Started Relative to Peak Vector Month", y = "Probability of Detection",
       title = "Sampling Based on Vector Timing")

prob_not_sampled_poisson_rearranged_rainfall <- array(data = NA, dim = c(65, 12, 10, 12))
for (i in 1:65) {
  index <- peak_rainfall_month[i]
  if (index == 1) {
    prob_not_sampled_poisson_rearranged_rainfall[i, , , ] <- prob_not_sampled_poisson[i, , , ]
  } else {
    prob_not_sampled_poisson_rearranged_rainfall[i, , , ] <- prob_not_sampled_poisson[i, , , c(index:12, 1:(index - 1))]
  }
}

x <- 1 - prob_not_sampled_poisson_rearranged_rainfall[, 3, 1, ] # 3 months sampled, 1 night per month
x <- data.frame(id = 1:65, cluster = overall$cluster, x) %>%
  pivot_longer(-c(id, cluster), names_to = "month_start", values_to = "prob_detect")
x$month_start <- factor(x$month_start, levels = paste0("X", 1:12))
x$id <- factor(x$id)
x$cluster <- factor(x$cluster)
y <- x %>%
  group_by(cluster, month_start) %>%
  summarise(mean = mean(prob_detect))

rain <- ggplot(x, aes(x = month_start, y = prob_detect, group = id, col = cluster)) +
  geom_line(alpha = 0.2) +
  geom_line(data = y, aes(x = month_start, y = mean, group = cluster), size = 1.5) +
  scale_x_discrete(labels =  0:11) +
  scale_colour_manual(values = c("#DF536B", "black")) +
  labs(x = "Month Sampling Started Relative to Peak Rainfall Month", y = "Probability of Detection",
       title = "Sampling Based on Rainfall Timing")


plot_grid(rain + theme(legend.position = "none"), 
          vec + theme(legend.position = "none"), nrow = 1)



# x <- 1 - prob_not_sampled_poisson[, 2, 1, ] # 3 months sampled, 1 night per month
# x <- data.frame(id = 1:65, cluster = overall$cluster, x) %>%
#   pivot_longer(-c(id, cluster), names_to = "month_start", values_to = "prob_detect")
# x$month_start <- factor(x$month_start, levels = paste0("X", 1:12))
# x$id <- factor(x$id)
# x$cluster <- factor(x$cluster)
# y <- x %>%
#   group_by(cluster, month_start) %>%
#   summarise(mean = mean(prob_detect))
# 
# ggplot(x, aes(x = month_start, y = prob_detect, group = id, col = cluster)) +
#   geom_line(alpha = 0.2) +
#   geom_line(data = y, aes(x = month_start, y = mean, group = cluster)) +
#   scale_x_discrete(labels =  month.abb) +
#   labs(x = "Month Sampling Started", y = "Probability of Detection")


