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
metadata <- raw_df[retain_index, ] %>%
  dplyr::select(Time.Series.ID, Country, Admin.1, Admin.2, City., Year.Start, Year.End) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2, city = City., 
         start = Year.Start, end = Year.End)
overall <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))

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

# Extracting Air Temp and Rainfall Data (also sort out leap year stuff)
leap_years <- c(1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008, 2016,	2020)
months_length <- c(15, 16, 14, 14, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16, 15, 16, 15, 15, 15, 16, 15, 15, 15, 16)
rainfall_storage <- matrix(nrow = dim(overall)[1], ncol = length(months_length)) # this is 24 long rather than 25 like ento fits - why is this?? Need to check
airtemp_storage <- matrix(nrow = dim(overall)[1], ncol = length(months_length)) # this is 24 long rather than 25 like ento fits - why is this?? Need to check
for (i in 1:length(overall$id)) {
  index <- overall$id[i]
  rf_temp <- c()
  airtemp_temp <- c()
  
  rf <- read.csv(here(paste0("data/location_specific_rainfall/rainfall_ts", index, ".csv")))
  rf <- rf %>%
    group_by(daymonth_id) %>%
    summarise(rainfall = mean(rainfall))
  
  airtemp <- read.csv(here(paste0("data/location_specific_airtemp/airtemp_ts", index, ".csv")))
  airtemp <- airtemp %>%
    group_by(daymonth_id) %>%
    summarise(airtemp = mean(mean_2m_air_temperature_celsius))
  
  counter <- 1
  count_vec <- counter
  for (j in 1:length(months_length)) {
    indices <- counter:(counter + months_length[j] - 1)
    rf_temp <- c(rf_temp, sum(rf$rainfall[indices]))
    airtemp_temp <- c(airtemp_temp, mean(airtemp$airtemp[indices]))
    counter <- counter + months_length[j]
  }
  rainfall_storage[i, ] <- rf_temp
  airtemp_storage[i, ] <- airtemp_temp
  print(i)
}

airtemp_df <- data.frame(id = overall$id, airtemp_storage) %>%
  pivot_longer(-id, values_to = "airtemp", names_to = "timepoint") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint)))
rainfall_df <- data.frame(id = overall$id, rainfall_storage) %>%
  pivot_longer(-id, values_to = "rainfall", names_to = "timepoint") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint)))

overall_pv2 <- overall_pv %>%
  left_join(airtemp_df, by = c("id", "timepoint")) %>%
  left_join(rainfall_df, by = c("id", "timepoint"))
norm_vector <- t(apply(mean_realisation, 1, normalise_total)) 
norm_dens_df <- data.frame(id = overall$id, norm_vector) %>%
  pivot_longer(-id, values_to = "norm_vector_dens", names_to = "timepoint") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint)))
overall_pv3 <- overall_pv2 %>%
  left_join(norm_dens_df, by = c("id", "timepoint")) %>%
  mutate(bin_airtemp = cut(airtemp, breaks = seq(0, 45, 5)),
         bin_airtemp2 = as.numeric(bin_airtemp),
         airtemp_thresh = ifelse(airtemp >= 35, "too_high", 
                                 ifelse(airtemp <= 15, "too_low", "fine"))) %>%
  mutate(lagged_timepoint = ifelse(timepoint + 2 <= 25, timepoint + 2, (timepoint + 2) - 25))

lag <- 2
overall_pv3_subset <- overall_pv3 %>%
  filter(id %in% c(6, 19, 22, 36, 37, 52, 55, 63, 64, 65, 73)) 

coef <- 0.0043
int <- 12
ggplot(data = overall_pv3_subset) +
  geom_line(aes(x = timepoint, y = norm_vector_dens, col = factor(country))) +
  facet_wrap(~id, scales = "free_y") +
  scale_y_continuous(limits=c(0, 0.11), position = "left",
                     sec.axis = sec_axis(trans=~int + ./coef, name = "Air Temp")) +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(2, 24, length.out = 12)) +
  geom_line(aes(x = lagged_timepoint, y = coef * (airtemp - int))) +
  scale_fill_manual(values = gg_color_hue(6)) +
  scale_colour_manual(values = gg_color_hue(6)) +
  ylab("Monthly Catch") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

cluster <- readRDS("data/systematic_review_results/cluster_membership.rds") %>%
  select(id, cluster)
overall_pv3_subset <- overall_pv3 %>%
  left_join(cluster, by = "id") %>%
  filter(cluster == 1)
  #filter(id %in% c(6, 19, 22, 36, 37, 52, 55, 63, 64, 65, 73))

overall_pv4_subset <- overall_pv3_subset %>%
  select(id, timepoint, norm_vector_dens)
overall_pv4_subset2 <- overall_pv3_subset %>%
  select(id, lagged_timepoint, airtemp, rainfall)
overall_pv5 <- overall_pv4_subset %>%
  left_join(overall_pv4_subset2, by = c("timepoint" = "lagged_timepoint"))

ggplot(data = overall_pv5) +
  #geom_point(aes(x = airtemp, y = norm_vector_dens)) +
  geom_smooth(aes(x = airtemp, y = norm_vector_dens))
ggplot(data = overall_pv5) +
  #geom_point(aes(x = rainfall, y = norm_vector_dens)) +
  geom_smooth(aes(x = rainfall, y = norm_vector_dens)) 



##############################################
ggplot(data = overall_pv3_subset) +
  geom_line(aes(x = timepoint, y = norm_vector_dens, col = factor(id)), size = 2)



plot(overall_pv3_subset$airtemp,
     overall_pv3_subset$norm_vector_dens)

# ggplot(data = overall_pv3_subset) +
#   geom_path(aes(x = timepoint, y = mean, col = factor(country)), size = 2) +
#   geom_point(aes(x = timepoint, y = raw_catch, group = factor(id)), col = "black") +
#   facet_wrap(~id, scales = "free_y") +
#   scale_y_continuous(limits=c(0, NA), position = "left",
#                      sec.axis = sec_axis(trans=~int + ./coef, name = "Surface Temp")) +
#   scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
#                                 "J", "A", "S", "O", "N", "D"),
#                      breaks = seq(2, 24, length.out = 12)) +
#   geom_line(aes(x = timepoint, y = coef * (airtemp - int))) +
#   scale_fill_manual(values = gg_color_hue(6)) +
#   scale_colour_manual(values = gg_color_hue(6)) +
#   ylab("Monthly Catch") +
#   theme_bw() +
#   theme(legend.position = "bottom",
#     axis.title.x = element_blank(),
#     strip.background = element_blank(),
#     strip.text = element_blank(),
#     legend.justification = c(1, 0)) 

ggplot(data = overall_pv2_subset) +
  geom_path(aes(x = timepoint, y = airtemp), size = 2) +
  facet_wrap(~id, scales = "free_y")



summary <- overall_pv3 %>%
  group_by(airtemp_thresh) %>%
  summarise(mean = mean(norm_vector_dens, na.rm = TRUE),
            n = n())

mean_airtemp_dens <- overall_pv3 %>%
  filter(!is.na(airtemp)) %>%
  group_by(bin_airtemp, bin_airtemp2) %>%
  summarise(mean_norm_vd = mean(norm_vector_dens))
ggplot(mean_airtemp_dens, aes(x = bin_airtemp2, y = mean_norm_vd)) +
  geom_line()


plot(overall_pv2$airtemp, overall_pv2$norm_vector_dens)
  
norm_rainfall_storage <- normalise_total(rainfall_storage)
smoothed_rainfall <- t(apply(norm_rainfall_storage, 1, raster::movingFun, n = 4, circular = TRUE))
rainfall_start_index <- apply(smoothed_rainfall, 1, function(x) which(x == max(x)))
rainfall_seas_3 <- apply(norm_rainfall_storage, 1, percent_incidence, 3, 2)
rainfall_seas_4 <- apply(norm_rainfall_storage, 1, percent_incidence, 4, 2)




# Reordering mean fitted time-series to start at max
reordered_mean_realisation <- matrix(nrow = length(retain_index), ncol = (12 * interpolating_points + 1))
start_index <- apply(mean_realisation, 1, function(x) which(x == max(x)))
end_index <- dim(reordered_mean_realisation)[2]
for (i in 1:length(retain_index)) {
  reordered_mean_realisation[i, ] <- mean_realisation[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
}


# Reordering rainfall time-series to start at max of vector density (note NOT max rainfall, which is what I calculate above)
reordered_rainfall <- matrix(nrow = length(retain_index), ncol = length(months_length))
end_rain_index <- dim(reordered_rainfall)[2]
for (i in 1:length(retain_index)) {
  reordered_rainfall[i, ] <- norm_rainfall_storage[i, c(start_index[i]:end_rain_index, 1:(start_index[i]-1))]
}