# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(scales); library(patchwork); library(here); library(zoo);
library(tictoc); library(ggnewscale)

# Loading in custom functions
options(scipen = 999)
source(here("functions", "time_series_characterisation_functions.R"))
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
source(here::here("functions", "time_series_characterisation_functions.R"))

# Loading functions and mosquito bionomics data
overall <- readRDS(file = here::here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

# Loading and processing different IRS parameters
actellic <- read.csv(here::here("data", "IRS_parameters", "Actellic_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

bendiocarb <- read.csv(here::here("data", "IRS_parameters", "Bendiocarb_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

# Extracting Mean Realisation for Each Time-Series
set.seed(10)
urban_rural <- overall$city
interpolating_points <- 2
mean_realisation <- matrix(nrow = dim(overall)[1], ncol = (12 * interpolating_points + 1))
prior <- "informative"
for (i in 1:length(overall$id)) {
  # Indexing 
  index <- overall$id[i]
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here::here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here::here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
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
}
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Converting stephensi profiles to be 365 days in length and generating estimates of seasonality
scaled_year_output <- t(apply(normalised_output, 1, function(x) {
  temp <- approx(x, n = 365)$y
  temp <- 365 * temp/sum(temp)
  temp[temp < 0.1] <- 0.1
  temp }))
steph_seasonality_list <- vector(mode = "list", length = dim(scaled_year_output)[1])
for (i in 1:(dim(scaled_year_output)[1])) {
  steph_seasonality_list[[i]] <- scaled_year_output[i, ]
}
saveRDS(steph_seasonality_list, file = here::here("data/steph_seasonality_list.rds"))

seasonality <- c()
for (i in 1:(length(steph_seasonality_list))) {
  seasonality <- c(seasonality, calc_incidence_seasonality(steph_seasonality_list[[i]], 3))
}
seas_df <- data.frame(seasonal_profile = 1:length(steph_seasonality_list), seasonality = seasonality)

# Extracting Rainfall Data (also sort out leap year stuff)
timings <- seq(0, 365, 15)
months_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
monthly_rainfall_storage <- matrix(nrow = dim(overall)[1], ncol = length(months_length)) 
daily_rainfall_storage <- matrix(nrow = dim(overall)[1], ncol = 365) 
for (i in 1:length(overall$id)) {
  index <- overall$id[i]
  temp <- c()
  rf <- read.csv(here(paste0("data/location_specific_rainfall/rainfall_ts", index, ".csv")))
  rf <- rf %>%
    dplyr::group_by(daymonth_id) %>%
    dplyr::summarise(rainfall = mean(rainfall))
  daily_rainfall_storage[i, ] <- rf$rainfall[1:365]
  counter <- 1
  count_vec <- counter
  for (j in 1:length(months_length)) {
    indices <- counter:(counter + months_length[j] - 1)
    temp <- c(temp, sum(rf$rainfall[indices]))
    counter <- counter + months_length[j]
  }
  monthly_rainfall_storage[i, ] <- temp
}
peak_rainfall_month <- apply(monthly_rainfall_storage, 1, function(x) which(x == max(x)))
peak_rainfall_fortnight <- peak_rainfall_month * 2
peak_rainfall_fortnight <- data.frame(seasonal_profile = 1:65, rainfall_peak_timing = peak_rainfall_fortnight)

rainfall_list <- vector(mode = "list", length = dim(scaled_year_output)[1])
for (i in 1:(dim(scaled_year_output)[1])) {
  rainfall_list[[i]] <- daily_rainfall_storage[i, ]
}

# Setting Average Mosquito Density
density <- 20
years <- 20
density_vec <- rep(density, 365 * years)

# Intervention Timings
timings <- seq(0, 365, 15)

# Loading in IRS results and adding in id for seasonal profile
act_mod <- readRDS(here::here("outputs", "malaria_model_IRS_running_actelic.rds"))
ben_mod <- readRDS(here::here("outputs", "malaria_model_IRS_running_bendiocarb.rds"))
for (i in 1:length(steph_seasonality_list)) {
  act_mod[[i]]$seasonal_profile <- i
  act_mod[[i]] <- act_mod[[i]][-1, ]
  ben_mod[[i]]$seasonal_profile <- i
  ben_mod[[i]] <- ben_mod[[i]][-1, ]
}

# Loading in counterfactual results (i.e. no IRS)
counterfactual <- readRDS(here("outputs", "malaria_model_IRS_running_counterfactuals.rds"))
counterfactual <- bind_rows(counterfactual) %>%
  dplyr::rename(counterfactual_incidence = incidence, counterfactual_prevalence = prevalence)

# Joining IRS and no IRS results
start <- ((years - 2) * 365)
end <- max(act_mod[[1]]$t)
act_mod2 <- bind_rows(act_mod) %>%
  left_join(counterfactual, by = c("seasonal_profile", "t")) %>%
  filter(t >= start)
ben_mod2 <- bind_rows(ben_mod) %>%
  left_join(counterfactual, by = c("seasonal_profile", "t"))  %>%
  filter(t >= start)

# Generating annual reducation in incidence for each IRS, timing of delivery and seasonal profile combo
reduction <- data.frame(insecticide = "x", seasonal_profile = -1, timing = -1, count_inc = -1, inc = -1, inc_red = -1, inc_prop_red = -1)
for (i in 1:length(steph_seasonality_list)) {
  
  for (j in 1:length(timings)) {
    
    act_temp <- act_mod2[act_mod2$seasonal_profile == i & act_mod2$timing == timings[j], ]
    ben_temp <- ben_mod2[act_mod2$seasonal_profile == i & ben_mod2$timing == timings[j], ]

    sum_start <- ifelse(timings[j] == 0, 1, timings[j])
    sum_end <- sum_start + 365
    
    counterfactual_incidence <- sum(act_temp$counterfactual_incidence[sum_start:sum_end])
    act_incidence <- sum(act_temp$incidence[sum_start:sum_end])
    ben_incidence <- sum(ben_temp$incidence[sum_start:sum_end])

    red_temp <- data.frame(insecticide = c("pirim", "ben"), 
                           seasonal_profile = i, 
                           timing = j, 
                           count_inc = counterfactual_incidence,
                           inc = c(act_incidence, ben_incidence),
                           inc_red = counterfactual_incidence - c(act_incidence, ben_incidence), 
                           inc_prop_red = (counterfactual_incidence - c(act_incidence, ben_incidence))/counterfactual_incidence)
    
    reduction <- rbind(reduction, red_temp)
    
  }
  print(i)
}
reduction <- reduction[-1, ]

# Extracting the maximum reduction in incidence achieved for each insecticide and temporal profile
reduction$insecticide <- factor(reduction$insecticide)
max_reduction <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  mutate(max_inc_red = max(inc_red),
         max_inc_prop_red = max(inc_prop_red)) %>%
  ungroup() %>%
  filter(inc_red == max_inc_red) %>%
  dplyr::select(-max_inc_red, -max_inc_prop_red)
max_reduction$timing_type <- "maximum"

# Extracting the reduction in incidence achieved if timed based on rainfall peak
rainfall_reduction <- reduction %>%
  dplyr::left_join(peak_rainfall_fortnight, by = "seasonal_profile") %>%
  dplyr::group_by(insecticide, seasonal_profile) %>%
  dplyr::filter(timing == rainfall_peak_timing) %>%
  dplyr::select(-rainfall_peak_timing)
rainfall_reduction$timing_type <- "rainfall"

# Calculating the average incidence reduction across all timed deliveries
average_reduction <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  dplyr::summarise(timing = 100,
                   count_inc = mean(count_inc),
                   inc = mean(inc),
                   inc_red = count_inc - inc,
                   inc_prop_red = inc_red/count_inc)
average_reduction$timing_type <- "annual_average"

# Generating overall dataframe with all the results 
raw_overall_df <- rbind(max_reduction, rainfall_reduction, average_reduction)

# Combining seasonality and reduction in malaria incidence together, and plotting
overall_df <- raw_overall_df %>%
  left_join(seas_df, by = "seasonal_profile")
overall_df$timing_type <- factor(overall_df$timing_type, levels = c("annual_average", "rainfall", "maximum"))

# Plotting Maximum (Optimal Timing) Overall Reduction vs Seasonality and IRS Profiles
time_since_IRS <- 1:365
actellic_output <- 1 / (1 + exp(-(actellic$irs_decay_mort1 + actellic$irs_decay_mort2 * time_since_IRS)))
bendiocarb_output <- 1 / (1 + exp(-(bendiocarb$irs_decay_mort1 + bendiocarb$irs_decay_mort2 * time_since_IRS)))

irs_killing <- data.frame(time_since_IRS = time_since_IRS, 
                          Pirimiphos_Methyl = actellic_output,
                          Bendiocarb = bendiocarb_output) %>%
  pivot_longer(-time_since_IRS, names_to = "insecticide", values_to = "mortality_effect")
insecticide_names <- c(`Pirimiphos_Methyl` = "Pirimiphos Methyl",
                       `Bendiocarb` = "Bendiocarb")
irs_killing_plot <- ggplot(irs_killing, aes(x = time_since_IRS, y = mortality_effect, colour = insecticide)) +
  geom_path(size = 2) +
  theme_bw() +
  scale_colour_manual(values = c("#ae337a", "#209e88"),
                      labels = c("Bendiocarb", "Pirimiphos Methyl")) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~insecticide,
             labeller = as_labeller(insecticide_names)) +
  labs(x = "Time Since IRS Spraying (Days)", y = "Probability of Mosquitoes Dying") +
  theme_bw() +
  theme(legend.justification = c(1, 0), #
        strip.background = element_rect(fill = "white"),
        axis.title.y = element_text(size = 8),
        legend.position = "none") + #c(0.97, 0.65)) +
  guides(col = guide_legend(ncol = 1, title = "Insecticide"))

data_overall_bend <- overall_df[overall_df$timing_type == "maximum" & overall_df$insecticide == "ben", ]
data_overall_pirim <- overall_df[overall_df$timing_type == "maximum" & overall_df$insecticide == "pirim", ]
data_overall_bend_rand <- overall_df[overall_df$timing_type == "annual_average" & overall_df$insecticide == "ben", ]
data_overall_pirim_rand <- overall_df[overall_df$timing_type == "annual_average" & overall_df$insecticide == "pirim", ]

irs_overall_reduction <- ggplot() +
  geom_point(data = data_overall_bend, aes(x = seasonality, y = 100 * inc_prop_red, col = seasonality)) +
  geom_segment(aes(x = data_overall_bend$seasonality,
                   xend = data_overall_bend$seasonality,
                   y = 100 * data_overall_bend_rand$inc_prop_red,
                   yend = 100 * data_overall_bend$inc_prop_red - 0.5, 
                   col = data_overall_bend$seasonality),
               arrow=arrow(type="closed", length = unit(0.02, "npc"))) +
  geom_point(data = data_overall_bend_rand, aes(x = seasonality, y = 100 * inc_prop_red), col = "grey") +
  theme_bw() +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  new_scale_color() +
  geom_point(data = data_overall_pirim, aes(x = seasonality, y = 100 * inc_prop_red, col = seasonality)) +
  geom_segment(aes(x = data_overall_pirim$seasonality,
                   xend = data_overall_pirim$seasonality,
                   y = 100 * data_overall_pirim_rand$inc_prop_red,
                   yend = 100 * data_overall_pirim$inc_prop_red - 0.5, 
                   col = data_overall_pirim$seasonality),
               arrow=arrow(type="closed", length = unit(0.02, "npc"))) +
  geom_point(data = data_overall_pirim_rand, aes(x = seasonality, y = 100 * inc_prop_red), col = "grey") +
  theme_bw() +
  scale_colour_viridis(option = "viridis", limits = c(0.25, 1), name = "Seasonality") +
  theme(legend.position = "none") +
  labs(x = "% of Annual Incidence in 3 Months", y = "% Reduction in Incidence")

irs_overall <- plot_grid(irs_killing_plot, irs_overall_reduction, nrow = 2, rel_heights = c(0.3, 1))

# Plotting Bendiocarb Comparative Impact Depending On Timing
ben_data <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type != "maximum", ] %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
ben_data$timing_type_new <- ifelse(ben_data$timing_type == "rainfall", "rainfall_comp", "average_comp")
ben_data2 <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type == "maximum", ]  %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
ben_data2$timing_type_new <- "rainfall_comp"
ben_data3 <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type == "maximum", ]  %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
ben_data3$timing_type_new <- "average_comp"
ben_data_overall <- rbind(ben_data, ben_data2, ben_data3)

set.seed(10)
ben_impact <- ggplot() +
  geom_boxplot(data = ben_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, group = new_timing),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = ben_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  facet_wrap(~timing_type_new) +
  geom_line(data = ben_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, group = seasonal_profile, 
                             col = seasonality), position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "right", #plot.margin = unit(c(1,1,1,1), "cm"),
        strip.background = element_blank(), strip.text = element_blank()) +
  lims(y = c(9, 37.3)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Random\nMonth", "_____\nTiming")) +
  scale_y_continuous(limits = c(9, 37.3), position = "right")

# Plotting Pirimithos Methyl Comparative Impact Depending On Timing
pirim_data <- overall_df[overall_df$insecticide == "pirim" & overall_df$timing_type != "maximum", ] %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
pirim_data$timing_type_new <- ifelse(pirim_data$timing_type == "rainfall", "rainfall_comp", "average_comp")
pirim_data2 <- overall_df[overall_df$insecticide == "pirim" & overall_df$timing_type == "maximum", ]  %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
pirim_data2$timing_type_new <- "rainfall_comp"
pirim_data3 <- overall_df[overall_df$insecticide == "pirim" & overall_df$timing_type == "maximum", ]  %>%
  dplyr::select(seasonal_profile, inc_prop_red, timing_type, seasonality) %>%
  mutate(new_timing = ifelse(timing_type == "maximum", 2, 1)) 
pirim_data3$timing_type_new <- "average_comp"
pirim_data_overall <- rbind(pirim_data, pirim_data2, pirim_data3)

set.seed(10)
pirim_impact <- ggplot() +
  geom_boxplot(data = pirim_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, group = new_timing),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = pirim_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  facet_wrap(~timing_type_new) +
  geom_line(data = pirim_data_overall, aes(x = new_timing, y = 100 * inc_prop_red, group = seasonal_profile, 
                                         col = seasonality), position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "viridis", limits = c(0.25, 1), name = "Seasonality") +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "right", #plot.margin = unit(c(1,1,1,1), "cm"),
        strip.background = element_blank(), strip.text = element_blank()) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Random\nMonth", "_____\nTiming")) +
  scale_y_continuous(limits = c(32.5, 50), position = "right")

second_irs_figure <- plot_grid(pirim_impact, ben_impact, nrow = 2)
final_irs_figure <- plot_grid(irs_overall, second_irs_figure, ncol = 2, rel_widths = c(1, 1))
final_irs_figure

ggsave(filename = here("figures", "Fig4_Overall.pdf"), 
                       plot = final_irs_figure, 
                       height = 9.1, 
                       width = 11.48)
