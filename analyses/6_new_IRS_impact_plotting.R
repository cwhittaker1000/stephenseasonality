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

# Plotting Overall Impact, Annual Average vs Optimal Timing - Bendiocarb
set.seed(10)
data <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type != "rainfall", ]
a <- ggplot() +
  geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  scale_x_discrete(labels = c("Random\nMonth", "Optimal\nTiming")) +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  lims(y = c(18.7, 37.3))

data2 <- data %>%
  dplyr::select(seasonal_profile, seasonality, inc_prop_red, timing_type) %>%
  pivot_wider(id_cols = c("seasonal_profile", "seasonality"), names_from = timing_type, values_from = inc_prop_red) %>%
  mutate(diff = maximum - annual_average,
         diff_prop = diff/annual_average,
         inc_prop_red = maximum,
         start = 0)
b <- ggplot() +
  geom_bar(data = data2, aes(x = inc_prop_red * 100, y = diff_prop * 100, fill = seasonality,
                             group = seasonal_profile), alpha = 0.5, stat = "identity", width = 0.4, col = "black", size = 0.2) +
  scale_fill_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  coord_flip() +
  theme_bw() +
  labs(y = "% Increased Impact") +
  theme(legend.position = "right",
        axis.line = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = +7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,-0.1), "cm")) +
  scale_x_continuous(position = "bottom", limits = c(18.7, 37.3),
                     breaks = c(20, 25, 30, 35))
c <- plot_grid(a, b, nrow = 1, rel_widths = c(0.8, 0.6), align = "h") # 5.8 * 5.2 looks nice

# Plotting Overall Impact, Annual Average vs Optimal Timing - Bendiocarb
set.seed(10)
data <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type != "annual_average", ]
d <- ggplot() +
  geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  scale_x_discrete(labels = c("Rainfall\nTiming", "Optimal\nTiming")) +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  lims(y = c(9, 37.3))

data2 <- data %>%
  dplyr::select(seasonal_profile, seasonality, inc_prop_red, timing_type) %>%
  pivot_wider(id_cols = c("seasonal_profile", "seasonality"), names_from = timing_type, values_from = inc_prop_red) %>%
  mutate(diff = maximum - rainfall,
         diff_prop = diff/rainfall,
         inc_prop_red = maximum,
         start = 0)
e <- ggplot() +
  geom_bar(data = data2, aes(x = inc_prop_red * 100, y = diff_prop * 100, fill = seasonality,
                             group = seasonal_profile), alpha = 0.5, stat = "identity", width = 0.4, col = "black", size = 0.2) +
  scale_fill_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  coord_flip() +
  theme_bw() +
  labs(y = "% Increased Impact") +
  theme(legend.position = "right",
        axis.line = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = +7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,-0.1), "cm")) +
  scale_x_continuous(position = "bottom", limits = c(9, 37.3),
                     breaks = c(20, 25, 30, 35))
f <- plot_grid(d, e, nrow = 1, rel_widths = c(0.8, 0.6), align = "h") # 5.8 * 5.2 looks nice

g <- plot_grid(c, f, nrow = 1)

# Plotting Overall Impact, Annual Average vs Optimal Timing - Pirimithos
set.seed(10)
data <- overall_df[overall_df$insecticide == "pirim" & overall_df$timing_type != "rainfall", ]
h <- ggplot() +
  geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  scale_x_discrete(labels = c("Random\nMonth", "Optimal\nTiming")) +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  lims(y = c(40.5, 50.2))

data2 <- data %>%
  dplyr::select(seasonal_profile, seasonality, inc_prop_red, timing_type) %>%
  pivot_wider(id_cols = c("seasonal_profile", "seasonality"), names_from = timing_type, values_from = inc_prop_red) %>%
  mutate(diff = maximum - annual_average,
         diff_prop = diff/annual_average,
         inc_prop_red = maximum,
         start = 0)
i <- ggplot() +
  geom_bar(data = data2, aes(x = inc_prop_red * 100, y = diff_prop * 100, fill = seasonality,
                             group = seasonal_profile), alpha = 0.5, stat = "identity", width = 0.4, col = "black", size = 0.2) +
  scale_fill_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  coord_flip() +
  theme_bw() +
  labs(y = "% Increased Impact") +
  theme(legend.position = "right",
        axis.line = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = +7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,-0.1), "cm")) +
  scale_x_continuous(position = "bottom", limits = c(40.5, 50.2),
                     breaks = c(20, 25, 30, 35))
j <- plot_grid(h, i, nrow = 1, rel_widths = c(0.8, 0.6), align = "h") # 5.8 * 5.2 looks nice

# Plotting Overall Impact, Annual Average vs Optimal Timing - Pirimithos
set.seed(10)
data <- overall_df[overall_df$insecticide == "pirim" & overall_df$timing_type != "annual_average", ]
k <- ggplot() +
  geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  scale_x_discrete(labels = c("Rainfall\nTiming", "Optimal\nTiming")) +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  lims(y = c(32.5, 50.2))

data2 <- data %>%
  dplyr::select(seasonal_profile, seasonality, inc_prop_red, timing_type) %>%
  pivot_wider(id_cols = c("seasonal_profile", "seasonality"), names_from = timing_type, values_from = inc_prop_red) %>%
  mutate(diff = maximum - rainfall,
         diff_prop = diff/rainfall,
         inc_prop_red = maximum,
         start = 0)
l <- ggplot() +
  geom_bar(data = data2, aes(x = inc_prop_red * 100, y = diff_prop * 100, fill = seasonality,
                             group = seasonal_profile), alpha = 0.5, stat = "identity", width = 0.4, col = "black", size = 0.2) +
  scale_fill_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  coord_flip() +
  theme_bw() +
  labs(y = "% Increased Impact") +
  theme(legend.position = "right",
        axis.line = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = +7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,-0.1), "cm")) +
  scale_x_continuous(position = "bottom", limits = c(32.5, 50.2),
                     breaks = c(20, 25, 30, 35))
m <- plot_grid(k, l, nrow = 1, rel_widths = c(0.8, 0.6), align = "h") # 5.8 * 5.2 looks nice

n <- plot_grid(j, m, nrow = 1)

o <- plot_grid(g, n, nrow = 2)

#####

set.seed(10)
IRS_max_red_boxplot <- ggplot(overall_df, aes(x = insecticide, y = 100 * inc_prop_red, col = insecticide)) +
  geom_boxplot(fill = NA, outlier.shape = NA) + 
  scale_colour_manual(values = c("#FA9F42", "#0B6E4F", "#B2B2B3"),
                      labels = c("Bendiocarb", "Pirimiphos Methyl")) +
  geom_jitter(aes(x = insecticide, y = 100 * inc_prop_red), size = 1, width = 0.25) +
  scale_x_discrete(labels = c("Bendiocarb", "Pirimiphos Methyl")) +
  facet_grid(~timing_type) +
  theme_bw() +
  labs(y = "% Reduction In Incidence") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")




ggplot() +
  geom_boxplot(data = data2, aes(y = 100 * diff_prop), fill = NA, alpha = 0.5, col = "grey") +
  geom_point(data = data2, aes(x = 0, y = 100 * diff_prop, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") 




ggplot() +
  geom_boxplot(data = data2, aes(y = 100 * diff_prop), fill = NA, alpha = 0.5, col = "grey") +
  geom_point(data = data2, aes(x = 0, y = 100 * diff_prop, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") 



data <- overall_df[overall_df$insecticide == "ben" & overall_df$timing_type != "annual_average", ]
ggplot() +
  geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
               fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
  new_scale_color() +
  geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
             size = 2, position = position_dodge(0.3)) +
  geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            position = position_dodge(0.3), alpha = 0.7) +
  scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
  scale_x_discrete(labels = c("Rainfall\nTiming", "Optimal\nTiming")) +
  #facet_wrap(~insecticide, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "% Reduction in Incidence")


ggplot() +
  geom_boxplot(data = data2, aes(y = 100 * diff, col = seasonality), fill = NA, outlier.shape = NA) +
  geom_point(data = data2, aes(x = 0, y = 100 * diff, col = seasonality)) 
plot(data2$seasonality, data2$diff)


data <- overall_df[overall_df$insecticide == "clothi" &
                     overall_df$timing_type != "annual_average", ]
data2 <- data %>%
  dplyr::select(seasonal_profile, seasonality, inc_prop_red, timing_type) %>%
  pivot_wider(id_cols = c("seasonal_profile", "seasonality"), names_from = timing_type, values_from = inc_prop_red) %>%
  mutate(diff = maximum - rainfall)

ggplot() +
  geom_boxplot(data = data2, aes(y = 100 * diff, col = seasonality), fill = NA, outlier.shape = NA) +
  geom_point(data = data2, aes(x = 0, y = 100 * diff, col = seasonality)) 
plot(data2$seasonality, data2$diff)

summary(lm(data2$diff ~ 0 + data2$seasonality))

ggplot(overall_df[overall_df$insecticide == "clothi" & overall_df$timing_type != "annual_average", ]) +
  geom_line(aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            alpha = 0.2, position = position_dodge(0.1)) +
  geom_boxplot(aes(x = timing_type, y = 100 * inc_prop_red, fill = timing_type, col = timing_type),
               fill = NA, outlier.shape = NA) + 
  geom_point(aes(x = timing_type, y = 100 * inc_prop_red, col = timing_type, group = seasonal_profile), 
             size = 1, position = position_dodge(0.1)) 

ggplot(overall_df[overall_df$insecticide == "clothi" & overall_df$timing_type != "annual_average", ]) +
  geom_line(aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
            alpha = 0.2, position = position_dodge(0.1)) +
  geom_point(aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, fill = timing_type), 
             size = 1, position = position_dodge(0.1)) +
  geom_boxplot(aes(x = timing_type, y = 100 * inc_prop_red, fill = timing_type),
               fill = NA, outlier.shape = NA) 



# ggplot() +
#   geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
#                fill = NA, outlier.shape = NA, alpha = 0.5, col = "grey") + 
#   new_scale_color() +
#   geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
#              size = 2, position = position_dodge(0.3)) +
#   geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
#             position = position_dodge(0.3), alpha = 0.7) +
#   scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
#   scale_x_discrete(labels = c("Random\nMonth", "Optimal\nTiming")) +
#   #facet_wrap(~insecticide, scales = "free_y") +
#   theme_bw() +
#   labs(x = "", y = "% Reduction in Incidence")
# data <- overall_df[overall_df$insecticide == "ben", ]
# data$timing_type <- factor(data$timing_type, levels = c("annual_average", "maximum", "rainfall"))
# ggplot() +
#   geom_boxplot(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = timing_type),
#                fill = NA, outlier.shape = NA, alpha = 0.25, col = "grey") + 
#   new_scale_color() +
#   geom_point(data = data, aes(x = timing_type, y = 100 * inc_prop_red, col = seasonality, group = seasonal_profile),
#              size = 2, position = position_dodge(0.3)) +
#   geom_line(data = data, aes(x = timing_type, y = 100 * inc_prop_red, group = seasonal_profile, col = seasonality),
#             position = position_dodge(0.3), alpha = 0.7) +
#   scale_colour_viridis(option = "magma", limits = c(0.25, 1), name = "Seasonality") +
#   scale_x_discrete(labels = c("Random\nMonth", "Optimal\nTiming")) +
#   #facet_wrap(~insecticide, scales = "free_y") +
#   theme_bw() +
#   labs(x = "", y = "% Reduction in Incidence")

# b <- ggplot() +
#   geom_point(data = data2, aes(x = diff_prop * 100, y = inc_prop_red * 100,
#                                col = seasonality)) +
#   geom_segment(data = data2, aes(y = inc_prop_red * 100, yend = inc_prop_red * 100, x = 0, xend = diff_prop * 100, 
#                                  col = seasonality)) +
#   scale_y_continuous(position = "left", limits = c(18.7, 37.1)) +
#   scale_colour_viridis(option = "magma", limits = c(0.25, 1), 
#                        name = "Seasonality") +
#   theme_bw() +
#   labs(x = "% Increase\nin Impact") +
#   theme(legend.position = "none",
#         axis.line = element_blank(), #element_line(colour = "black"),
#         axis.line.y = element_line(size = 0.25),
#         axis.line.x = element_line(size = 0.25),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(), 
#         axis.title.y = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         #panel.border = element_blank(),
#         #panel.background = element_blank(),
#         plot.margin = unit(c(1,1,0,-0.11), "cm")) +
#   scale_x_continuous(expand = c(0, 0))
# plot_grid(a, b, nrow = 1, rel_widths = c(1.5, 0.8), align = "h")

