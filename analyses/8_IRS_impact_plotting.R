# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(scales); library(patchwork); library(here); library(zoo);
library(tictoc)

# Helper function 
calc_incidence_seasonality <- function(input, num_months) {
  days_in_year <- length(input)
  period <- round(days_in_year * num_months/12)
  incidence_vector <- vector(mode = "numeric", length = days_in_year)
  for (i in 1:days_in_year) {
    index <- i:(i+period-1)
    if (sum(index > days_in_year) > 0) {
      index[index>days_in_year] <- index[index>days_in_year] - days_in_year
    }
    incidence_vector[i] <- sum(input[index])/sum(input)
  }
  return(max(incidence_vector))
}

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
source(here::here("functions", "time_series_characterisation_functions.R"))
overall <- readRDS(file = here::here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

# Loading and processing different IRS parameters
actellic <- read.csv(here("data", "IRS_parameters", "Actellic_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

bendiocarb <- read.csv(here("data", "IRS_parameters", "Bendiocarb_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

sumishield <- read.csv(here("data", "IRS_parameters", "Sumishield_uncertainty.csv")) %>%
  dplyr::rename(irs_decay_mort1 = Mort1, irs_decay_mort2 = Mort2, irs_decay_succ1 = Succ1, 
                irs_decay_succ2 = Succ2, irs_decay_det1 = Det1, irs_decay_det2 = Det2) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

# Plotting different IRS killing effects over time
time_since_IRS <- 1:365
actellic_output <- 1 / (1 + exp(-(actellic$irs_decay_mort1 + actellic$irs_decay_mort2 * time_since_IRS)))
bendiocarb_output <- 1 / (1 + exp(-(bendiocarb$irs_decay_mort1 + bendiocarb$irs_decay_mort2 * time_since_IRS)))
sumishield_output <- 1 / (1 + exp(-(sumishield$irs_decay_mort1 + sumishield$irs_decay_mort2 * time_since_IRS)))

irs_killing <- data.frame(time_since_IRS = time_since_IRS, Pirimiphos_Methyl = actellic_output,
                          Bendiocarb = bendiocarb_output, Clothiandin = sumishield_output) %>%
  pivot_longer(-time_since_IRS, names_to = "insecticide", values_to = "mortality_effect")

irs_killing_plot <- ggplot(irs_killing, aes(x = time_since_IRS, y = mortality_effect, colour = insecticide)) +
  geom_path(size = 2) +
  theme_bw() +
  scale_colour_manual(values = c("#FA9F42", "#0B6E4F", "#B2B2B3"),
                      labels = c("Bendiocarb", "Clothiandin", "Pirimiphos Methyl")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time Since IRS Spraying (Days)", y = "Probability of Mosquitoes Dying") +
  theme_bw() +
  theme(legend.justification = c(1, 0), #
        legend.position = "none") +
        #legend.position = c(1, 0.75)) +
  guides(col = guide_legend(ncol = 1, title = "Insecticide"))
  
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

# Converting stephensi profiles to be 365 days in length
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

# Setting Average Mosquito Density
density <- 20
years <- 20
density_vec <- rep(density, 365 * years)

# Intervention Timings
timings <- seq(0, 365, 15)

# Loading in IRS results and adding in id for seasonal profile
act_mod <- readRDS(here::here("outputs", "malaria_model_IRS_running_actelic.rds"))
ben_mod <- readRDS(here::here("outputs", "malaria_model_IRS_running_bendiocarb.rds"))
sumi_mod <- readRDS(here::here("outputs", "malaria_model_IRS_running_sumishield.rds"))
for (i in 1:length(steph_seasonality_list)) {
  act_mod[[i]] <- act_mod[[i]][-1, ]
  act_mod[[i]]$seasonal_profile <- i
  ben_mod[[i]]$seasonal_profile <- i
  ben_mod[[i]] <- ben_mod[[i]][-1, ]
  sumi_mod[[i]]$seasonal_profile <- i
  sumi_mod[[i]] <- sumi_mod[[i]][-1, ]
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
sumi_mod2 <- bind_rows(sumi_mod) %>%
  left_join(counterfactual, by = c("seasonal_profile", "t")) %>%
  filter(t >= start)

# Extracting maximum reduction achieved for each IRS and seasonal profile combo
reduction <- data.frame(insecticide = "x", seasonal_profile = -1, timing = -1, count_inc = -1, inc = -1, inc_red = -1, inc_prop_red = -1)
for (i in 1:length(steph_seasonality_list)) {
  
  for (j in 1:length(timings)) {
    
    act_temp <- act_mod2[act_mod2$seasonal_profile == i & act_mod2$timing == timings[j], ]
    ben_temp <- ben_mod2[act_mod2$seasonal_profile == i & ben_mod2$timing == timings[j], ]
    sumi_temp <- sumi_mod2[act_mod2$seasonal_profile == i & sumi_mod2$timing == timings[j], ]
    
    sum_start <- ifelse(timings[j] == 0, 1, timings[j])
    sum_end <- sum_start + 365
    
    counterfactual_incidence <- sum(act_temp$counterfactual_incidence[sum_start:sum_end])
    act_incidence <- sum(act_temp$incidence[sum_start:sum_end])
    ben_incidence <- sum(ben_temp$incidence[sum_start:sum_end])
    sumi_incidence <- sum(sumi_temp$incidence[sum_start:sum_end])
    
    red_temp <- data.frame(insecticide = c("pirim", "ben", "clothi"), 
                           seasonal_profile = i, 
                           timing = j, 
                           count_inc = counterfactual_incidence,
                           inc = c(act_incidence, ben_incidence, sumi_incidence),
                           inc_red = counterfactual_incidence - c(act_incidence, ben_incidence, sumi_incidence), 
                           inc_prop_red = (counterfactual_incidence - c(act_incidence, ben_incidence, sumi_incidence))/counterfactual_incidence)
    
    reduction <- rbind(reduction, red_temp)
    
  }
  print(i)
}
reduction <- reduction[-1, ]
reduction$insecticide <- factor(reduction$insecticide)
max_red <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  mutate(max_inc_red = max(inc_red),
         max_inc_prop_red = max(inc_prop_red)) %>%
  ungroup() %>%
  filter(inc_red == max_inc_red)

# Generating estimates of profile seasonality 
seasonality <- c()
for (i in 1:(length(steph_seasonality_list))) {
  seasonality <- c(seasonality, calc_incidence_seasonality(steph_seasonality_list[[i]], 3))
}
seas_df <- data.frame(seasonal_profile = 1:length(steph_seasonality_list), seasonality = seasonality)

# Combining seasonality and reduction in malaria incidence together, and plotting
max_red2 <- max_red %>%
  left_join(seas_df, by = "seasonal_profile")
set.seed(10)
IRS_max_red_boxplot <- ggplot(max_red2, aes(x = insecticide, y = 100 * max_inc_prop_red, col = insecticide)) +
  geom_boxplot(fill = NA, outlier.shape = NA) + 
  scale_colour_manual(values = c("#FA9F42", "#0B6E4F", "#B2B2B3"),
                      labels = c("Bendiocarb", "Clothiandin", "Pirimiphos Methyl")) +
  geom_jitter(aes(x = insecticide, y = 100 * max_inc_prop_red), size = 1, width = 0.25) +
  scale_x_discrete(labels = c("Bendiocarb", "Clothiandin", "Pirimiphos Methyl")) +
  theme_bw() +
  labs(y = "Maximum % Reduction In Incidence") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        strip.background = element_blank(), strip.placement = "outside")

IRS_max_red_season_scatterplot <- ggplot(max_red2, aes(x = seasonality, y = 100 * max_inc_prop_red, col = insecticide)) +
  geom_point() + 
  scale_colour_manual(values = c("#FA9F42", "#0B6E4F", "#B2B2B3"),
                      labels = c("Bendiocarb", "Clothiandin", "Pirimiphos Methyl")) +
  theme_bw() +
  labs(y = "Maximum % Reduction In Incidence",
       x = "% of Annual Incidence In 3 Months") +
  theme(legend.position = "none",
        legend.justification = c(1, 0)) +#
        #legend.position = c(1, 0.05)) +
  guides(col = guide_legend(ncol = 1, title = "Insecticide", override.aes = list(size = 3)))

cowplot::plot_grid(irs_killing_plot, IRS_max_red_boxplot, IRS_max_red_season_scatterplot, nrow = 1,
                   align = "h", axis = "tb")

max_red_act <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  mutate(max_inc_red = max(inc_red),
         max_inc_prop_red = max(inc_prop_red)) %>%
  ungroup() %>%
  filter(inc_red == max_inc_red) %>%
  filter(insecticide == "pirim")
max_red_ben <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  mutate(max_inc_red = max(inc_red),
         max_inc_prop_red = max(inc_prop_red)) %>%
  ungroup() %>%
  filter(inc_red == max_inc_red) %>%
  filter(insecticide == "ben")
max_red_sumi <- reduction %>%
  group_by(insecticide, seasonal_profile) %>%
  mutate(max_inc_red = max(inc_red),
         max_inc_prop_red = max(inc_prop_red)) %>%
  ungroup() %>%
  filter(inc_red == max_inc_red) %>%
  filter(insecticide == "clothi")

act_timings <- unique(act_mod2$timing)[max_red_act$timing[c(7, 11)]]
act_disp <- act_mod2 %>%
  filter(seasonal_profile %in% c(7, 11)) %>%
  filter(ifelse(seasonal_profile == 7, timing == act_timings[1], timing == act_timings[2])) %>%
  dplyr::select(seasonal_profile, timing, t, incidence, counterfactual_incidence)
act_disp$insecticide <- "Pirimiphos Methyl"

ben_timings <- unique(ben_mod2$timing)[max_red_ben$timing[c(7, 11)]]
ben_disp <- ben_mod2 %>%
  filter(seasonal_profile %in% c(7, 11)) %>%
  filter(ifelse(seasonal_profile == 7, timing == ben_timings[1], timing == ben_timings[2])) %>%
  dplyr::select(seasonal_profile, timing, t, incidence, counterfactual_incidence)
ben_disp$insecticide <- "Bendiocarb"

sumi_timings <- unique(sumi_mod2$timing)[max_red_sumi$timing[c(7, 11)]]
sumi_disp <- sumi_mod2 %>%
  filter(seasonal_profile %in% c(7, 11)) %>%
  filter(ifelse(seasonal_profile == 7, timing == sumi_timings[1], timing == sumi_timings[2])) %>%
  dplyr::select(seasonal_profile, timing, t, incidence, counterfactual_incidence)
sumi_disp$insecticide <- "Clothiandin"

overall_eg_plot <- rbind(act_disp, ben_disp, sumi_disp)
hline_dat <- data.frame(seasonal_profile = c(7, 11), timing = c(act_timings, ben_timings, sumi_timings),
                        insecticide = c(rep("Pirimiphos Methyl", 2), rep("Bendiocarb", 2), rep("Clothiandin", 2)))

inc_eg_plot <- ggplot() +
  geom_line(data = overall_eg_plot, aes(x = t - 6570, y = 1000*counterfactual_incidence, 
                                        group = seasonal_profile), col = "black") +
  geom_vline(data = hline_dat, aes(xintercept = timing, col = insecticide, group = seasonal_profile), linetype = "dashed") +
  facet_wrap(~seasonal_profile, nrow = 1) +
  geom_line(data = overall_eg_plot, aes(x = t - 6570, y = 1000*incidence, col = insecticide)) +
  geom_line(data = overall_eg_plot, aes(x = t - 6570, y = 1000*counterfactual_incidence, 
                                        group = seasonal_profile), col = "black") +
  geom_vline(xintercept = c(0, 365, 730), linetype = "dashed") +
  scale_colour_manual(values = c("#FA9F42", "#0B6E4F", "#B2B2B3"),
                      labels = c("Bendiocarb", "Clothiandin", "Pirimiphos Methyl")) +
  labs(x = "Time Since IRS Spraying (Days)", y = "Malaria Incidence Per 1000") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank())

x <- cowplot::plot_grid(irs_killing_plot, IRS_max_red_boxplot, 
                        IRS_max_red_season_scatterplot, nrow = 1,
                   align = "h", axis = "tb", labels = c("A", "B", "C"),
                   label_size = 18)

y <- cowplot::plot_grid(x, inc_eg_plot, nrow = 2, rel_heights = c(1.2, 1),
                   labels = c("", "D"))

ggsave(y, file = "figures/Fig5_Overall.pdf", width = 12, height = 7)
# for (i in 1:length(steph_seasonality_list)) {
#   
#   timing_index_act <- unique(act_mod2$timing)[max_red_act$timing[i]]
#   temp_act <- act_mod2[act_mod2$seasonal_profile == max_red_act$seasonal_profile[i] &
#                          act_mod2$timing == timing_index_act, ]
#   timing_index_ben <- unique(ben_mod2$timing)[max_red_ben$timing[i]]
#   temp_ben <- ben_mod2[ben_mod2$seasonal_profile == max_red_ben$seasonal_profile[i] &
#                          act_mod2$timing == timing_index_ben, ]
#   
#   plot(temp_act$counterfactual_incidence, type = "l",
#        main = paste0("Seasonal Profile ", i, " - Seasonality = ", round(seasonality[i], 2), 
#                      "\n Timing = ", timing_index), ylab = "", xlab = "",
#        ylim = c(0, max(temp_act$counterfactual_incidence)))
#   abline(v = c(0, 365, 730), lty = 2)
#   abline(v = timing_index_act, lty = 1, col = "red")
#   abline(v = timing_index_ben, lty = 1, col = "blue")
#   lines(temp_act$incidence, col = "red")  
#   lines(temp_ben$incidence, col = "blue")  
#   lines(temp_act$counterfactual_incidence)  
#   browser()
# }
