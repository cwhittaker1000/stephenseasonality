# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(ICDMM); #devtools::install_github("https://github.com/jhellewell14/ICDMM", dependencies = TRUE)
library(scales); library(patchwork); library(here); library(zoo)

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
source(here("functions", "time_series_characterisation_functions.R"))
overall <- readRDS(file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)
features_df <- readRDS(file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
cluster_output <- readRDS(file = here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_output$cluster

# Functions
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

summary_function <- function(x) {
  temp <- x %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(prev_mean = mean(prev),
                     prev_lower = quantile(prev, 0.05),
                     prev_upper = quantile(prev, 0.95),
                     inc_mean = mean(Incidence),
                     inc_lower = quantile(Incidence, 0.05),
                     inc_upper = quantile(Incidence, 0.95))
    # dplyr::summarise(prev_mean = mean(prev),
    #                  prev_lower = min(prev),
    #                  prev_upper = max(prev),
    #                  inc_mean = mean(Incidence),
    #                  inc_lower = min(Incidence),
    #                  inc_upper = max(Incidence))
  temp$prev_lower[temp$prev_lower < 0] <- 0
  temp$inc_lower[temp$inc_lower < 0] <- 0
  return(temp)
}

inset_theme <- theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 7))

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
}
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Generating the vector of how mosquito density changes over time 
density_start <- 0.1
density_end <- 10
values <- sigmoid(seq(-10, 10, length.out = 365 * 4)) # How long you want introduction to last - here, 10 years
density_vec <- c(rep(density_start, 365 * 1), 
                 pmin((values * (density_end - density_start)) + density_start, density_end),
                 rep(density_end, 365 * 23))
dens_vec_df <- data.frame(vector_density = density_vec, time = seq(1:length(density_vec))/365)

# Generating seasonal profiles based on actual stephensi data
scaled_year_output <- t(apply(normalised_output, 1, function(x) {
  temp <- approx(x, n = 365)$y
  temp <- 365 * temp/sum(temp)
  temp[temp < 0.1] <- 0.1
  temp }))

# Running the deterministic malaria model with increasing vector density over time and with/without seasonality
steph_seasonality_list <- vector(mode = "list", length = dim(scaled_year_output)[1])
for (i in 1:(dim(scaled_year_output)[1])) {
  steph_seasonality_list[[i]] <- scaled_year_output[i, ]
}

fresh_run <- FALSE
if (fresh_run) {
  multi_outputs <- sapply(1:(length(steph_seasonality_list)), function(x){
    
    set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                            #Model parameters to work out transmission
                                            init_EIR = 0.000001, # initial EIR from which the endemic equilibria solution is created
                                            #These are the mosquito parameters (bionomics)
                                            Q0 = stephensi_data$Q0, 
                                            chi = stephensi_data$chi, 
                                            bites_Bed = stephensi_data$bites_Bed,
                                            bites_Indoors = stephensi_data$bites_Indoors, 
                                            #These are our seasonal variations
                                            scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                            #custom_seasonality = if(x == 1) loc_seasonality else NA,
                                            custom_seasonality = steph_seasonality_list[[x]], # NA for perennial
                                            #This sets up how long we want to run for and the density vec
                                            time_length = length(density_vec),
                                            density_vec = density_vec)
    
    #Process model formulation
    set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
    mod_run <- set_up_model$run(t = 1:length(density_vec))
    out <- set_up_model$transform_variables(mod_run)
    model_ran <- as.data.frame(out)
    model_ran <- data.frame(id = x, model_ran)

  }, simplify = FALSE)
  
  saveRDS(multi_outputs, file = "outputs/malaria_model_running_stephensi_profiles.rds")
} else {
  multi_outputs <- readRDS("outputs/malaria_model_running_stephensi_profiles.rds")
}

# Loading In Saved Runs and Combining Into Overall Data Frame
for (i in 1:(length(steph_seasonality_list))) {
  if (i == 1) {
    temp <- multi_outputs[[i]][, c("id", "t", "prev", "Incidence", "mv")]
  } else {
    temp <- rbind(temp, multi_outputs[[i]][, c("id", "t", "prev", "Incidence", "mv")])
  }
}
temp$id <- as.factor(temp$id)
#rm(multi_outputs)

# Calculating Aggregate Quantities E.g. Degree of Seasonality, Time to 2% Etc
time_to_2_percent <- c()
time_to_2_percent_yearly_average <- c()
seasonality <- c()
for (i in 1:(length(steph_seasonality_list))) {
  temp_calc <- temp$prev[temp$id == i]
  above_2_percent <- which(temp_calc > 0.02)[1]/365
  yearly_averages <- colMeans(matrix(temp_calc, nrow = 365))
  yearly_average_above_2_percent <- which(yearly_averages > 0.02)[1]
  
  time_to_2_percent <- c(time_to_2_percent, above_2_percent)
  time_to_2_percent_yearly_average <- c(time_to_2_percent_yearly_average, yearly_average_above_2_percent)
  
  seasonality <- c(seasonality, calc_incidence_seasonality(steph_seasonality_list[[i]], 3))
}


####### CLUSTER FIGURE EXAMPLE ####### 
one_results <- temp[temp$id %in% which(cluster_membership == 1), ]
one_summary <- summary_function(one_results)
one_summary$id <- 1
one_vectors <- matrix(unlist(steph_seasonality_list[which(cluster_membership == 1)]), nrow = 365)
mean_one_vec <- data.frame(t = 1:365, density = apply(one_vectors, 1, mean)/sum(apply(one_vectors, 1, mean)))
mean_one_vec$lower <- apply(one_vectors, 1, quantile, 0.05)/sum(apply(one_vectors, 1, mean))
mean_one_vec$upper <- apply(one_vectors, 1, quantile, 0.95)/sum(apply(one_vectors, 1, mean))

two_results <- temp[temp$id %in% which(cluster_membership == 2), ]
two_summary <- summary_function(two_results)
two_summary$id <- 2
two_vectors <- matrix(unlist(steph_seasonality_list[which(cluster_membership == 2)]), nrow = 365)
mean_two_vec <- data.frame(t = 1:365, density = apply(two_vectors, 1, mean)/sum(apply(one_vectors, 1, mean)))
mean_two_vec$lower <- apply(two_vectors, 1, quantile, 0.05)/sum(apply(two_vectors, 1, mean))
mean_two_vec$upper <- apply(two_vectors, 1, quantile, 0.95)/sum(apply(two_vectors, 1, mean))

three_results <- temp[temp$id %in%  which(cluster_membership == 3), ]
three_summary <- summary_function(three_results)
three_summary$id <- 3
three_vectors <- matrix(unlist(steph_seasonality_list[which(cluster_membership == 3)]), nrow = 365)
mean_three_vec <- data.frame(t = 1:365, density = apply(three_vectors, 1, mean)/sum(apply(one_vectors, 1, mean)))
mean_three_vec$lower <- apply(three_vectors, 1, quantile, 0.05)/sum(apply(three_vectors, 1, mean))
mean_three_vec$upper <- apply(three_vectors, 1, quantile, 0.95)/sum(apply(three_vectors, 1, mean))

four_results <- temp[temp$id %in% which(cluster_membership == 4), ]
four_summary <- summary_function(four_results)
four_summary$id <- 4
four_vectors <- matrix(unlist(steph_seasonality_list[which(cluster_membership == 4)]), nrow = 365)
mean_four_vec <- data.frame(t = 1:365, density = apply(four_vectors, 1, mean)/sum(apply(one_vectors, 1, mean)))
mean_four_vec$lower <- apply(four_vectors, 1, quantile, 0.05)/sum(apply(four_vectors, 1, mean))
mean_four_vec$upper <- apply(four_vectors, 1, quantile, 0.95)/sum(apply(four_vectors, 1, mean))

mean_cluster_seasonalities <- c(round(mean(seasonality[which(cluster_membership == 1)]), 2), 
                                round(mean(seasonality[which(cluster_membership == 2)]), 2), 
                                round(mean(seasonality[which(cluster_membership == 3)]), 2), 
                                round(mean(seasonality[which(cluster_membership == 4)]), 2))
cols <- c("#E0521A", "#3F88C5", "#44BBA4", "#393E41")
cluster_summary <- rbind(one_summary, two_summary, three_summary, four_summary)
malaria_plots <- ggplot(data = cluster_summary, aes(fill = factor(id))) +
  geom_ribbon(aes(x = t, ymin = inc_lower, ymax = inc_upper), 
              alpha = 0.2) +
  facet_wrap(~id, nrow = 4) +
  lims(x = c(0, 10000)) +
  scale_x_continuous(limits = c(0,10000), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.x = element_line(colour="black", size=0.5),
        axis.line.x.bottom = element_line(colour="black", size=0.5),
        axis.line.y.left = element_line(colour="black", size=0.5),
        panel.spacing = unit(1, "lines")) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, colour = "dark grey")

test1 <- ggplot(mean_one_vec, aes(x = t, y = density)) +
  geom_line(col = "#E0521A", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#E0521A", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_one_vec$upper))) +
  inset_theme
test2 <- ggplot(mean_two_vec, aes(x = t, y = density)) +
  geom_line(col = "#3F88C5", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#3F88C5", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_two_vec$upper))) +
  inset_theme
test3 <- ggplot(mean_three_vec, aes(x = t, y = density)) +
  geom_line(col = "#44BBA4", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#44BBA4", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_three_vec$upper))) + 
  inset_theme
test4 <- ggplot(mean_four_vec, aes(x = t, y = density)) +
  geom_line(col = "#393E41", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#393E41", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_four_vec$upper))) +
  inset_theme

cluster_seas_plotdf <- data.frame(seasonality = seasonality, 
                                  rel_time = time_to_2_percent/min(time_to_2_percent),
                                  id = cluster_membership)
cluster_seas_meandf <- data.frame(id = c(1, 2, 3, 4),
                                  seasonality = c(mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == 1]),
                                                  mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == 2]),
                                                  mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == 3]),
                                                  mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == 4])))
cluster_malaria_plots <- malaria_plots + 
  geom_label(data = cluster_seas_meandf, x = 8200, y = 0.0015, 
             label = paste0("Mean Seasonality = ", mean_cluster_seasonalities),
             fill = "white", label.size = NA) +
  inset_element(test1, 0.05, 0.84, 0.22, 0.99) +
  inset_element(test2, 0.05, 0.585, 0.22, 0.735) +
  inset_element(test3, 0.05, 0.335, 0.22, 0.485) +
  inset_element(test4, 0.05, 0.08, 0.22, 0.23) 

## Cluster Seasonality vs Time to 2%
seasonality_timing <- ggplot(cluster_seas_plotdf, aes(x = seasonality, y = rel_time)) +
  geom_smooth(col = "black") +
  geom_point(size = 2, aes(col = factor(id))) +
  scale_y_continuous(position = "right") +
  scale_colour_manual(values = cols) +
  labs(x = "% Annual Catch In 3 Months", y = "Relative Time to 2%") +
  theme_bw() +
  theme(legend.position = "none")

## Peak Seasonal Incidence
max_incidence <- temp %>%
  group_by(id) %>%
  dplyr::summarise(max_inc = max(Incidence))
mean_max_incidence <- data.frame(id = factor(c(1, 2, 3, 4)),
                                 max_inc = c(mean(max_incidence$max_inc[which(cluster_membership == 1)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 2)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 3)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 4)])),
                                 low_inc = c(min(max_incidence$max_inc[which(cluster_membership == 1)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 2)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 3)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 4)])),
                                 high_inc = c(max(max_incidence$max_inc[which(cluster_membership == 1)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 2)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 3)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 4)])))
seasonality_max_incidence <- ggplot(mean_max_incidence, aes(x = id, y = max_inc, fill = id)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low_inc, ymax = high_inc), width = 0.5) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Maximum Annual Incidence") +
  theme_bw() +
  theme(legend.position = "none") 

rhs <- plot_grid(seasonality_timing, seasonality_max_incidence, nrow = 2, align = "v", axis = "b")
overall <- plot_grid(cluster_malaria_plots, rhs, rel_widths = c(2, 1))


####### URBAN RURAL FIGURE EXAMPLE ####### 
urban_results <- temp[temp$id %in% which(urban_rural == "Urban"), ]
urban_summary <- summary_function(urban_results)
urban_summary$id <- "Urban"
urban_vectors <- matrix(unlist(steph_seasonality_list[which(urban_rural == "Urban")]), nrow = 365)
mean_urban_vec <- data.frame(t = 1:365, density = apply(urban_vectors, 1, mean)/sum(apply(urban_vectors, 1, mean)))
mean_urban_vec$lower <- apply(urban_vectors, 1, quantile, 0.05)/sum(apply(urban_vectors, 1, mean))
mean_urban_vec$upper <- apply(urban_vectors, 1, quantile, 0.95)/sum(apply(urban_vectors, 1, mean))

rural1_results <- temp[temp$id %in% which(urban_rural == "Rural" & features_df$peaks == 1), ]
rural1_summary <- summary_function(rural1_results)
rural1_summary$id <- "Rural One Peak"
rural1_vectors <- matrix(unlist(steph_seasonality_list[which(urban_rural == "Rural" & features_df$peaks == 1)]), nrow = 365)
mean_rural1_vec <- data.frame(t = 1:365, density = apply(rural1_vectors, 1, mean)/sum(apply(rural1_vectors, 1, mean)))
mean_rural1_vec$lower <- apply(rural1_vectors, 1, quantile, 0.05)/sum(apply(rural1_vectors, 1, mean))
mean_rural1_vec$upper <- apply(rural1_vectors, 1, quantile, 0.95)/sum(apply(rural1_vectors, 1, mean))

rural2_results <- temp[temp$id %in% which(urban_rural == "Rural" & features_df$peaks == 2), ]
rural2_summary <- summary_function(rural2_results)
rural2_summary$id <- "Rural Two Peak"
rural2_vectors <- matrix(unlist(steph_seasonality_list[which(urban_rural == "Rural" & features_df$peaks == 2)]), nrow = 365)
mean_rural2_vec <- data.frame(t = 1:365, density = apply(rural2_vectors, 1, mean)/sum(apply(rural2_vectors, 1, mean)))
mean_rural2_vec$lower <- apply(rural2_vectors, 1, quantile, 0.05)/sum(apply(rural2_vectors, 1, mean))
mean_rural2_vec$upper <- apply(rural2_vectors, 1, quantile, 0.95)/sum(apply(rural2_vectors, 1, mean))

mean_cluster_seasonalities <- c(round(mean(seasonality[which(urban_rural == "Urban")]), 2), 
                                round(mean(seasonality[which(urban_rural == "Rural" & features_df$peaks == 1)]), 2), 
                                round(mean(seasonality[which(urban_rural == "Rural" & features_df$peaks == 2)]), 2))
cols <- c("#447604", "#6EA65D","#807A85")
rur_urb_summary <- rbind(urban_summary, rural1_summary, rural2_summary)
malaria_plots <- ggplot(data = rur_urb_summary, aes(fill = factor(id))) +
  geom_ribbon(aes(x = t, ymin = inc_lower, ymax = inc_upper), 
              alpha = 0.2) +
  facet_wrap(~id, nrow = 3) +
  lims(x = c(0, 10000)) +
  scale_x_continuous(limits = c(0,10000), expand = c(0, 0)) +
  scale_fill_manual(values = cols) +
  labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
  theme_bw() + 
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.x = element_line(colour="black", size=0.5),
        axis.line.x.bottom = element_line(colour="black", size=0.5),
        axis.line.y.left = element_line(colour="black", size=0.5),
        panel.spacing = unit(1, "lines")) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, colour = "dark grey")

test1 <- ggplot(mean_urban_vec, aes(x = t, y = density)) +
  geom_line(col = "#807A85", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#807A85", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_urban_vec$upper))) +
  inset_theme
test2 <- ggplot(mean_rural1_vec, aes(x = t, y = density)) +
  geom_line(col = "#447604", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#447604", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_rural1_vec$upper))) +
  inset_theme
test3 <- ggplot(mean_rural2_vec, aes(x = t, y = density)) +
  geom_line(col = "#6EA65D", size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#6EA65D", alpha = 0.2) +
  labs(x = "Time (Days)", y = "Vector Density") +
  lims(y = c(0, max(mean_rural2_vec$upper))) + 
  inset_theme

cluster_seas_plotdf <- data.frame(seasonality = seasonality, 
                                  rel_time = time_to_2_percent/min(time_to_2_percent),
                                  id = overall$city)
cluster_seas_meandf <- data.frame(id = c("Urban", "Rural One Peak", "Rural Two Peak"),
                                  seasonality = c(mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == "Urban"]),
                                                  mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == "Rural" & features_df$peaks == 1]),
                                                  mean(cluster_seas_plotdf$seasonality[cluster_seas_plotdf$id == "Rural" & features_df$peaks == 2])))
cluster_malaria_plots <- malaria_plots + 
  geom_label(data = cluster_seas_meandf, x = 8200, y = 0.0015, 
             label = paste0("Mean Seasonality = ", mean_cluster_seasonalities),
             fill = "white", label.size = NA) +
  inset_element(test2, 0.05, 0.82, 0.22, 0.97) +
  inset_element(test3, 0.05, 0.48, 0.22, 0.63) +
  inset_element(test1, 0.05, 0.13, 0.22, 0.28) 

# add in rural 1 and rural 2
seasonality_timing <- ggplot(cluster_seas_plotdf[cluster_seas_plotdf$id != "Mixture/Unclear", ], aes(x = seasonality, y = rel_time)) +
  geom_smooth(col = "black") +
  geom_point(size = 2, aes(col = factor(id))) +
  scale_y_continuous(position = "right") +
  scale_colour_manual(values = cols) +
  labs(x = "% Annual Catch In 3 Months", y = "Relative Time to 2%") +
  theme_bw() +
  theme(legend.position = "none")

max_incidence <- temp %>%
  group_by(id) %>%
  dplyr::summarise(max_inc = max(Incidence))
mean_max_incidence <- data.frame(id = factor(c(1, 2, 3, 4)),
                                 max_inc = c(mean(max_incidence$max_inc[which(cluster_membership == 1)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 2)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 3)]),
                                             mean(max_incidence$max_inc[which(cluster_membership == 4)])),
                                 low_inc = c(min(max_incidence$max_inc[which(cluster_membership == 1)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 2)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 3)]),
                                             min(max_incidence$max_inc[which(cluster_membership == 4)])),
                                 high_inc = c(max(max_incidence$max_inc[which(cluster_membership == 1)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 2)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 3)]),
                                              max(max_incidence$max_inc[which(cluster_membership == 4)])))
seasonality_max_incidence <- ggplot(mean_max_incidence, aes(x = id, y = max_inc, fill = id)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = low_inc, ymax = high_inc), width = 0.5) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")) +
  scale_y_continuous(position = "right") +
  labs(x = "", y = "Maximum Annual Incidence") +
  theme_bw() +
  theme(legend.position = "none") 

########################

# Subsets of the Data According to Unimodal/Bimodal, Rural/Urban & Cluster Membership

## Unimodal/Bimodal
unimodal <- which(features_df$peaks == 1)
bimodal <- which(features_df$peaks == 2)
unimodal_most <- unimodal[order(seasonality[unimodal])[26:50]]
unimodal_least <- unimodal[order(seasonality[unimodal])[1:25]]

## Urban/Rural
urban <- which(features_df$peaks == 1 & features_df$cit == "Urban")
rural_one <- which(features_df$peaks == 1 & features_df$cit == "Rural")
rural_two <- which(features_df$peaks == 2 & features_df$cit == "Rural")

## Cluster Membership
one <- which(cluster_membership == 1)
two <- which(cluster_membership == 2)
three <- which(cluster_membership == 3)
four <- which(cluster_membership == 4)

## Mean Seasonality
mean(seasonality[unimodal_most])
mean(seasonality[unimodal_least])
mean(seasonality[bimodal])

mean(seasonality[urban])
mean(seasonality[rural_one])
mean(seasonality[rural_two])

mean(seasonality[one])
mean(seasonality[two])
mean(seasonality[three])
mean(seasonality[four])

## Mean Time to 2 Percent
mean(time_to_2_percent[unimodal_most])
mean(time_to_2_percent[unimodal_least])
mean(time_to_2_percent[bimodal])

mean(time_to_2_percent[urban])
mean(time_to_2_percent[rural_one])
mean(time_to_2_percent[rural_two])

mean(time_to_2_percent[one])
mean(time_to_2_percent[two])
mean(time_to_2_percent[three])
mean(time_to_2_percent[four])

plot(seasonality, time_to_2_percent/min(time_to_2_percent), pch = 20, xlab = "% Annual Catch In 3 Months", ylab = "Relative Time to 2%")

# Plotting Seasonal Vector Profiles
unimodal_most_vectors <- steph_seasonality_list[unimodal_most]

## Overall Plotting
overall_summary <- summary_function(temp)
# overall_summary <- temp %>%
#   dplyr::group_by(t) %>%
#   dplyr::summarise(prev_mean = mean(prev),
#                    prev_lower = quantile(prev, 0.01),
#                    prev_upper = quantile(prev, 0.99),
#                    inc_mean = mean(Incidence),
#                    inc_lower = quantile(Incidence, 0.01),
#                    inc_upper = quantile(Incidence, 0.99))
# overall_summary$prev_lower[temp$prev_lower < 0] <- 0
# overall_summary$inc_lower[temp$inc_lower < 0] <- 0

ggplot() +
  geom_ribbon(data = overall_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2)
ggplot() +
  geom_ribbon(data = overall_summary, aes(x = t, ymin = inc_lower, ymax = inc_upper), 
              alpha = 0.2)

## Bimodal vs Unimodal Plotting
unimodal_most_results <- temp[temp$id %in% unimodal_most, ]
unimodal_most_summary <- summary_function(unimodal_most_results)

unimodal_least_results <- temp[temp$id %in% unimodal_least, ]
unimodal_least_summary <- summary_function(unimodal_least_results)

bimodal_results <- temp[temp$id %in% bimodal, ]
bimodal_summary <- summary_function(bimodal_results)

ggplot() +
  geom_ribbon(data = bimodal_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2) +
  geom_ribbon(data = unimodal_least_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2, fill = "blue") +
  geom_ribbon(data = unimodal_most_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2, fill = "red") +
  labs(x = "Time (Days)", y = "Prevalence (%)")

## Rural/Urban Plotting 
urban_results <- temp[temp$id %in% urban, ]
urban_summary <- summary_function(urban_results)

rural_one_results <- temp[temp$id %in% rural_one, ]
rural_one_summary <- summary_function(rural_one_results)

rural_two_results <- temp[temp$id %in% rural_two, ]
rural_two_summary <- summary_function(rural_two_results)

ggplot() +
  geom_ribbon(data = urban_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2) +
  geom_ribbon(data = rural_one_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2, fill = "blue") +
  geom_ribbon(data = rural_two_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
              alpha = 0.2, fill = "red")

## Cluster Plotting

#############################################################



unimodal_seasonality <- seasonality[unimodal]
unimodal_most <- unimodal[order(unimodal_seasonality)[26:50]]
unimodal_least <- unimodal[order(unimodal_seasonality)[1:25]]
bimodal <- which(features_df$peaks == 2)


urban_results <- temp[temp$id %in% urban, ]
urban_summary <- urban_results %>%
  dplyr::group_by(t) %>%
  dplyr::summarise(prev_mean = mean(prev),
                   prev_lower = prev_mean - 1.96 * sd(prev),
                   prev_upper = prev_mean + 1.96 * sd(prev),
                   inc_mean = mean(Incidence),
                   inc_lower = inc_mean - 1.96 * sd(Incidence),
                   inc_upper = inc_mean + 1.96 * sd(Incidence))
urban_summary$prev_lower[urban_summary$prev_lower < 0] <- 0
urban_summary$inc_lower[urban_summary$inc_lower < 0] <- 0

rural_one_results <- temp[temp$id %in% rural_one, ]
rural_one_summary <- rural_one_results %>%
  dplyr::group_by(t) %>%
  dplyr::summarise(prev_mean = mean(prev),
                   prev_lower = prev_mean - 1.96 * sd(prev),
                   prev_upper = prev_mean + 1.96 * sd(prev),
                   inc_mean = mean(Incidence),
                   inc_lower = inc_mean - 1.96 * sd(Incidence),
                   inc_upper = inc_mean + 1.96 * sd(Incidence))
rural_one_summary$prev_lower[rural_one_summary$prev_lower < 0] <- 0
rural_one_summary$inc_lower[rural_one_summary$inc_lower < 0] <- 0

rural_two_results <- temp[temp$id %in% rural_two, ]
rural_two_summary <- rural_two_results %>%
  dplyr::group_by(t) %>%
  dplyr::summarise(prev_mean = mean(prev),
                   prev_lower = prev_mean - 1.96 * sd(prev),
                   prev_upper = prev_mean + 1.96 * sd(prev),
                   inc_mean = mean(Incidence),
                   inc_lower = inc_mean - 1.96 * sd(Incidence),
                   inc_upper = inc_mean + 1.96 * sd(Incidence))
rural_two_summary$prev_lower[rural_two_summary$prev_lower < 0] <- 0
rural_two_summary$inc_lower[rural_two_summary$inc_lower < 0] <- 0

ggplot() +
  geom_ribbon(data = urban_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), alpha = 0.2) +
  geom_ribbon(data = rural_one_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), alpha = 0.2, fill = "blue") +
  geom_ribbon(data = rural_two_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), alpha = 0.2, fill = "red")




  lims(x = c(0, 10000)) +
  scale_x_continuous(limits = c(0,10000), expand = c(0, 0)) +
  scale_color_manual(values = cols) +
  facet_wrap(~factor(id), nrow = 4) +
  labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
  theme_bw() + 
  scale_y_continuous(limits = c(-1,14), expand = c(0, 0), position = "left",
                     breaks = c(0, 5, 10)) +
  theme(legend.position = "none", 
        #panel.grid.major.x = element_line(colour="grey", size=0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.x = element_line(colour="black", size=0.5),
        axis.line.x.bottom = element_line(colour="black", size=0.5),
        axis.line.y.left = element_line(colour="black", size=0.5),
        panel.spacing = unit(1, "lines")) +
  #axis.line.y.right = element_line(colour="dark grey", size=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, colour = "dark grey")


# Plotting the Individual Malaria Profiles
ex_ind <- c(1, 7, 14, 21)
example_malaria <- ggplot(temp[temp$id %in% ex_ind, ], aes(x = t, y = 10000 * Incidence, col = factor(id))) +
  geom_line(size = 1.5) +
  lims(x = c(0, 10000)) +
  scale_x_continuous(limits = c(0,10000), expand = c(0, 0)) +
  scale_color_manual(values = cols) +
  facet_wrap(~factor(id), nrow = 4) +
  labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
  theme_bw() + 
  scale_y_continuous(limits = c(-1,14), expand = c(0, 0), position = "left",
                     breaks = c(0, 5, 10)) +
  theme(legend.position = "none", 
        #panel.grid.major.x = element_line(colour="grey", size=0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.x = element_line(colour="black", size=0.5),
        axis.line.x.bottom = element_line(colour="black", size=0.5),
        axis.line.y.left = element_line(colour="black", size=0.5),
        panel.spacing = unit(1, "lines")) +
  #axis.line.y.right = element_line(colour="dark grey", size=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, colour = "dark grey")

test1 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[1]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[1]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test2 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[2]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[2]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test3 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[3]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[3]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test4 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == 0, ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[4]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
malaria_plots <- example_malaria + inset_element(test1, 0.05, 0.82, 0.22, 0.97) +
  inset_element(test2, 0.05, 0.565, 0.22, 0.715) +
  inset_element(test3, 0.05, 0.315, 0.22, 0.465) +
  inset_element(test4, 0.05, 0.06, 0.22, 0.21)

# Plotting Seasonality vs Time to 2%
seasonality_dynamics_df <- data.frame(mal_seas = seasonality, rel_time = time_to_2_percent/min(time_to_2_percent))
seasonal_establishment_plot <- ggplot(seasonality_dynamics_df, aes(x = mal_seas, y = rel_time, col = mal_seas)) +
  scale_colour_gradientn(colours=cols) +
  geom_line(size = 2) +
  theme_bw() +
  labs(x = "% Malaria Incidence In Any 3 Month period",
       y = "Relative Increase In Time Taken\nTo Reach 2% Prevalence") +
  theme(legend.position = "none")

# Overall Plotting
left <- plot_grid(vector_plot, seasonal_establishment_plot, nrow = 2, align = "v", axis = "lrtl")
overall <- plot_grid(left, malaria_plots, rel_widths = c(1, 2))
#7.5 x 15

# # Plotting Seasonal Vector Profiles
# cols <- c("#A0CFD3", "#8D94BA", "#9A7AA0", "#87677B")
# ex_ind2 <- c(1, 3, 7, 11, 14, 17)
# cols_21 <- colorRampPalette(c(cols))
# seasonal_vector_plot <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar %in% c(scalar_values[ex_ind2], 0), ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
#   geom_line(size = 2) +
#   scale_colour_manual(values = rev(cols_21(length(ex_ind2) + 1))) +
#   theme_bw() +
#   labs(x = "Time (Days)", y = "Vector Density") +
#   theme(legend.position = "none")
# annual_vector_plot <- ggplot(dens_vec_df, aes(x = time, y = vector_density)) +
#   geom_line(size = 1) +
#   theme_bw() +
#   labs(x = "Time (Years)", y = "Annual Vector\nDensity") +
#   geom_line() +
#   lims(x = c(0, 6)) +
#   theme(axis.title = element_text(size = 10),
#         panel.background = element_blank())
# vector_plot <- seasonal_vector_plot + inset_element(annual_vector_plot, 0.01, 0.52, 0.50, 0.97)









# (optional) Specifies whether graphs in the grid should be horizontally ("h") 
# or vertically ("v") aligned. Options are "none" (default), "hv" (align in both
# directions), "h", and "v".

# (optional) Specifies whether graphs should be aligned by the left ("l"), 
# right ("r"), top ("t"), or bottom ("b") margins. Options are "none" 
# (default), or a string of any combination of l, r, t, and b in any order 
# (e.g. "tblr" or "rlbt" for aligning all margins). Must be specified if 
# any of the graphs are complex (e.g. faceted) and alignment is specified 
# and desired. See align_plots() for details.

# par(mfrow = c(1, 1))
# plot(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 2), ylab = "Seasonality factor", xlab = "Time (days)", bty = "l", type = "l", col = "red", lwd = 2)
# lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 1), col = "blue", lwd = 2) # scalar controls how flat
# lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 0.7), col = "black", lwd = 2)

# plot(density_vec, bty = "l", ylab = "Vector density", xlab = "Time (days)")
# abline(v = seq(365, length(density_vec), by = 365), lty = 2)

# a <- ggplot(dens_vec_df, aes(x = time, y = vector_density)) +
#   geom_line() +
#   theme_bw() +
#   geom_rect(aes(xmin = 5, xmax = 6, ymin = 9, ymax = 11), 
#             fill = "grey", alpha = 0.03) +
#   labs(x = "Time (Years)", y = "Average Annual Vector Density") +
#   geom_line() +
#   lims(x = c(0, 6))

# d <- ggplot(temp, aes(x = t, y = 10000 * Incidence, col = factor(id))) +
#   geom_line() +
#   lims(x = c(0, 11000)) +
#   scale_color_manual(values = cols) +
#   theme_bw() + 
#   labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
#   theme(legend.position = "right")
# d

# palette(rev(c("#C0E1CB", "#8FD5A6", "#61BA81", "#329F5B", "#1F9151", "#0C8346", "#0D5D56")))
# palette(gray(seq(0,.9,len = 25)))
# par(mfrow = c(3, 1))
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$prev, type = "l", bty = "l", ylab = "Prevalence", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0,0.055), xlim = c(0, 10000))
#   } else {
#     lines(multi_outputs[[i]]$prev, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$Incidence, type = "l", bty = "l", ylab = "Vector Density", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0, max(multi_outputs[[i]]$Incidence)), xlim = c(0, 7500))
#   } else {
#     lines(multi_outputs[[i]]$Incidence, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$mv, type = "l", bty = "l", ylab = "Vector Density", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0, max(multi_outputs[[i]]$mv)), xlim = c(0, 7500))
#   } else {
#     lines(multi_outputs[[i]]$mv, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)), type = "l", bty = "l", ylab = "Mean Annual Vector Density", xlab = "Time (Years)", 
#          col = palette()[i], ylim = c(0, max(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)))))
#   } else {
#     lines(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)), col = palette()[i])
#   }
# }
# 
# 
# ggplot() +
#   geom_ribbon(data = one_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
#               alpha = 0.2, fill = "#E0521A") +
#   geom_ribbon(data = two_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
#               alpha = 0.2, fill = "#3F88C5") +
#   geom_ribbon(data = three_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
#               alpha = 0.2, fill = "#44BBA4") +
#   geom_ribbon(data = three_summary, aes(x = t, ymin = prev_lower, ymax = prev_upper), 
#               alpha = 0.2, fill = "#393E41")