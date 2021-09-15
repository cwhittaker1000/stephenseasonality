#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(sf); library(ggplot2);
library(zoo); library(cowplot)

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))

# Load metadata and admin unit geometries
metadata <- readRDS(here("data", "processed", "metadata_and_processed_counts.rds"))
admin2 <- readRDS(here("data", "processed", "simplified_admin2.rds"))
admin1 <- readRDS(here("data", "processed", "simplified_admin1.rds"))
admin0 <- readRDS(here("data", "processed", "complex_admin0.rds"))

# Helper function to get ggplot2 default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Collating centroids for all of the locations
for (i in 1:nrow(metadata)) {
  
  # Get relevant admin unit and convert to sf object
  if(!is.na(metadata$admin2[i])) {
    metadata_admin <- metadata$admin2[i]
    admin_2_data <- admin2[admin2$NAME_2 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_2_data)
    centroid <- st_centroid(admin_geometry)
    centroid <- centroid[c("continent", "ISO", "NAME_0", "geometry", "NAME_1", "NAME_2")]
  } else {
    metadata_admin <- metadata$admin1[i]
    admin_1_data <- admin1[admin1$NAME_1 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_1_data)
    centroid <- st_centroid(admin_geometry)
    centroid <- centroid[c("continent", "ISO", "NAME_0", "geometry", "NAME_1")]
    centroid$NAME_2 <- NA
  }
  
  if (i == 1) {
    centroid_metadata <- centroid
  } else {
    centroid_metadata <- rbind(centroid_metadata, centroid)
  }
  
  print(i)
  
}

# Plotting the centroids of study locations on a background of the below countries
countries <- c("Afghanistan", "Bangladesh", "Bhutan", "Djibouti", "Egypt", "Ethiopia",
               "India", "Iran", "Iraq", "Israel", "Jordan", "Kenya", "Lebanon", "Myanmar", 
               "Nepal", "Oman", "Pakistan", "Saudi Arabia", "Sudan", "Somalia", "South Sudan", "Yemen")
red_admin0 <- admin0[admin0$NAME_0 %in% countries, ]
red_admin0$stephensi <- c("Yes", "No", "No", "Yes", "No", "No", 
                          "Yes", "Yes", "No", "No", "No", "No", "No", "Yes", 
                          "No", "No", "Yes", "No", "No", "No", "No", "No")
a <- ggplot() + 
  geom_sf(data = red_admin0) +
  geom_sf(data = st_jitter(centroid_metadata, factor = .008), 
          aes(fill = NAME_0), col = "black", size = 2, shape = 21) +
  scale_fill_manual(values = gg_color_hue(6))

ggsave2(file = here("figures", "Fig1A.pdf"), plot = a, width = 11, height = 6, dpi = 500)

# Plotting some representative time-series
retain_index <- metadata$id 
for_plot_index <- c(5, 33, 38, 42, 63, 70)
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
lower_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
upper_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
counter <- 1
for (i in 1:length(for_plot_index)) {
  
  index <- which(retain_index %in% for_plot_index[i])
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", for_plot_index[i], ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", for_plot_index[i], ".rds"))
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
  mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  lower <- as.numeric(exp(f_lower)[order(all_timepoints)])
  upper <- as.numeric(exp(f_upper)[order(all_timepoints)])
  mean_realisation[counter, ] <- mean
  lower_realisation[counter, ] <- lower
  upper_realisation[counter, ] <- upper
  counter <- counter + 1
  # plot(mean_realisation[i, ], main = paste0(metadata$id[index], "  ", metadata$country[index]))
  # browser()
}

# Plotting all of the outputs to see which to feature in Fig 1A plot
mean <- as.data.frame(mean_realisation)
mean$country <- metadata$country[fit_index]
mean_pv <- mean %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "mean")

lower <- as.data.frame(lower_realisation)
lower$country <- metadata$country[fit_index]
lower_pv <- lower %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "lower")

upper <- as.data.frame(upper_realisation)
upper$country <- metadata$country[fit_index]
upper_pv <- upper %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "upper")

catch_data <- metadata %>%
  filter(id %in% for_plot_index) %>%
  select(country, Jan:Dec)
colnames(catch_data) <- c("country", seq(2, 24, length.out = 12)) 
catch_data <- catch_data %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
  left_join(total_raw_catch, by = "country") %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  select(country, timepoint, raw_catch)

overall <- mean_pv %>%
  left_join(lower_pv, by = c("country", "timepoint")) %>%
  left_join(upper_pv, by = c("country", "timepoint")) %>%
  left_join(total_per_country, by = c("country")) %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  left_join(selected, by = c("country", "timepoint")) 

b <- ggplot(data = overall_plus) +
  geom_path(aes(x = timepoint, y = mean, col = country), size = 2) +
  geom_point(aes(x = timepoint, y = raw_catch)) +
  geom_ribbon(aes(x = timepoint, ymin = lower, ymax = upper, fill = country), alpha = 0.2) +
  facet_wrap(~country, scales = "free_y") +
  scale_y_continuous(limits=c(0, NA)) +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(2, 24, length.out = 12))

ggsave2(file = here("figures", "Fig1B.pdf"), plot = b, width = 11, height = 6, dpi = 500)

# Scrap Code
# ggplot() + 
#   geom_sf(data = red_admin0, aes(fill = stephensi)) +
#   scale_fill_manual(values = c("grey60", "white")) +
#   geom_sf(data = st_jitter(centroid_metadata, factor = .005), col = "red", size = 1) #+
#   #coord_sf(xlim = c(20, 105), ylim = c(-5, 40), expand = FALSE)
# 
# 
# plot(st_geometry(admin0))
# plot(st_geometry(centroid_metadata), pch = 3, col = 'red', add = TRUE)
# 
# ggplot() + 
#   geom_sf(data = admin0) +
#   geom_sf(data = st_jitter(centroid_metadata, factor = .005), col = "red", size = 1) +
#   coord_sf(xlim = c(30, 105), ylim = c(-5, 40), expand = FALSE)
# 
# fit_index <- retain_index %in% c(5, 33, 42, 38, 63, 70)
# example_country_fits <- as.data.frame(mean_realisation[fit_index, ])
# example_country_fits$country <- metadata$country[fit_index]
# total_per_country <- example_country_fits %>%
#   pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
#   group_by(country) %>%
#   summarise(total = sum(catch))
# example_country_fits <- example_country_fits %>%
#   pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
#   left_join(total_per_country, by = "country") %>%
#   mutate(norm_catch = catch/total) %>%
#   mutate(timepoint = as.numeric(gsub("V", "", timepoint)))
# 
# ggplot(data = example_country_fits) +
#   geom_path(aes(x = timepoint, y = norm_catch, col = country), size = 2) +
#   facet_wrap(~country, scales = "free_y") +
#   scale_y_continuous(limits=c(0, NA))

# for (i in 1:length(retain_index)) {
#   
#   index <- retain_index[i]
#   
#   # Loading in and processing the fitted time-series
#   if (prior == "informative") {
#     STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
#   } else if (prior == "uninformative") {
#     STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
#   }  
#   
#   # Extracting the mean fitted time-series
#   timepoints <- STAN_output$timepoints
#   all_timepoints <- STAN_output$all_timepoints
#   ordered_timepoints <- all_timepoints[order(all_timepoints)]
#   MCMC_output <- STAN_output$chain_output
#   f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
#   f_mean <- apply(f, 2, mean)
#   mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
#   plot(mean, main = paste0(i, "   ", index, "    ", metadata$id[i], "  ", metadata$country[i]))
#   browser()
# }
# total_raw_catch <- catch_data %>%
#   pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
#   group_by(country) %>%
#   summarise(total_raw = sum(raw_catch))
# total_per_country <- mean %>%
#   pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
#   group_by(country) %>%
#   summarise(total = sum(catch))
