# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(zoo)

# Loading in custom functions
source(here("functions", "time_series_characterisation_functions.R"))

# Loading in ecological covariate data for each country
dji <- read.csv("data/environmental_covariates/DJI_environment_data.csv") %>%
  mutate(country = "aOther")
eth <- read.csv("data/environmental_covariates/ETH_environment_data.csv") %>%
  mutate(country = "aOther")
eri <- read.csv("data/environmental_covariates/ERI_environment_data.csv") %>%
  mutate(country = "aOther")
sdn <- read.csv("data/environmental_covariates/SDN_environment_data.csv") %>%
  mutate(country = "aOther")
som <- read.csv("data/environmental_covariates/SOM_environment_data.csv") %>%
  mutate(country = "aOther")
hoa_covariates_raw <- data.table::rbindlist(list(dji, eth, eri, sdn, som), fill = TRUE)
country_vector <- c(rep("dji", dim(dji)[1]), rep("eth", dim(eth)[1]), rep("eri", dim(eri)[1]), rep("sdn", dim(sdn)[1]), rep("som", dim(som)[1]))

# Loading in geographical data (shapefiles) for each country
dji_shp <- st_as_sf(getData('GADM', country = 'DJI', level = 0, path = here("data/admin_units")))
eth_shp <- st_as_sf(getData('GADM', country = 'ETH', level = 0, path = here("data/admin_units")))
eri_shp <- st_as_sf(getData('GADM', country = 'ERI', level = 0, path = here("data/admin_units")))
sdn_shp <- st_as_sf(getData('GADM', country = 'SDN', level = 0, path = here("data/admin_units")))
som_shp <- st_as_sf(getData('GADM', country = 'SOM', level = 0, path = here("data/admin_units")))
hoa_outline <- rbind(dji_shp, eth_shp, eri_shp, sdn_shp, som_shp)

# Loading in Metadata and Location Data
ts_metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE)))
cluster_membership <- readRDS(here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_membership[, c("id", "cluster")]
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))
overall <- overall %>%
  left_join(cluster_membership, by = "id")

# Loading In Environmental Covariates 
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(cluster, country, population_per_1km:mean_temperature_driest_quarter, -LC_190, -LC_210, -mean_temperature_driest_quarter) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country_peaks = paste0(cluster, "_", country))
data$cluster <- ifelse(data$cluster == 1, "one", "two")
data$cluster <- as.factor(data$cluster)
data$population_per_1km <- log(data$population_per_1km)

# Processing Arran's Horn of Africa Country Rasters
lon <- hoa_covariates_raw$lon
lat <- hoa_covariates_raw$lat
hoa_envt_covariates <- hoa_covariates_raw %>%
  dplyr::select(-lon, -lat, -ISO, -cell, -goat, -sheep, -cattle, -irrigation, -nrow, -ncol, -worldclim_9) %>%
  mutate(country_peaks = paste0(country, "_1")) %>%
  mutate(cluster = 3) %>%
  rename(population_per_1km = population_density) %>%
  mutate(LC_121 = 0, LC_71 = 0) %>%
  mutate(country = as.factor(country), cluster = as.factor(cluster))
hoa_envt_covariates$population_per_1km <- log(hoa_envt_covariates$population_per_1km + 1)

# Checking where the NAs Are for Ecological Covariates and Replacing With 0s (NA = 0 In This Context As Means Country Lacks That Landcover Type)
for (i in 1:ncol(hoa_envt_covariates)) {
  print(c(i, sum(is.na(hoa_envt_covariates[, i]))))
}
colnames(hoa_envt_covariates)[43:52]
hoa_envt_covariates[is.na(hoa_envt_covariates)] <- 0

# Predicting Horn of Africa Using Each of the Random Forests 
iterations_ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
raster_iterations <- vector(mode = "list", length = 25L)
for (i in 1:25) {
  
  # Load In The Random Forest Objects and Check Variables
  final_random_forest_fit_ups <- iterations_ups$model[[i]]
  rf_recipe <- iterations_ups$recipe[[i]]
  vars <- rf_recipe$var_info$variable
  if (!identical(vars[!(vars %in% colnames(hoa_envt_covariates))], character(0))) {
    stop("Not identical, variables are fucked")
  }
  
  # Generating Normalised and Centred Covariates Using Previous Recipe
  juiced_dji <- bake(object = rf_recipe, new_data = hoa_envt_covariates[country_vector == "dji", ])
  juiced_eth <- bake(object = rf_recipe, new_data = hoa_envt_covariates[country_vector == "eth", ])
  juiced_eri <- bake(object = rf_recipe, new_data = hoa_envt_covariates[country_vector == "eri", ])
  juiced_sdn <- bake(object = rf_recipe, new_data = hoa_envt_covariates[country_vector == "sdn", ])
  juiced_som <- bake(object = rf_recipe, new_data = hoa_envt_covariates[country_vector == "som", ])

  # Making Country-Specific Rasters of Model Predictions to Then Stitch Together
  dji_predictions <- predict(final_random_forest_fit_ups, juiced_dji, "prob")
  dji_raster <- raster(ncol = dji$ncol[1], nrow = dji$nrow[1], ext = extent(dji_shp))
  dji_raster[dji$cell] <- dji_predictions$.pred_one
  dji_raster_plot <- as.data.frame(as(dji_raster, "SpatialPixelsDataFrame"))
  dji_raster_plot$lat <- dji$lat
  dji_raster_plot$lon <- dji$lon
  
  eth_predictions <- predict(final_random_forest_fit_ups, juiced_eth, "prob")
  eth_raster <- raster(ncol = eth$ncol[1], nrow = eth$nrow[1], ext = extent(eth_shp))
  eth_raster[eth$cell] <- eth_predictions$.pred_one
  eth_raster_plot <- as.data.frame(as(eth_raster, "SpatialPixelsDataFrame"))
  eth_raster_plot$lat <- eth$lat
  eth_raster_plot$lon <- eth$lon
  
  eri_predictions <- predict(final_random_forest_fit_ups, juiced_eri, "prob")
  eri_raster <- raster(ncol = eri$ncol[1], nrow = eri$nrow[1], ext = extent(eri_shp))
  eri_raster[eri$cell] <- eri_predictions$.pred_one
  eri_raster_plot <- as.data.frame(as(eri_raster, "SpatialPixelsDataFrame"))
  eri_raster_plot$lat <- eri$lat
  eri_raster_plot$lon <- eri$lon
  
  sdn_predictions <- predict(final_random_forest_fit_ups, juiced_sdn, "prob")
  sdn_raster <- raster(ncol = sdn$ncol[1], nrow = sdn$nrow[1], ext = extent(sdn_shp))
  sdn_raster[sdn$cell] <- sdn_predictions$.pred_one
  sdn_raster_plot <- as.data.frame(as(sdn_raster, "SpatialPixelsDataFrame"))
  sdn_raster_plot$lat <- sdn$lat
  sdn_raster_plot$lon <- sdn$lon
  
  som_predictions <- predict(final_random_forest_fit_ups, juiced_som, "prob")
  som_raster <- raster(ncol = som$ncol[1], nrow = som$nrow[1], ext = extent(som_shp))
  som_raster[som$cell] <- som_predictions$.pred_one
  som_raster_plot <- as.data.frame(as(som_raster, "SpatialPixelsDataFrame"))
  som_raster_plot$lat <- som$lat
  som_raster_plot$lon <- som$lon
  
  # Generating Empty Snap Raster To Stitch All The Countries Together
  snap <- raster(resolution = c(0.09748167,0.09746483), xmn = 21, xmx = 60, ymn = -2, ymx = 25)
  snap$layer <- NA
  dji_file <- projectRaster(dji_raster, snap, method = "bilinear")
  eth_file <- projectRaster(eth_raster, snap, method = "bilinear")
  eri_file <- projectRaster(eri_raster, snap, method = "bilinear")
  sdn_file <- projectRaster(sdn_raster, snap, method = "bilinear")
  som_file <- projectRaster(som_raster, snap, method = "bilinear")
  x <- list(dji_file, eth_file, eri_file, sdn_file, som_file)
  x$fun <- mean
  x$na.rm <- TRUE
  x$tolerance <- 1
  combined_raster <- do.call(merge, x)
  
  # Saving the Raster as a Dataframe and Into the Overall list
  raster_plot <- as.data.frame(as(combined_raster, "SpatialPixelsDataFrame"))
  raster_plot$id <- i
  raster_iterations[[i]] <- raster_plot
  
  print(i)
  
}
saveRDS(raster_iterations, file = here("outputs", "seasonality_prediction", "simple_random_forest_seasonality_prediction.rds"))

# Plotting the Output
ken_shp <- st_as_sf(getData('GADM', country = 'KEN', level = 0, path = here("data/admin_units")))
egy_shp <- st_as_sf(getData('GADM', country = 'EGY', level = 0, path = here("data/admin_units")))
uga_shp <- st_as_sf(getData('GADM', country = 'UGA', level = 0, path = here("data/admin_units")))
drc_shp <- st_as_sf(getData('GADM', country = 'COD', level = 0, path = here("data/admin_units")))
ssd_shp <- st_as_sf(getData('GADM', country = 'SSD', level = 0, path = here("data/admin_units")))
car_shp <- st_as_sf(getData('GADM', country = 'CAF', level = 0, path = here("data/admin_units")))
chad_shp <- st_as_sf(getData('GADM', country = 'TCD', level = 0, path = here("data/admin_units")))
lib_shp <- st_as_sf(getData('GADM', country = 'LBY', level = 0, path = here("data/admin_units")))
rwa_shp <- st_as_sf(getData('GADM', country = 'RWA', level = 0, path = here("data/admin_units")))
tan_shp <- st_as_sf(getData('GADM', country = 'TZA', level = 0, path = here("data/admin_units")))
bdi_shp <- st_as_sf(getData('GADM', country = 'BDI', level = 0, path = here("data/admin_units")))
sda_shp <- st_as_sf(getData('GADM', country = 'SAU', level = 0, path = here("data/admin_units")))
yem_shp <- st_as_sf(getData('GADM', country = 'YEM', level = 0, path = here("data/admin_units")))
uae_shp <- st_as_sf(getData('GADM', country = 'ARE', level = 0, path = here("data/admin_units")))
oman_shp <- st_as_sf(getData('GADM', country = 'OMN', level = 0, path = here("data/admin_units")))
qat_shp <- st_as_sf(getData('GADM', country = 'QAT', level = 0, path = here("data/admin_units")))
hoa_neighbours <- rbind(ken_shp, egy_shp, uga_shp, drc_shp, ssd_shp, car_shp, uae_shp, oman_shp,
                        chad_shp, lib_shp, rwa_shp, tan_shp, bdi_shp, sda_shp, yem_shp, qat_shp)

# Average Over the Random Forest Outputs
all_rasters <- bind_rows(raster_iterations) %>% # alt: all_rasters <- do.call("rbind", raster_iterations)
  group_by(x, y) %>%
  summarise(layer = mean(layer, na.rm = TRUE))
fig4a <- ggplot() +
  geom_tile(data = all_rasters, aes(x = x, y = y, fill = layer)) +
  # scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"),
  #                      midpoint = 0.5, limits = c(0, 1), space = "Lab",
  #                      na.value = "grey50", guide = "colourbar") +
  scale_fill_gradient2(low = palette()[1], mid = "white", high = palette()[2],
                       midpoint = 0.5, limits = c(0, 1), space = "Lab",
                       na.value = "grey50", guide = "colourbar") +
  geom_sf(data = hoa_neighbours, col = gray(0.6), fill = gray(0.90), size = 1) +
  geom_sf(data = hoa_outline, col = "black", fill = NA, size = 1) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  lims(x = c(22, 52), y = c(-2, 25)) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold")) +
  guides(fill = guide_colourbar(title = "Seasonality", ticks = FALSE))

# Annual Catch Figures Describing Seasonality 
cluster <- readRDS(here("data", "systematic_review_results", "cluster_membership.rds"))
mean_one <- mean(ts_metadata$per_ind_4_months[cluster$cluster == 1])
mean_two <- mean(ts_metadata$per_ind_4_months[cluster$cluster == 2])
fig4b <- ggplot(ts_metadata, aes(x = 100 * per_ind_4_months, fill = factor(peaks))) +
  geom_histogram(bins = 10, colour = "dark grey", position = "identity", alpha = 0.4) +
  theme_bw() +
  geom_segment(x = 100 * mean_one, xend = 100 * mean_one, y = 0, yend = 17, size = 1, col = palette()[2], linetype = "dashed") +
  geom_segment(x = 100 * mean_two, xend = 100 * mean_two, y = 0, yend = 17, size = 1, col = palette()[1], linetype = "dashed") + 
  scale_fill_manual(values = palette()[2:1]) + 
  scale_x_continuous(breaks = seq(30, 100, 10)) +
  labs(fill = "Cluster", x = "% of Annual Catch In 4 Months", y = "Number of Studies") +
  theme(legend.position = "left")

# Extracting Mean Realisations
id <- overall$id
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(id), ncol = (12 * interpolating_points + 1))
overdisp <- vector(mode = "numeric",length = length(id))
for (i in 1:length(id)) {
  
  index <- id[i]
  
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
  overdisp[i] <- mean(MCMC_output[, "overdispersion"])
}

# Annual Surveillance Considerations - Probability of Detecting Stephensi If Catches Only Done In X Months
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Reordering mean fitted time-series to start at max
reordered_mean_realisation <- matrix(nrow = length(id), ncol = (12 * interpolating_points + 1))
start_index <- apply(normalised_output, 1, function(x) which(x == max(x)))
end_index <- dim(reordered_mean_realisation)[2]
for (i in 1:length(id)) {
  reordered_mean_realisation[i, ] <- normalised_output[i, c(start_index[i]:end_index, 1:(start_index[i]-1))]
}
one_output <- reordered_mean_realisation[cluster_membership$cluster == 1, ]
mean_one_output <- apply(one_output, 2, mean)

two_output <- reordered_mean_realisation[cluster_membership$cluster == 2, ]
mean_two_output <- apply(two_output, 2, mean)

aug_norm_output <- rbind(normalised_output, mean_one_output, mean_two_output)

summary_probs <- array(dim = c(length(id) + 2, 12, 3))
p_overall <- 1
for (k in 1:67) {
  
  # Generating Monthly Relative Probabilities of Finding Stephensi and Multiplying By Defined Overall Prob Detection Given Present
  time_series <- aug_norm_output[k, ]
  avg_month_prob <- vector(mode = "numeric", length = 12L)
  counter <- 1
  for (i in 1:length(avg_month_prob)) {
    if (i < 12) {
      avg_month_prob[i] <- mean(time_series[c(counter, counter + 1)])
      counter <- counter + 2
    } else {
      avg_month_prob[i] <- mean(time_series[c(counter, counter + 1, counter + 2)])
      counter <- counter + 2
    }
  }
  p_zero_caught <- 1 - (avg_month_prob/max(avg_month_prob) * p_overall) # option1 <- Alt = 1 - (test * 0.5)
  
  # Iterating Over All Possible Start Months and All Possible Number of Consecutive Months Sampled
  #   rows are the month we start at
  #   columns are the number of months we consider
  prob_matrix <- matrix(data = NA, nrow = 12, ncol = 12)
  for (i in 1:12) {
    for (j in 1:12) {
      if ((i+j-1) > 12) {
        indices <- c(i:12, 1:(j - length(i:12)))
      } else {
        indices <- i:(i+j-1)
      }
      temp_p <- prod(p_zero_caught[indices])
      prob_matrix[i, j] <- temp_p
    }
  }
  temp_median_prob <- apply(prob_matrix, 2, quantile, prob = c(0.25, 0.5, 0.75))
  summary_probs[k, , ] <- t(temp_median_prob)
  print(k)
}

lower_probs <- summary_probs[1:65, , 1]
median_probs <- summary_probs[1:65, , 2]
upper_probs <- summary_probs[1:65, , 3]

cluster_one_median <- data.frame(time = seq(1, 12), median = summary_probs[66, , 2], lower = summary_probs[66, , 1], upper = summary_probs[66, , 3])
cluster_two_median <- data.frame(time = seq(1, 12), median = summary_probs[67, , 2], lower = summary_probs[67, , 1], upper = summary_probs[67, , 3])

median_summary <- data.frame(id = rep(id, 3), stat = c(rep("med", 65), rep("low", 65), rep("high", 65)),
                             cluster = rep(cluster_membership$cluster, each = 3), 
                             rbind(median_probs, lower_probs, upper_probs)) %>%
  pivot_longer(cols = X1:X12, names_to = "time", values_to = "prob") %>%
  mutate(time = as.numeric(gsub("X", "", time))) %>%
  group_by(time, stat, cluster) %>%
  summarise(median = median(prob),
            lower = min(prob),
            upper = max(prob)) %>%
  pivot_wider(names_from = "stat", 
              values_from = median:upper)#
medians <- data.frame(id = id, median_probs) %>%
  pivot_longer(cols = X1:X12, names_to = "time", values_to = "prob") %>%
  mutate(time = as.numeric(gsub("X", "", time)))

fig4c <- ggplot(median_summary) +
  geom_ribbon(aes(x = time, ymin = lower_med, ymax = upper_med, fill = factor(cluster)), alpha = 0.1) +
  geom_path(aes(x = time, y = median_med, colour = factor(cluster)), size = 1) +
  scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12)) +
  scale_fill_manual(values = palette()[2:1]) + 
  scale_colour_manual(values = palette()[2:1]) + 
  labs(y = "Probability of Missing Anopheles Stephensi", x = "Number of Consecutive Months Sampled") +
  theme_bw() +
  theme(legend.position = "none")

#fig4c <- ggplot(median_summary) +
  #geom_ribbon(aes(x = time, ymin = lower_med, ymax = upper_med), fill = "dark grey", alpha = 0.2) +
  #geom_ribbon(aes(x = time, ymin = lower_med, ymax = upper_med), fill = "#70798C", alpha = 0.2) +
  #geom_path(aes(x = time, y = median_med), size = 1, colour = "#70798C") +
  # geom_path(data = cluster_one_median, aes(x = time, y = median), colour = palette()[2], size = 1) +
  # geom_ribbon(data = cluster_one_median, aes(x = time, ymin = lower, ymax = upper), fill = palette()[2], alpha = 0.2) +
  # geom_path(data = cluster_two_median, aes(x = time, y = median), colour = palette()[1], size = 1) +
  # geom_ribbon(data = cluster_two_median, aes(x = time, ymin = lower, ymax = upper), fill = palette()[1], alpha = 0.2) +
  # geom_path(data = medians, aes(x = time, y = prob, group = id), alpha = 0.2) +
  # scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12)) +
  # scale_y_continuous(limits = c(0, 1)) +
  # labs(y = "Probability of Missing Anopheles Stephensi", x = "Number of Consecutive Months Sampled") +
  # theme_bw()

fig4subset <- cowplot::plot_grid(fig4b, fig4c, ncol = 1)
fig4overall <- cowplot::plot_grid(fig4subset, fig4a, rel_widths = c(1, 1.3)) ## 9.5 * 6.25 (h x w)

fig4b <- fig4b +
  theme(legend.position = "right")
fig4a <- fig4a +
  theme(legend.position = c(0.85, 0.80))
fig4subset <- cowplot::plot_grid(fig4b, fig4c, ncol = 1, align = "v", axis = "l")
fig4overall <- cowplot::plot_grid(fig4a, fig4subset, rel_widths = c(1.25, 1)) ## 9.5 * 6.25 (h x w)
ggsave(fig4overall, file = here("figures", "Fig4_Overall.pdf"), width = 9.5, height = 6.15)

######################################################################
######################################################################

# alt_summary <- data.frame(id = rep(id, 3), rbind(median_probs, lower_probs, upper_probs)) %>%
#   pivot_longer(cols = X1:X12, names_to = "time", values_to = "prob") %>%
#   mutate(time = as.numeric(gsub("X", "", time))) %>%
#   group_by(time) %>%
#   summarise(median = median(prob),
#             lower = min(prob),
#             upper = max(prob)) 
# alt_summary_plot <- ggplot(alt_summary) +
#   geom_path(aes(x = time, y = median)) +
#   geom_ribbon(aes(x = time, ymin = lower, ymax = upper), alpha = 0.2) +
#   scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12)) +
#   labs(y = "Probability of Missing Anopheles Stephensi",
#        x = "Number of Months Sampled")
# 
# # Monthly Surveillance Considerations
# unsmoothed_counts <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds")) %>%
#   dplyr::select(id, Jan:Dec) %>%
#   pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "catch") %>%
#   group_by(id) %>%
#   summarise(num_months = sum(!is.na(catch)),
#             total_catch = sum(catch, na.rm = TRUE),
#             avg_catch = total_catch/num_months)
# 
# overdisp_df <- data.frame(overdisp = overdisp, catch = log(unsmoothed_counts$avg_catch))
# 
# linear_with_int <- glm(overdisp ~ catch, data = overdisp_df)
# linear_int_coefs <- coef(linear_with_int)
# 
# linear_no_int <- glm(overdisp ~ 0 + catch, data = overdisp_df)
# linear_no_coefs <- coef(linear_no_int)
# 
# test_catches <- seq(0, 10, 1)
# plot(log(unsmoothed_counts$avg_catch), overdisp, ylim = c(0, 6), xlim = c(0, 6.5))
# lines(test_catches, linear_int_coefs[1] + linear_int_coefs[2] * test_catches)
# lines(test_catches, linear_no_coefs[1] * test_catches, col = "red")
# 
# perc <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1)
# catch <- c(5, 10, 25, 50, 100, 200, 500)
# 
# mean(unsmoothed_counts$avg_catch)
# median(unsmoothed_counts$avg_catch)
# 
# inputs <- data.frame(catch = log(catch))
# with_int_preds <- predict(linear_with_int, inputs)
# no_int_preds <- predict(linear_no_int, inputs)
# probs <- array(dim = c(length(perc), length(catch), 2)) # rows = catch size, cols = % pop = stephensi
# for (i in 1:length(perc)) {
#   for (j in 1:length(catch)) {
#     probs[i, j, 1] <- 1 - dnbinom(0, mu = perc[i] * catch[j], size = no_int_preds[j])
#     probs[i, j, 2] <- 1 - dnbinom(0, mu = perc[i] * catch[j], size = with_int_preds[j])
#   }
# }
# 
# dnbinom(0, mu = 1, size = 3)
# dnbinom(0, mu = 1, size = 0.5)
# 
# x <- probs[, , 1]
# y <- probs[, , 2]
# 
# temp_rest <- rnbinom(25000, mu = 99, size = no_int_preds[1])
# 
# 
# inputs <- data.frame(catch = c(1, 2, 3, 5))
# with_int_preds <- predict(linear_with_int, inputs)
# no_int_preds <- predict(linear_no_int, inputs)
# 
# prop_overall_population <- 0.01
# 
# x <- data.frame(catch = seq(1, 100000, 1))
# y <- predict(linear_with_int, log(x))
# preds <- dnbinom(0, mu = x$catch, size = y)
# plot(log(x$catch), preds, type = "l", ylim = c(0, 1))
# 
# x <- data.frame(catch = seq(1, 100000, 1))
# y <- predict(linear_no_int, log(x))
# preds <- dnbinom(0, mu = x$catch, size = y)
# lines(log(x$catch), preds, col = "red")
# 
# 
# 
# 
# 
# 
# num <- 25000
# mu <- c(10, 100, 1000, 100000)
# ind <- vector(mode = "numeric", length = 3L)
# par(mfrow = c(1, 2))
# for (i in 1:4) {
#   temp <- rnbinom(num, mu = mu[i], size = no_int_preds[i])
#   ordered_temp <- temp[rev(order(temp))]
#   summed_temp <- vector(mode = "numeric", length = num)
#   for (j in 1:num) {
#     summed_temp[j] <- sum(ordered_temp[1:j])/sum(ordered_temp)
#   }
#   if (i == 1) {
#     plot(summed_temp, type = "l")
#   } else {
#     lines(summed_temp, type = "l")
#   }
#   temp_ind <- which(abs(summed_temp - 0.8) == min(abs((summed_temp - 0.8))))
#   ind[i] <- temp_ind/num
# }
# 
# test <- rnbinom(num, mu = 100, size = 3.34)
# #hist(test)
# ordered_test <- test[rev(order(test))]
# summed <- vector(mode = "numeric", length = num)
# for (i in 1:num) {
#   summed[i] <- sum(ordered_test[1:i])/sum(ordered_test)
# }
# plot(summed)
# 
# 
# x <- seq(0, 1, 0.001)
# plot(qnbinom(x, mu = 100, size = 1.93), x)
# 
# test <- rnbinom(10000, mu = 100, size = 3.34)
# hist(test)
# ordered_test <- test[rev(order(test))]
# 
# a <- ecdf(ordered_test)
# plot(a)
# quantile(a, 0.8)
# 
# par(mfrow = c(1, 2))
# num <- 10000
# test <- rnbinom(num, mu = 100, size = 3.34)
# #hist(test)
# ordered_test <- test[rev(order(test))]
# summed <- vector(mode = "numeric", length = num)
# for (i in 1:num) {
#   summed[i] <- sum(ordered_test[1:i])/sum(ordered_test)
# }
# plot(summed)
# 
# test <- rnbinom(num, mu = 100, size = 0.5)
# #hist(test)
# ordered_test <- test[rev(order(test))]
# summed <- vector(mode = "numeric", length = num)
# for (i in 1:num) {
#   summed[i] <- sum(ordered_test[1:i])/sum(ordered_test)
# }
# ind <- which((summed - 0.8) == min(abs((summed - 0.8))))
# 
# ind/num
# 
# 
# pnorm(1.96, 0, 1)
# qnorm(0.975, 0, 1) # qnorm = CFDF
# 
# rnbinom(1, mu = 1, size = 4)
# 
# hist(rnbinom(10000, mu = 10, size = 2))
