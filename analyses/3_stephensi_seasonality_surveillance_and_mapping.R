# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels); library(cowplot)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(zoo); library(RColorBrewer)

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
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2")) %>%
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
  dplyr::group_by(x, y) %>%
  dplyr::summarise(layer = mean(layer, na.rm = TRUE))
fig3a <- ggplot() +
  geom_tile(data = all_rasters, aes(x = x, y = y, fill = layer)) +
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

# Simulating Strategies With Sequential Monthly Sampling 
EIR <- 1
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

# Heatmaps Exploring Range of Sampling Effort
#     1st dim is the time series
#     2nd dim is the number of months sampled
#     3rd dim is the number of nights per month sampled
#     4th dim is which month you start at 
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

max <- max(apply(not_sampled_vector_peak_absolute_summary, c(2, 3), mean))
size <- 12

temp1 <- apply(not_sampled_annual_avg_absolute_summary, c(2, 3), mean)
temp1 <- data.frame(temp1)
temp1$Y <- paste0("Y", 1:dim)
temp1$Y <- factor(temp1$Y, levels = paste0("Y", 1:dim))
temp1$sampling_start <- "Annual Average"

temp2 <- apply(not_sampled_rainfall_peak_absolute_summary, c(2, 3), mean)
temp2 <- data.frame(temp2)
temp2$Y <- paste0("Y", 1:dim)
temp2$Y <- factor(temp2$Y, levels = paste0("Y", 1:dim))
temp2$sampling_start <- "Rainfall Peak"

temp3 <- apply(not_sampled_vector_peak_absolute_summary, c(2, 3), mean)
temp3 <- data.frame(temp3)
temp3$Y <- paste0("Y", 1:dim)
temp3$Y <- factor(temp3$Y, levels = paste0("Y", 1:dim))
temp3$sampling_start <- "Vector Peak"

temp3[, -c(6:7)]/temp1[, -c(6:7)]
temp3[, -c(6:7)]/temp2[, -c(6:7)]
temp2[, -c(6:7)]/temp1[, -c(6:7)]

temp_long <- rbind(temp1, temp2, temp3) %>%
  pivot_longer(cols = -c(Y, sampling_start), names_to = "X")
temp_long$X  <- factor(temp_long$X, levels = unique(temp_long$X))
temp_long$sampling_start  <- factor(temp_long$sampling_start, 
                                    levels = c("Annual Average", "Rainfall Peak", "Vector Peak"))

sampling_heatmaps <- ggplot(temp_long) +
  geom_tile(aes(x = Y, y = X, fill = value)) +
  scale_y_discrete(position = "left", labels = 1:dim) +
  scale_fill_viridis_c(option = "magma", name = "Prob. of\nDetection",
                       limits = c(0, max)) +
  theme(plot.title = element_text(size=size)) +
  scale_x_discrete(labels =  1:dim) +
  facet_wrap(~sampling_start) +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = 2.5, ymax = 3.5),
            fill = NA, col = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(y = "Sampling Days Per Month", x = "Number of Months Sampled") 

# Cumulative Probability Plots
#   1st dim is the time series
#   2nd dim is the number of months sampled
#   3rd dim is the number of nights per month sampled
#   4th dim is which month you start at 
not_sampled_vector_peak_absolute_summary <- array(data = NA, dim = c(65, 12, 3)) # 1-12 consecutive months sampled, 1-3 nights per month
not_sampled_rainfall_peak_absolute_summary <- array(data = NA, dim = c(65, 12, 3))
not_sampled_annual_avg_absolute_summary <- array(data = NA, dim = c(65, 12, 3))
for (i in 1:65) {
  ts_peak_vector_month <- peak_density_month[i] 
  ts_peak_rainfall_month <- peak_rainfall_month[i] 
  
  not_sampled_vector_peak_absolute_summary[i, 1:12, 1:3] <- 1- prob_not_sampled_poisson[i, 1:12, 1:3, ts_peak_vector_month]
  not_sampled_rainfall_peak_absolute_summary[i, 1:12, 1:3] <- 1 - prob_not_sampled_poisson[i, 1:12, 1:3, ts_peak_rainfall_month]
  not_sampled_annual_avg_absolute_summary[i, 1:12, 1:3] <- 1 - apply(prob_not_sampled_poisson[i, 1:12, 1:3, 1:12], c(1, 2), mean)
}

x <- not_sampled_vector_peak_absolute_summary[, , 3]
mean <- apply(x, 2, mean)
x <- data.frame(id = 1:65, cluster = overall$cluster, x) %>%
  pivot_longer(-c(id, cluster), names_to = "months_sampled", values_to = "prob_detect")
x$months_sampled <- factor(x$months_sampled, levels = paste0("X", 1:12))
x$id <- factor(x$id)
x$cluster <- factor(x$cluster)
x$months_sampled <- as.numeric(gsub("X", "", x$months_sampled))
temp <- data.frame(id = 1:65, cluster = overall$cluster, months_sampled = 0, prob_detect = 0)
x <- rbind(x, temp)
meanx <- x %>%
  group_by(cluster, months_sampled) %>%
  summarise(mean = mean(prob_detect))

y <- not_sampled_rainfall_peak_absolute_summary[, , 3]
mean <- apply(y, 2, mean)
y <- data.frame(id = 1:65, cluster = overall$cluster, y) %>%
  pivot_longer(-c(id, cluster), names_to = "months_sampled", values_to = "prob_detect")
y$months_sampled <- factor(y$months_sampled, levels = paste0("X", 1:12))
y$id <- factor(y$id)
y$cluster <- factor(y$cluster)
y$months_sampled <- as.numeric(gsub("X", "", y$months_sampled))
temp <- data.frame(id = 1:65, cluster = overall$cluster, months_sampled = 0, prob_detect = 0)
y <- rbind(y, temp)
meany <- y %>%
  group_by(cluster, months_sampled) %>%
  summarise(mean = mean(prob_detect))

x$sampling_start <- "Vector\nPeak"
y$sampling_start <- "Rainfall\nPeak"
meanx$sampling_start <- "Vector\nPeak"
meany$sampling_start <- "Rainfall\nPeak"
z <- rbind(x, y)
z$cluster <- ifelse(z$cluster == 1, "Cluster 1", "Cluster 2")
meanz <- rbind(meanx, meany)
meanz$cluster <- ifelse(meanz$cluster == 1, "Cluster 1", "Cluster 2")

meanz2 <- meanz %>%
  filter(months_sampled == 3) %>%
  group_by(cluster, sampling_start) 

cum_prob_plots <- ggplot() +
  geom_line(data = z, aes(x = months_sampled, y = prob_detect, group = interaction(id, sampling_start), 
                          col = sampling_start), alpha = 0.2) +
  geom_line(data = meanz, aes(x = months_sampled, y = mean, col = sampling_start), alpha = 1, size = 2) +
  facet_wrap(~cluster) +
  scale_colour_manual(values = c("#2589BD", "#FC9F5B")) +
  lims(x = c(0, 12), y = c(0, 1)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  labs(x = "Number of Months Sampled Following Peak",
       y = "Probability of Successful Detection") +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "right") + 
  guides(col = guide_legend(title = "Sampling\nStarts", nrow = 2))

fig3bc <- plot_grid(sampling_heatmaps, cum_prob_plots, nrow = 2, axis = "lr", align = "v")

fig3 <- plot_grid(fig3a, fig3bc, rel_widths = c(1, 1.5))
fig3
ggsave(here("figures", "Fig3_Overall.pdf"), fig4, width = 12.5, height = 6)
