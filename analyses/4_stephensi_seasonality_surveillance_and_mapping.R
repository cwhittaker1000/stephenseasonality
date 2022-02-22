# Loading Required Libraries
library(tidyverse); library(sf); library(tidymodels)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here) 

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
  scale_fill_gradient2(low = muted("red"), mid = "white", high = muted("blue"),
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
        legend.position = "right", 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold")) +
  guides(fill = guide_colourbar(title = "Prob.Single\nSeasonal Peak", ticks = FALSE))