# Loading Required Libraries
library(tidyverse); library(sf); library(data.table); library(tidymodels)
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here)

# Loading in geographical data
dji <- read.csv("data/environmental_covariates/DJI_environment_data.csv") 
eth <- read.csv("data/environmental_covariates/ETH_environment_data.csv")
eri <- read.csv("data/environmental_covariates/ERI_environment_data.csv")
sdn <- read.csv("data/environmental_covariates/SDN_environment_data.csv")
som <- read.csv("data/environmental_covariates/SOM_environment_data.csv")

dji_shp <- st_as_sf(getData('GADM', country = 'DJI', level = 0, path = here("data/admin_units")))
eth_shp <- st_as_sf(getData('GADM', country = 'ETH', level = 0, path = here("data/admin_units")))
eri_shp <- st_as_sf(getData('GADM', country = 'ERI', level = 0, path = here("data/admin_units")))
sdn_shp <- st_as_sf(getData('GADM', country = 'SDN', level = 0, path = here("data/admin_units")))
som_shp <- st_as_sf(getData('GADM', country = 'SOM', level = 0, path = here("data/admin_units")))
hoa_outline <- rbind(dji_shp, eth_shp, eri_shp, sdn_shp, som_shp)

# Loading in Location Specific Data and Training the Recipe
ts_metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE)))
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))

data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(peaks, country, population_per_1km:mean_temperature_driest_quarter, -LC_190, -LC_210, -mean_temperature_driest_quarter) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country_peaks = paste0(peaks, "_", country))
data$peaks <- ifelse(data$peaks == 1, "one", "two")
data$peaks <- as.factor(data$peaks)
data$population_per_1km <- log(data$population_per_1km)

# Processing Arran's Horn of Africa Country Rasters
vars_from_rf <- c("LC_10", "LC_11", "LC_110", "LC_120", "LC_122", "LC_130", "LC_150", "LC_180", "LC_20", "LC_30",
                  "population_per_1km", "precipitation_coldest_quarter", # "country_India", "country_Iran", # don't forget about these
                  "precipitation_seasonality_cv", "temperature_seasonality")
dji_dat <- dji %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol) %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0, LC_122 = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
eth_dat <- eth %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol) %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
eri_dat <- eri %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
sdn_dat <- sdn %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
som_dat <- som %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
hoa <- rbind(dji_dat, eth_dat, eri_dat, sdn_dat, som_dat)
hoa$population_per_1km <- log(hoa$population_per_1km)

# Regenerating the same recipe used in the Random Forest Fitting to Process this New Data
data_for_recipe_recapit <- data %>%
  dplyr::select(peaks, vars_from_rf) %>%
  mutate(lat = 0, lon = 0)
set.seed(915)
envt_recipe_ups <- recipe(peaks  ~ ., data = data_for_recipe_recapit) %>% # go and manually check that this gives same results on *data* object as the full recipe does!!
  step_center(all_numeric(), -lat, -lon) %>% 
  step_scale(all_numeric(), -lat, -lon) 
envt_prepped_ups <- prep(envt_recipe_ups, training = data_for_recipe_recapit, verbose = TRUE)
juiced_HOA <- bake(envt_prepped_ups, hoa) %>%
  dplyr::select(lat, lon, everything()) %>%
  mutate(country_India = 0, country_Iran = 0)

# Load In The Random Forest Objects and Generate Predictions
juiced_HOA$country_peaks <- "bloop" # add dummy category
iterations_ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
final_random_forest_fit_ups <- iterations_ups$model[[1]]
predictions <- predict(final_random_forest_fit_ups, juiced_HOA, "prob")
peaks <- ifelse(predictions$.pred_one >= 0.5, "one", "two")

# Plotting as Non-Raster
df_spat <- data.table(peaks = peaks, lon = juiced_HOA$lon, lat = juiced_HOA$lat, pred_one = predictions$.pred_one, log_pop = juiced_HOA$population_per_1km)
DT_sf <- st_as_sf(df_spat, coords = c("lon", "lat"), crs = 4326)
ggplot() +
  geom_sf(data = DT_sf, aes(col = pred_one)) +
  geom_sf(data = hoa_outline, col = "black", fill = NA, size = 1.2) +
  scale_colour_gradient2(low = muted("red"),
                         mid = "white",
                         high = muted("blue"),
                         midpoint = 0.5,
                         limits = c(0, 1),
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "colour")

# Plotting as Rasters
dji_predictions <- predict(final_random_forest_fit_ups, juiced_HOA[1:length(dji$ISO), ], "prob")
dji_raster <- raster(ncol = dji$ncol[1], nrow = dji$nrow[1], ext = extent(dji_shp))
dji_raster[dji$cell] <- dji_predictions$.pred_one
dji_raster_plot <- as.data.frame(as(dji_raster, "SpatialPixelsDataFrame"))
dji_raster_plot$lat <- dji$lat
dji_raster_plot$lon <- dji$lon

eth_predictions <- predict(final_random_forest_fit_ups, juiced_HOA[(length(dji$ISO) + 1):(length(dji$ISO) + length(eth$ISO)), ], "prob")
eth_raster <- raster(ncol = eth$ncol[1], nrow = eth$nrow[1], ext = extent(eth_shp))
eth_raster[eth$cell] <- eth_predictions$.pred_one
eth_raster_plot <- as.data.frame(as(eth_raster, "SpatialPixelsDataFrame"))
eth_raster_plot$lat <- eth$lat
eth_raster_plot$lon <- eth$lon

eri_predictions <- predict(final_random_forest_fit_ups, juiced_HOA[(length(dji$ISO) + length(eth$ISO) + 1):(length(dji$ISO) + length(eth$ISO) + length(eri$ISO)), ], "prob")
eri_raster <- raster(ncol = eri$ncol[1], nrow = eri$nrow[1], ext = extent(eri_shp))
eri_raster[eri$cell] <- eri_predictions$.pred_one
eri_raster_plot <- as.data.frame(as(eri_raster, "SpatialPixelsDataFrame"))
eri_raster_plot$lat <- eri$lat
eri_raster_plot$lon <- eri$lon

sdn_predictions <- predict(final_random_forest_fit_ups, juiced_HOA[(length(dji$ISO) + length(eth$ISO) + length(eri$ISO) + 1):(length(dji$ISO) + length(eth$ISO) + length(eri$ISO) + length(sdn$ISO)), ], "prob")
sdn_raster <- raster(ncol = sdn$ncol[1], nrow = sdn$nrow[1], ext = extent(sdn_shp))
sdn_raster[sdn$cell] <- sdn_predictions$.pred_one
sdn_raster_plot <- as.data.frame(as(sdn_raster, "SpatialPixelsDataFrame"))
sdn_raster_plot$lat <- sdn$lat
sdn_raster_plot$lon <- sdn$lon

som_predictions <- predict(final_random_forest_fit_ups, juiced_HOA[(length(dji$ISO) + length(eth$ISO) + length(eri$ISO) + length(sdn$ISO) + 1):(length(dji$ISO) + length(eri$ISO) + length(eth$ISO) + length(sdn$ISO) + length(som$ISO)), ], "prob")
som_raster <- raster(ncol = som$ncol[1], nrow = som$nrow[1], ext = extent(som_shp))
som_raster[som$cell] <- som_predictions$.pred_one
som_raster_plot <- as.data.frame(as(som_raster, "SpatialPixelsDataFrame"))
som_raster_plot$lat <- som$lat
som_raster_plot$lon <- som$lon

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

# lat and lon are wrong way round currently
raster_plot <- as.data.frame(as(combined_raster, "SpatialPixelsDataFrame"))
ggplot() +
  geom_tile(data = raster_plot, aes(x = x, y = y, fill = layer)) +
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
