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

dji_shp <- st_as_sf(getData('GADM', country = 'DJI', level = 0))
eth_shp <- st_as_sf(getData('GADM', country = 'ETH', level = 0))
eri_shp <- st_as_sf(getData('GADM', country = 'ERI', level = 0))
sdn_shp <- st_as_sf(getData('GADM', country = 'SDN', level = 0))
som_shp <- st_as_sf(getData('GADM', country = 'SOM', level = 0))
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
predictions <- predict(final_random_forest_fit_ups, juiced_HOA[1:length(dji$ISO), ], "prob")
peaks <- ifelse(predictions$.pred_one >= 0.5, "one", "two")

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


extent(dji_shp)
extent(eth_shp)
extent(eri_shp)
extent(sdn_shp)
extent(som_shp)

xmin <- 21
xmax <- 48
ymin <- -2
ymax <- 15

overall_raster <- as.data.frame(as(c(dji_raster, eth_raster, eri_raster, sdn_raster, som_raster), "SpatialPixelsDataFrame"))

x <- list(dji_raster, eth_raster)
x$fun <- mean
x$na.rm <- TRUE
x$tolerance <- 100
combined_raster <- do.call(merge, x)

snap <- raster(resolution = c(0.09748167,0.09746483), xmn = 21, xmx = 60, ymn = -2, ymx = 25)
snap$layer <- NA
dji_file <- projectRaster(dji_raster, snap, method = "bilinear")
eth_file <- projectRaster(eth_raster, snap, method = "bilinear")
eri_file <- projectRaster(eri_raster, snap, method = "bilinear")
sdn_file <- projectRaster(sdn_raster, snap, method = "bilinear")
som_file <- projectRaster(som_raster, snap, method = "bilinear")

plot(dji_file)
plot(eth_file)
plot(eri_file)
plot(sdn_file)
plot(som_file)

x <- list(dji_file, eth_file, eri_file, sdn_file, som_file)
x$fun <- mean
x$na.rm <- TRUE
x$tolerance <- 1
combined_raster <- do.call(merge, x)
plot(combined_raster)

raster_plot <- as.data.frame(as(combined_raster, "SpatialPixelsDataFrame"))
ggplot() +
  geom_tile(data = raster_plot, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = hoa_outline, col = "black", fill = NA, size = 1.2) 

dji_test <- dji_raster
dji_test$layer[is.na(dji_test$layer)] <- 0
dji_file <- projectRaster(dji_test, snap, method = "bilinear")
plot(dji_file)

snap[1:48918]
dji_raster$layer[1:306]
dji_test$layer[1:306]
dji_file[1:48198]

sum(dji_raster$layer[1:306], na.rm = TRUE)
sum(dji_file$layer[1:48198], na.rm = TRUE)

eth_file <- projectRaster(eth_raster, snap, method = "ngb")
x <- list(dji_file, eth_file)
names(x) <- NULL
x$fun <- mean
mos <- do.call(raster::mosaic, x)
plot(mos)

mos_plot <- as.data.frame(as(mos, "SpatialPixelsDataFrame"))
mos_plot$lat <- som$lat
mos_plot$lon <- som$lon

ggplot() +
  geom_tile(data = overall_raster, aes(x = lon, y = lat, fill = layer))


combined_raster <- Reduce("+", x)

origin(dji_raster) <- origin(eth_raster)


overall_raster <- rbind(dji_raster_plot, eth_raster_plot, eri_raster_plot, sdn_raster_plot, som_raster_plot)

as(dji_raster, "SpatialPixelsDataFrame")
as(rbind(dji_raster), "SpatialPixelsDataFrame")

typeof(dji_raster_plot)
class(dji_raster_plot)
typeof(rbind(dji_raster_plot))
class()

ggplot() +
  geom_tile(data = overall_raster, aes(x = lon, y = lat, fill = layer))

ggplot() +
  geom_tile(data = dji_raster_plot, aes(x = lon, y = lat, fill = layer)) +
  geom_tile(data = eth_raster_plot, aes(x = lon, y = lat, fill = layer)) +
  geom_tile(data = eri_raster_plot, aes(x = lon, y = lat, fill = layer)) +
  geom_tile(data = sdn_raster_plot, aes(x = lon, y = lat, fill = layer)) +
  geom_tile(data = som_raster_plot, aes(x = lon, y = lat, fill = layer)) +
  scale_fill_gradient2(low = muted("red"),
                       mid = "white",
                       high = muted("blue"),
                       midpoint = 0.5,
                       limits = c(0, 1),
                       space = "Lab",
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "colour") +
  geom_sf(data = hoa_outline, col = "black", fill = NA, size = 1.2) 


df_spat <- data.table(peaks = peaks, lon = juiced_HOA$lon, lat = juiced_HOA$lat,
                      pred_one = predictions$.pred_one, log_pop = juiced_HOA$population_per_1km)
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
  

ggplot() +
  geom_sf(data = DT_sf, aes(col = log_pop))

