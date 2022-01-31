# Loading Required Libraries
library(tidyverse); library(sf); library(data.table)

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
dji <- read.csv("data/environmental_covariates/DJI_environment_data.csv") %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol) %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0, LC_122 = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
eth <- read.csv("data/environmental_covariates/ETH_environment_data.csv") %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol) %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
eri <- read.csv("data/environmental_covariates/ERI_environment_data.csv") %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
sdn <- read.csv("data/environmental_covariates/SDN_environment_data.csv") %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
som <- read.csv("data/environmental_covariates/SOM_environment_data.csv") %>%
  dplyr::select(-elevation, -goat, -sheep, -irrigation, -nrow, -ncol)  %>%
  rename(population_per_1km = population_density,
         precipitation_coldest_quarter = worldclim_19,
         precipitation_seasonality_cv = worldclim_15,
         temperature_seasonality = worldclim_4) %>%
  mutate(country_India = 0, country_Iran = 0) %>%
  dplyr::select(lat, lon, vars_from_rf)
hoa <- rbind(dji, eth, eri, sdn, som)
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

## Load In The Random Forest Objects
## load em in
## fit_final_random_forest_fit_ups <- _________
juiced_HOA$country_peaks <- "bloop" # add dummy category

iterations_ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
final_random_forest_fit_ups <- iterations_ups$model[[1]]
predictions <- predict(final_random_forest_fit_ups, juiced_HOA[4030, ], "prob")

peaks <- ifelse(predictions$.pred_one >= 0.5, "one", "two")
df_spat <- data.table(peaks = peaks, lon = juiced_HOA$lon, lat = juiced_HOA$lat,
                      pred_one = predictions$.pred_one, log_pop = juiced_HOA$population_per_1km)
DT_sf <- st_as_sf(df_spat, coords = c("lon", "lat"), crs = 4326)
ggplot(DT_sf) +
  geom_sf(aes(col = pred_one))

ggplot(DT_sf) +
  geom_sf(aes(col = log_pop))

