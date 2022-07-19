#######################################################################################################
##                                                                                                   ##
##                            Loading Required Libraries and Functions                               ##
##                                                                                                   ##
#######################################################################################################
library(dplyr); library(raster); library(rgdal); library(sf); library(raster); library(sp); 
library(tidyverse); library(here); library(rgee);  library(reticulate); library(rgeos); library(maptools); 
library(tidyr); library(maps); library(scales); library(ggplot2); library(zoo); library(cowplot); 
library(ggspatial); library(rnaturalearthdata); library(rnaturalearth); library(ggrepel)

# Loading in admin 0 units, simplifying and saving
admin0 <- readRDS(here("data", "admin_units", "raw_gadm_shapefiles", "level0.rds"))
admin0 <- admin0[!is.na(admin0$NAME_0), ]
countries <- c("Afghanistan", "Armenia", "Azerbaijan", "Bangladesh", "Bhutan", "Cambodia", "China", "Djibouti", "Egypt", 
               "Ethiopia", "India", "Israel", "Iraq", "Iran", "Jordan", "Kenya", "Kyrgyzstan",
               "Laos", "Lebanon", "Myanmar", "Nepal", "Oman", "Pakistan", "Saudi Arabia", 
               "Somalia", "South Sudan", "Sudan", "Syria", "Tajikistan", "Thailand", "Turkey", "Turkmenistan", "Uzbekistan",
               "Vietnam", "Yemen",
               "Mongolia", "Georgia", "Kazakhstan", "Turkey", "Ukraine", "Romania", "Bulgaria", "Libya",
               "Democratic Republic of the Congo", "Central African Republic", "Rwanda", "Burundi", "Tanzania",
               "Uganda", "Chad", "Malaysia", "Singapore", "Indonesia", "Uzbekistan", "Russia")

admin0 <- admin0[admin0$NAME_0 %in% countries, ]
saveRDS(admin0, here("data", "admin_units", "complex_admin0.rds"))
for (i in 1:nrow(admin0)) {
  temp <- admin0$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 5000)
  admin0$geometry[i] <- reduce
}
saveRDS(admin0, here("data", "admin_units", "simplified_admin0.rds"))

# Loading in admin 1 units, simplifying and saving
admin1 <- readRDS(here("data", "admin_units", "raw_gadm_shapefiles", "level1.rds"))
admin1 <- admin1[!is.na(admin1$NAME_0), ]
countries <- c("Afghanistan", "Djibouti", "Ethiopia", "India",  "Iran", "Myanmar", "Pakistan", "Saudi Arabia")
admin1 <- admin1[admin1$NAME_0 %in% countries, ]
saveRDS(admin1, here("data", "admin_units", "complex_admin1.rds"))
for (i in 1:nrow(admin1)) {
  temp <- admin1$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 3500)
  admin1$geometry[i] <- reduce
}
saveRDS(admin1, here("data", "admin_units", "simplified_admin1.rds"))

# Loading in admin 2 units, simplifying and saving
admin2 <- readRDS(here("data", "admin_units", "raw_gadm_shapefiles", "level2.rds"))
admin2 <- admin2[!is.na(admin2$NAME_0), ]
countries <- c("Afghanistan", "Djibouti", "Ethiopia", "India",  "Iran", "Myanmar", "Pakistan", "Saudi Arabia")
admin2 <- admin2[admin2$NAME_0 %in% countries, ]
saveRDS(admin2, here("data", "admin_units", "complex_admin2.rds"))
for (i in 1:nrow(admin2)) {
  temp <- admin2$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 1000)
  admin2$geometry[i] <- reduce
}
saveRDS(admin2, here("data", "admin_units", "simplified_admin2.rds"))

# Saving Admin 1 and Admin2 unit names for each country to assist with geolocation
admin1_csv <- data.frame(country = admin1$NAME_0, admin1 = admin1$NAME_1, 
                         id = admin1$GID_1, type = admin1$TYPE_1, alt_names = admin1$VARNAME_1)

admin2_csv <- data.frame(country = admin2$NAME_0, 
                         admin1 = admin2$NAME_1, id1 = admin2$GID_1,
                         admin2 = admin2$NAME_2, id2 = admin2$GID_2,
                         type = admin2$TYPE_2, alt_names = admin2$VARNAME_2)

write.csv(admin1_csv, file = here("data", "admin_units", "admin1_details.csv"), row.names = FALSE)
write.csv(admin2_csv, file = here("data", "admin_units", "admin2_details.csv"), row.names = FALSE)

# Initialise rgee -  See https://github.com/r-spatial/rgee for more details about installation and getting things working
reticulate::py_version()
reticulate::conda_version()

library(rgee); library(reticulate)
# ee_clean_pyenv()
# reticulate::conda_list()
# reticulate::conda_remove("rgee")
# ee_install()
ee_Initialize()
# rgee::ee_install_upgrade()
# py_install("geemap")
ee_check() # Check non-R dependencies
gm <- import("geemap")

# Load metadata and admin 2 units
metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
admin2 <- readRDS(here("data", "admin_units", "simplified_admin2.rds"))
admin1 <- readRDS(here("data", "admin_units", "simplified_admin1.rds"))

# Extracting CHIRPS rainfall data
for (i in 1:nrow(metadata)) {
  
  # Get years that the study spans
  years <- seq(metadata$start[i], metadata$end[i], 1)
  num_years <- length(years)
  for (j in 1:length(years)) {
    if (years[j] < 1981) {
      years[j] <- 1981
    }
  }
  years <- unique(years)
  start_date <- paste0(years[1], "-01-01")
  end_date <- paste0(years[length(years)] + 1, "-01-01")
  
  # Get relevant admin 2 unit and convert to sf object
  if(!is.na(metadata$admin2[i])) {
    metadata_admin <- metadata$admin2[i]
    admin_2_data <- admin2[admin2$NAME_2 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_2_data)
  } else {
    metadata_admin <- metadata$admin1[i]
    admin_1_data <- admin1[admin1$NAME_1 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_1_data)
  }
  
  # Extract rainfall for that location
  chirps <- ee$ImageCollection("UCSB-CHG/CHIRPS/DAILY") %>%
    ee$ImageCollection$filterDate(start_date, end_date) %>%
    ee$ImageCollection$map(function(x) x$select("precipitation")) %>% # Select only precipitation bands
    ee$ImageCollection$toBands() # from imagecollection to image
  chirps_extract <- ee_extract(x = chirps, y = admin_geometry, sf = FALSE)
  
  # Process rainfall into right format
  first_rainfall_entry <- min(grep("precipitation", colnames(chirps_extract)))
  x <- chirps_extract %>%
    pivot_longer(cols = starts_with("X"), names_to = "day", values_to = "rainfall")
  x$day <- gsub("X", "", x$day)
  x$day <- gsub("_precipitation", "", x$day)
  x$day <- as.Date(x$day, format = "%Y%m%d")
  x$daymonth_id <- format(x$day, "%m-%d")
  days <- x$day
  x <- x %>%
    group_by(daymonth_id) %>%
    summarise(rainfall = mean(rainfall, na.rm = TRUE))
  x$day <- days[1:nrow(x)]
  x$time_series_id <- metadata$id[i]
  
  write.csv(x, file = here("data", "location_specific_rainfall", paste0("rainfall_ts", metadata$id[i], ".csv")), row.names = FALSE)
  
  print(c(i, dim(x)[1]))
  
}

# Extracting ERA5 2m Air Temperature
for (i in 1:nrow(metadata)) {
  
  # Get years that the study spans
  years <- seq(metadata$start[i], metadata$end[i], 1)
  num_years <- length(years)
  for (j in 1:length(years)) {
    if (years[j] < 1981) {
      years[j] <- 1981
    }
  }
  years <- unique(years)
  start_date <- paste0(years[1], "-01-01")
  end_date <- paste0(years[length(years)] + 1, "-01-01")
  
  # Get relevant admin 2 unit and convert to sf object
  if(!is.na(metadata$admin2[i])) {
    metadata_admin <- metadata$admin2[i]
    admin_2_data <- admin2[admin2$NAME_2 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_2_data)
  } else {
    metadata_admin <- metadata$admin1[i]
    admin_1_data <- admin1[admin1$NAME_1 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_1_data)
  }
  
  # Extract rainfall for that location
  air_temp <- ee$ImageCollection("ECMWF/ERA5/DAILY") %>%
    ee$ImageCollection$filterDate(start_date, end_date) %>%
    ee$ImageCollection$map(function(x) x$select("mean_2m_air_temperature")) %>% # Select only precipitation bands
    ee$ImageCollection$toBands() # from imagecollection to image
  air_temp_extract <- ee_extract(x = air_temp, y = admin_geometry, sf = FALSE)
  
  # Process rainfall into right format
  first_entry <- min(grep("mean_2m_air_temperature", colnames(air_temp_extract)))
  x <- air_temp_extract %>%
    pivot_longer(cols = starts_with("X"), names_to = "day", values_to = "mean_2m_air_temperature")
  x$day <- gsub("X", "", x$day)
  x$day <- gsub("_mean_2m_air_temperature", "", x$day)
  x$day <- as.Date(x$day, format = "%Y%m%d")
  x$daymonth_id <- format(x$day, "%m-%d")
  days <- x$day
  x <- x %>%
    group_by(daymonth_id) %>%
    summarise(mean_2m_air_temperature = mean(mean_2m_air_temperature, na.rm = TRUE))
  x$day <- days[1:nrow(x)]
  x$time_series_id <- metadata$id[i]
  x$mean_2m_air_temperature_celsius <- x$mean_2m_air_temperature - 273.15 # convert kelvin to celsius
  
  write.csv(x, file = here("data", "location_specific_airtemp", paste0("airtemp_ts", metadata$id[i], ".csv")), row.names = FALSE)
  
  print(c(i, dim(x)[1]))
  
}

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))

# Load metadata and admin unit geometries
metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
admin2 <- readRDS(here("data", "admin_units", "simplified_admin2.rds"))
admin1 <- readRDS(here("data", "admin_units", "simplified_admin1.rds"))
admin0 <- readRDS(here("data", "admin_units", "simplified_admin0.rds"))
cluster <- readRDS(here("data/systematic_review_results/cluster_membership.rds")) %>%
  dplyr::select(id, cluster)
metadata <- metadata %>%
  left_join(cluster, by = "id")
metadata$cluster <- as.factor(metadata$cluster)
  
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
    centroid$city <- metadata$city[i]
    centroid$cluster <- metadata$cluster[i]
  } else {
    metadata_admin <- metadata$admin1[i]
    admin_1_data <- admin1[admin1$NAME_1 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_1_data)
    centroid <- st_centroid(admin_geometry)
    centroid <- centroid[c("continent", "ISO", "NAME_0", "geometry", "NAME_1")]
    centroid$NAME_2 <- NA
    centroid$city <- metadata$city[i]
    centroid$cluster <- metadata$cluster[i]
  }
  
  if (i == 1) {
    centroid_metadata <- centroid
  } else {
    centroid_metadata <- rbind(centroid_metadata, centroid)
  }
  
}

# Plotting the centroids of study locations on a background of the below countries
countries <- c("Afghanistan", "Bangladesh", "Bhutan", "Djibouti", "Egypt", "Ethiopia",
               "India", "Iran", "Iraq", "Israel", "Jordan", "Kenya", "Lebanon", "Myanmar", 
               "Nepal", "Oman", "Pakistan", "Saudi Arabia", "Sudan", "Somalia", "South Sudan", "Yemen",
               "China", "Thailand", "Laos", "Cambodia", "Vietnam", "Syria",
               "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Mongolia", "Armenia", "Azerbaijan",
               "Georgia", "Armenia", "Kazakhstan", "Turkey", "Ukraine", "Romania", "Bulgaria", "Libya",
               "Democratic Republic of the Congo", "Central African Republic", "Rwanda", "Burundi", "Tanzania", "Uganda", "Chad",
               "Malaysia", "Singapore", "Indonesia", "Uzbekistan", "Russia")
pres_countries <- c("Afghanistan", "Djibouti", "India", "Iran", "Myanmar", "Pakistan")

red_admin0 <- admin0[admin0$NAME_0 %in% countries, ]
red_admin0
red_admin0$stephensi <- ifelse(red_admin0$NAME_0 %in% pres_countries, "Yes", "No")

centroids <- st_centroid(red_admin0)
points <- as.data.frame(st_coordinates(centroids))
points$NAME_0 <- red_admin0$NAME_0
points <- points[admin0$NAME_0 %in% pres_countries, ]

set.seed(10)
N <- unname(table(metadata$country))
a <- ggplot() + 
  geom_sf(data = red_admin0, fill = ifelse(red_admin0$stephensi == "Yes", gray(0.98), gray(0.85)),
          col = ifelse(red_admin0$stephensi == "Yes", "black", gray(0.6))) +
  geom_sf(data = red_admin0[red_admin0$stephensi == "Yes", ], fill = gray(0.98),
          col = "black") +
  geom_sf(data = st_jitter(centroid_metadata, factor = .008), 
          aes(fill = NAME_0), col = "black", size = 2, shape = 21) +
  geom_label_repel(data = points, aes(x = X, y=Y, label=NAME_0),
                   color = "black", fontface = "bold", max.overlaps = 2,
                   force = 3, nudge_x = c(0, 0, -11, -1, 5, 9),
                   nudge_y = c(6, 4, -5, 7, 10, 9)) +
  theme_bw() +
  scale_fill_manual(values = gg_color_hue(6),
                    labels = c(paste0("Afghanistan (n = ", N[1],")"), 
                               paste0("Djibouti (n = ", N[2],")"),  
                               paste0("India (n = ", N[3],")"), 
                               paste0("Iran (n = ", N[4],")"), 
                               paste0("Myanmar (n = ", N[5],")"), 
                               paste0("Pakistan (n = ", N[6],")"))) +
  coord_sf(ylim = c(0, 45), xlim = c(30, 105), expand = FALSE)  + 
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold")) +
  guides(fill = guide_legend(title = "Country", override.aes = list(size = 5)))

urban_rural_plot <- ggplot() + 
  geom_sf(data = red_admin0, 
          fill = ifelse(red_admin0$stephensi == "Yes", gray(0.98), gray(0.85)),
          col = ifelse(red_admin0$stephensi == "Yes", "black", gray(0.6))) +
  geom_sf(data = red_admin0[red_admin0$stephensi == "Yes", ], fill = gray(0.98),
          col = "black") +
  geom_sf(data = st_jitter(centroid_metadata, factor = .008), 
          aes(fill = interaction(city, cluster)), col = "black", size = 2, shape = 21) +
  # geom_sf(data = st_jitter(centroid_metadata[centroid_metadata$city == "Mixture/Unclear", ], factor = .008), 
  #         fill = "black", col = "black", size = 2, shape = 21) +
  coord_sf(ylim = c(0, 45), xlim = c(30, 105), expand = FALSE) +
  #scale_fill_manual(values = c("#34DFF1", "#289AC1", "#FF8682", "#E2302A")) +
  scale_fill_discrete(labels = c("Unclear 1 Peak", "Rural 1 Peak", "Urban 1 Peak", "Rural 2 Peaks", "Urban 2 Peaks")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Type of Setting and Number of Peaks", nrow = 2))
ggsave2(file = here("figures", "Supp_Figure_UrbanLocations.pdf"), plot = urban_rural_plot, width = 14, height = 7, dpi = 500)

# Plotting some representative time-series
for_plot_index <- c(4, 57, 26, 53, 34, 31)
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
lower_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
upper_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
counter <- 1
for (i in 1:length(for_plot_index)) {
  
  index <- metadata$id[for_plot_index[i]]
  
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
  f_lower <- apply(f, 2, quantile, 0.125)
  f_upper <- apply(f, 2, quantile, 0.875)
  mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  lower <- as.numeric(exp(f_lower)[order(all_timepoints)])
  upper <- as.numeric(exp(f_upper)[order(all_timepoints)])
  mean_realisation[counter, ] <- mean
  lower_realisation[counter, ] <- lower
  upper_realisation[counter, ] <- upper
  counter <- counter + 1
  plot(mean_realisation[i, ], main = paste0(metadata$id[for_plot_index[i]], "  ", metadata$country[for_plot_index[i]]))
  #browser()
}

# Plotting all of the outputs to see which to feature in Fig 1A plot
mean <- as.data.frame(mean_realisation)
colnames(mean) <- paste0("V", seq(0.5, 12.5, 0.5))
mean$country <- metadata$country[for_plot_index]
mean_pv <- mean %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "mean")

lower <- as.data.frame(lower_realisation)
colnames(lower) <- paste0("V", seq(0.5, 12.5, 0.5))
lower$country <- metadata$country[for_plot_index]
lower_pv <- lower %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "lower")

upper <- as.data.frame(upper_realisation)
colnames(upper) <- paste0("V", seq(0.5, 12.5, 0.5))
upper$country <- metadata$country[for_plot_index]
upper_pv <- upper %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "upper")

catch_data <- metadata[for_plot_index, ] %>%
  dplyr::select(country, Jan:Dec)
total_raw_catch <- catch_data %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
  group_by(country) %>%
  summarise(total_raw = sum(raw_catch))
total_per_country <- mean %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
  group_by(country) %>%
  summarise(total = sum(catch))
colnames(catch_data) <- c("country", seq(1, 12, length.out = 12)) 
raw_catch_data <- catch_data %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
  left_join(total_raw_catch, by = "country") %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  dplyr::select(country, timepoint, raw_catch)

overall <- mean_pv %>%
  left_join(lower_pv, by = c("country", "timepoint")) %>%
  left_join(upper_pv, by = c("country", "timepoint")) %>%
  left_join(total_per_country, by = c("country")) %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  left_join(raw_catch_data, by =  c("country", "timepoint"))

# Extracting Rainfall Data (also sort out leap year stuff)
months_length <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
monthly_rainfall_storage <- matrix(nrow = length(for_plot_index), ncol = length(months_length)) 
for (i in 1:length(metadata$id[for_plot_index])) {
  index <- metadata$id[for_plot_index[i]]
  print(index)
  temp <- c()
  rf <- read.csv(here(paste0("data/location_specific_rainfall/rainfall_ts", index, ".csv")))
  rf <- rf %>%
    dplyr::group_by(daymonth_id) %>%
    dplyr::summarise(rainfall = mean(rainfall))
  counter <- 1
  count_vec <- counter
  for (j in 1:length(months_length)) {
    indices <- counter:(counter + months_length[j] - 1)
    temp <- c(temp, sum(rf$rainfall[indices]))
    counter <- counter + months_length[j]
  }
  monthly_rainfall_storage[i, ] <- temp
}

colnames(monthly_rainfall_storage) <- paste0("X", seq(1, 12, 1))
monthly_rainfall_storage <- data.frame(id = for_plot_index, 
                                       country = metadata$country[for_plot_index],
                                       monthly_rainfall_storage)
monthly_rainfall_storage <- monthly_rainfall_storage %>%
  pivot_longer(cols = -c(id, country), names_to = "timepoint", values_to = "rainfall") %>%
  mutate(timepoint = as.numeric(gsub("X", "", timepoint)))
  
overall_ft_rain <- overall %>%
  left_join(monthly_rainfall_storage, by = c("country", "timepoint"))

scaling_factor <- overall_ft_rain %>%
  group_by(country) %>%
  summarise(max_vector = max(upper, na.rm = TRUE),
            max_rainfall = max(rainfall, na.rm = TRUE),
            sf = max_vector/max_rainfall, na.rm = TRUE)
overall_ft_rain_sf <- overall_ft_rain %>%
  left_join(scaling_factor, by = c("country"))
overall_ft_rain_sf_filtered <- overall_ft_rain_sf %>%
  filter(!is.na(rainfall))

b <- ggplot() +
  geom_bar(data = overall_ft_rain_sf_filtered, aes(x = timepoint, y = rainfall * sf), stat = "identity", 
           fill = "#DAFBFF", col = "light grey") +
  geom_path(data = overall_ft_rain_sf, aes(x = timepoint, y = mean, col = country), size = 2) +
  geom_point(data = overall_ft_rain_sf, aes(x = timepoint, y = raw_catch)) +
  geom_ribbon(data = overall_ft_rain_sf, aes(x = timepoint, ymin = lower, ymax = upper, fill = country), alpha = 0.2) +
  facet_wrap(~country, scales = "free_y") +
  scale_y_continuous(limits=c(0, NA), position = "right") +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(1, 12, length.out = 12)) +
  ylab("Monthly Catch") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold")) 

c <- plot_grid(a, b, ncol = 2, align = "h", axis = "b")
ggsave2(file = here("figures", "Fig1_Overall.pdf"), plot = c, width = 15, height = 5.5, dpi = 500)

