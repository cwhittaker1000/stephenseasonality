#######################################################################################################
##                                                                                                   ##
##                            Loading Required Libraries and Functions                               ##
##                                                                                                   ##
#######################################################################################################
library(dplyr); library(raster); library(rgdal); library(sf); library(raster); library(sp); 
library(tidyverse); library(here); library(rgee);  library(reticulate)

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
admin1 <- readRDS(here("data", "raw_gadm_shapefiles", "level1.rds"))
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
admin2 <- readRDS(here("data", "raw_gadm_shapefiles/", "level2.rds"))
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
admin1_csv <- data.frame(country = simple_admin_1$NAME_0, admin1 = simple_admin_1$NAME_1, 
                         id = simple_admin_1$GID_1, type = simple_admin_1$TYPE_1, alt_names = simple_admin_1$VARNAME_1)

admin2_csv <- data.frame(country = simple_admin_2$NAME_0, 
                         admin1 = simple_admin_2$NAME_1, id1 = simple_admin_2$GID_1,
                         admin2 = simple_admin_2$NAME_2, id2 = simple_admin_2$GID_2,
                         type = simple_admin_2$TYPE_2, alt_names = simple_admin_2$VARNAME_2)

write.csv(admin1_csv, file = here("data", "admin_units", "admin1_details.csv"), row.names = FALSE)
write.csv(admin2_csv, file = here("data", "admin_units", "admin2_details.csv"), row.names = FALSE)

# Initialise rgee -  See https://github.com/r-spatial/rgee for more details about installation and getting things working
ee_Initialize()
#py_install("geemap")
gm <- import("geemap")
ee_check() # Check non-R dependencies

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
