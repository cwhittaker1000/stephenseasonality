# See https://github.com/r-spatial/rgee for more details about installation and getting things working

# Load required libraries
library(tidyverse); library(rgee); library(sf); library(raster); library(reticulate); library(here)

# Initialise rgee
ee_Initialize()

# Load metadata and admin 2 units
metadata <- readRDS(here("data", "processed", "metadata_and_processed_counts.rds"))
admin2 <- readRDS(here("data", "processed", "simplified_admin2.rds"))
admin1 <- readRDS(here("data", "processed", "simplified_admin1.rds"))

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

  write.csv(x, file = here("data", "processed", "location_specific_rainfall", paste0("rainfall_ts", metadata$id[i], ".csv")), row.names = FALSE)
    
  print(c(i, dim(x)[1]))
  
}


               