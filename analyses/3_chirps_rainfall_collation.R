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
for (i in 45:nrow(metadata)) {
  
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

# Generating monthly rainfall
monthly_sum_rainfall_storage <- matrix(nrow = nrow(metadata), ncol = 12)
rainfall_files <- list.files(here("data", "processed", "location_specific_rainfall"))
for (i in 1:nrow(metadata)) {
  temp_rainfall <- read.csv(here("data", "processed", "location_specific_rainfall", rainfall_files[i]), stringsAsFactors = FALSE)
  temp_rainfall$month <- format(as.Date(temp_rainfall$day), "%m")
  temp_rainfall <- temp_rainfall %>%
    group_by(month) %>%
    summarise(rainfall = sum(rainfall, na.rm = TRUE))
  monthly_sum_rainfall_storage[i, ] <- temp_rainfall$rainfall
  print(i)
}
monthly_sum_rainfall_storage <- cbind(metadata$id, monthly_sum_rainfall_storage)
colnames(monthly_sum_rainfall_storage) <- c("id", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
monthly_sum_rainfall_storage <- as.data.frame(monthly_sum_rainfall_storage) %>%
  tidyr::pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "rainfall")

ento_data <- metadata %>%
  dplyr::select(id, city, Jan:Dec) %>%
  pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "catch")

overall <- ento_data %>%
  left_join(monthly_sum_rainfall_storage, by = c("id", "month")) 

months <- overall$month

overall <- overall %>%
  group_by(id) %>%
  summarise(catch = catch/sum(catch),
            rainfall = rainfall/sum(rainfall))
overall$month <- months
overall$month <- factor(overall$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

par(mfrow = c(12, 6))
ggplot(data = overall[overall$id == 1, ]) +
  geom_point(aes(x = month, y = catch), colour = "#E0521A", size = 2) +
  geom_line(aes(x = month, y = catch, group = 1)) +
  geom_bar(aes(x = month, y = rainfall), stat = "identity", fill = "#6EB4D1", alpha = 0.5) +
  facet_wrap(~id, scales = "free") +
  scale_x_discrete(name = unique(new_overall$name),
                   labels = c("J", "", "", "F", "", "", "M", "", "", 
                              "A", "", "", "M", "", "", "J", "", "",
                              "J", "", "", "A", "", "", "S", "", "",
                              "O", "", "", "N", "", "", "D", "", ""))




mean_realisation <- matrix(nrow = nrow(metadata), ncol = (12 * interpolating_points + 1))
features <- matrix(nrow = length(retain_index), ncol = 7)
colnames(features) <- c("entropy", "period", "prop_points", "jan_dist", "peaks", "mean", "weight")
for (i in 1:length(retain_index)) {
  
  index <- retain_index[i]
  
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


ggplot(x, aes(x = day, y = rainfall, color = rainfall)) +
  geom_line(alpha = 0.4) +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme_minimal()

               