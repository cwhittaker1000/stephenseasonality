# See https://github.com/r-spatial/rgee for more details about installation and getting things working

# Load required libraries
library(tidyverse); library(rgee); library(sf); library(raster); library(reticulate)

# Initialise rgee
ee_Initialize()

adm1data <- getData('GADM', country = 'IND', level = 1)
example <- adm1data[11, ]

test <- st_as_sf(example)

nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

terraclimate <- ee$ImageCollection("UCSB-CHG/CHIRPS/DAILY") %>%
  ee$ImageCollection$filterDate("2001-01-01", "2002-01-01") %>%
  ee$ImageCollection$map(function(x) x$select("precipitation")) %>% # Select only precipitation bands
  ee$ImageCollection$toBands() # from imagecollection to image

ee_nc_rain <- ee_extract(x = terraclimate, y = nc[1, "NAME"], sf = FALSE)
ee_nc_rain <- ee_extract(x = terraclimate, y = test, sf = FALSE)

x <- ee_nc_rain %>%
  pivot_longer(-NAME, names_to = "day", values_to = "pr")

x <- ee_nc_rain %>%
  pivot_longer(X20010101_precipitation:X20011231_precipitation, names_to = "day", values_to = "pr")


x$day <- gsub("X", "", x$day)
x$day <- gsub("_precipitation", "", x$day)
x$day <- as.Date(x$day, format = "%Y%m%d")

ggplot(x, aes(x = day, y = pr, color = pr)) +
  geom_line(alpha = 0.4) +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme_minimal()



               