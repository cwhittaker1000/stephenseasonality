#######################################################################################################
##                                                                                                   ##
##                            Loading Required Libraries and Functions                               ##
##                                                                                                   ##
#######################################################################################################
library(dplyr); library(raster); library(rgdal); library(sf); library(raster); library(sp); 
library(tidyverse); library(here)

# Loading in admin 0 units, simplifying and saving
admin0 <- readRDS(here("data", "raw_gadm_shapefiles", "level0.rds"))
admin0 <- admin0[!is.na(admin0$NAME_0), ]
countries <- c("Afghanistan", "Armenia", "Azerbaijan", "Bangladesh", "Bhutan", "Cambodia", "China", "Djibouti", "Egypt", 
               "Ethiopia", "India", "Israel", "Iraq", "Iran", "Jordan", "Kenya", "Kyrgyzstan",
               "Laos", "Lebanon", "Myanmar", "Nepal", "Oman", "Pakistan", "Saudi Arabia", 
               "Somalia", "South Sudan", "Sudan", "Syria", "Tajikistan", "Thailand", "Turkey", "Turkmenistan", "Uzbekistan",
               "Vietnam", "Yemen")
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

# Checking everything has worked okay
# complex_admin_0 <- readRDS(here("data", "admin_units", "complex_admin0.rds"))
# simple_admin_0 <- readRDS(here("data", "admin_units", "simplified_admin0.rds"))
# plot(complex_admin_0$geometry[10])
# plot(simple_admin_0$geometry[10], add = TRUE, type = "l", col = "red")
# 
# complex_admin_1 <- readRDS(here("data", "adm_units", "complex_admin1.rds"))
# simple_admin_1 <- readRDS(here("data", "admin_units", "simplified_admin1.rds"))
# plot(complex_admin_1$geometry[10])
# plot(simple_admin_1$geometry[10], add = TRUE, type = "l", col = "red")
# 
# complex_admin_2 <- readRDS(here("data", "admin_units", "complex_admin2.rds"))
# simple_admin_2 <- readRDS(here("data", "admin_units", "simplified_admin2.rds"))
# plot(complex_admin_2$geometry[10])
# plot(simple_admin_2$geometry[10], add = TRUE, type = "l", col = "red")