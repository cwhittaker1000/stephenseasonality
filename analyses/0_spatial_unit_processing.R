#######################################################################################################
##                                                                                                   ##
##                            Loading Required Libraries and Functions                               ##
##                                                                                                   ##
#######################################################################################################
library(dplyr); library(raster); library(rgdal); library(sf); library(raster); library(sp); 
library(tidyverse); library(here)

# Loading in admin 0 units, simplifying and saving
admin0 <- readRDS(here("data", "raw", "raw_gadm_shapefiles", "level0.rds"))
admin0 <- admin0[!is.na(admin0$NAME_0), ]
countries <- c("Afghanistan", "Bangladesh", "Bhutan", "Djibouti", "Egypt", "Ethiopia", "India", "Israel", "Iraq", "Iran", "Jordan", "Kenya", 
               "Lebanon", "Myanmar", "Nepal", "Oman", "Pakistan", "Saudi Arabia", "Somalia", "South Sudan", "Sudan", "Yemen")
admin0 <- admin0[admin0$NAME_0 %in% countries, ]
saveRDS(admin0, here("data", "processed", "complex_admin0.rds"))
for (i in 1:nrow(admin0)) {
  temp <- admin0$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 5000)
  #print(c(object.size(temp), object.size(reduce)))
  #plot(temp)
  #plot(reduce, add = TRUE, type = "l", col = "red")
  admin0$geometry[i] <- reduce
  print(i)
  #browser()
}
saveRDS(admin0, here("data", "processed", "simplified_admin0.rds"))

# Loading in admin 1 units, simplifying and saving
admin1 <- readRDS(here("data", "raw", "raw_gadm_shapefiles", "level1.rds"))
admin1 <- admin1[!is.na(admin1$NAME_0), ]
countries <- c("Afghanistan", "Djibouti", "Ethiopia", "India",  "Iran", "Myanmar", "Pakistan", "Saudi Arabia")
admin1 <- admin1[admin1$NAME_0 %in% countries, ]
saveRDS(admin1, here("data", "processed", "complex_admin1.rds"))
for (i in 1:nrow(admin1)) {
  temp <- admin1$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 3500)
  #print(c(object.size(temp), object.size(reduce)))
  #plot(temp)
  #plot(reduce, add = TRUE, type = "l", col = "red")
  admin1$geometry[i] <- reduce
  print(i)
  #browser()
}
saveRDS(admin1, here("data", "processed", "simplified_admin1.rds"))


# Loading in admin 2 units, simplifying and saving
admin2 <- readRDS(here("data", "raw", "raw_gadm_shapefiles/", "level2.rds"))
admin2 <- admin2[!is.na(admin2$NAME_0), ]
countries <- c("Afghanistan", "Djibouti", "Ethiopia", "India",  "Iran", "Myanmar", "Pakistan", "Saudi Arabia")
admin2 <- admin2[admin2$NAME_0 %in% countries, ]
saveRDS(admin1, here("data", "processed", "complex_admin2.rds"))
#samples <- sample(seq(1:1447), 100)
#samples <- samples[order(samples)]
for (i in 1:nrow(admin2)) {
  temp <- admin2$geometry[i]
  reduce <- st_simplify(temp, dTolerance = 1000)
  #if (i %in% samples) {
  #  print(c(object.size(temp), object.size(reduce)))
  #  plot(temp)
  #  plot(reduce, add = TRUE, type = "l", col = "red")
  #  browser()
  #}
  admin2$geometry[i] <- reduce
  print(i)
}
saveRDS(admin1, here("data", "processed", "simplified_admin2.rds"))

# Checking everything has worked okay
complex_admin_0 <- readRDS(here("data", "processed", "complex_admin0.rds"))
simple_admin_0 <- readRDS(here("data", "processed", "simplified_admin0.rds"))
plot(complex_admin_0$geometry[10])
plot(simple_admin_0$geometry[10], add = TRUE, type = "l", col = "red")

complex_admin_1 <- readRDS(here("data", "processed", "complex_admin1.rds"))
simple_admin_1 <- readRDS(here("data", "processed", "simplified_admin1.rds"))
plot(complex_admin_1$geometry[10])
plot(simple_admin_1$geometry[10], add = TRUE, type = "l", col = "red")

complex_admin_2 <- readRDS(here("data", "processed", "complex_admin2.rds"))
simple_admin_2 <- readRDS(here("data", "processed", "simplified_admin2.rds"))
plot(complex_admin_2$geometry[10])
plot(simple_admin_2$geometry[10], add = TRUE, type = "l", col = "red")

