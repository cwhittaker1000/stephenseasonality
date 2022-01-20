
global_raster_cluster <- function(ISO, year = "2020"){
  
  #Load shapefiles
  message("Load shapefiles")
  SSA <- readOGR("data/shp/SSA_shapefile_adm0.shp")
  total_countries <- readOGR("data/shp/gadm36_0.shp")
  wanted_countries <- read.csv("data/ISO_use.csv", stringsAsFactors = FALSE)[, 1]
  
  #Now add in data on environment
  message("Load rasters")
  landcover_data <- raster(list.files("data/modelling_data/landcover", full.names = T, pattern = year))
  population_data <- raster(list.files("data/modelling_data/population", full.names = T, pattern = year))
  worldclim_data <- sapply(list.files("data/modelling_data/worldclim", full.names = T, pattern = ".tif"), function(x) raster(x), simplify = FALSE)
  EVI <- sapply(list.files("data/modelling_data/EVI/annual/", full.names = T, pattern = year), function(x) raster(x), simplify = FALSE)[[1]]
  elevation <- raster("data/modelling_data/landcover/alwdgg.tif")
  cow <- raster("data/modelling_data/livestock/5_Ct_2010_Da.tif")
  goat <- raster("data/modelling_data/livestock/5_Gt_2010_Da.tif")
  sheep <- raster("data/modelling_data/livestock/5_Sh_2010_Da.tif")
  irrigation <- raster("data/modelling_data/gmia_v5_aei_pct.asc")
  
  #Sort out extraction - fastest if we do it by country
  country_data_all <- as.data.frame(rbindlist(sapply(ISO, function(x){
    
    message("Population")
    country_population <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                                    crop(population_data, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                                    mask = T)
    
    template <- raster(extent(country_population), crs = crs(country_population), resolution = 0.1)
    pts <- as(country_population, "SpatialPoints")
    vals <- raster::extract(country_population, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_population_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
    
    writeRaster(country_population_resampled, paste0("output/global_raster/country_level/raster/", x, "_population_", year, ".tif"), overwrite = T)
    
    #Cattle raster
    message("Cow")
    country_cow <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                             crop(cow, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                             mask = T)
    
    template <- raster(extent(country_cow), crs = crs(country_cow), resolution = 0.1)
    pts <- as(country_cow, "SpatialPoints")
    vals <- raster::extract(country_cow, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_cow_resampled <- resample(rasterize(pts2, template, field = "vals", fun = mean, na.rm = T),
                                      country_population_resampled, method = "ngb")
    
    writeRaster(country_cow_resampled, paste0("output/global_raster/country_level/raster/", x, "_cow_", year, ".tif"), overwrite = T)
    
    message("Goat")
    country_goat <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                              crop(goat, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                              mask = T)
    
    template <- raster(extent(country_goat), crs = crs(country_goat), resolution = 0.1)
    pts <- as(country_goat, "SpatialPoints")
    vals <- raster::extract(country_goat, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_goat_resampled <- resample(rasterize(pts2, template, field = "vals", fun = mean, na.rm = T),
                                       country_population_resampled, method = "ngb")    
    
    writeRaster(country_goat_resampled, paste0("output/global_raster/country_level/raster/", x, "_goat_", year, ".tif"), overwrite = T)
    
    message("Sheep")
    country_sheep <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                               crop(sheep, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                               mask = T)
    
    template <- raster(extent(country_sheep), crs = crs(country_sheep), resolution = 0.1)
    pts <- as(country_sheep, "SpatialPoints")
    vals <- raster::extract(country_sheep, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_sheep_resampled <- resample(rasterize(pts2, template, field = "vals", fun = mean, na.rm = T),
                                        country_population_resampled, method = "ngb")    
    writeRaster(country_sheep_resampled, paste0("output/global_raster/country_level/raster/", x, "_sheep_", year, ".tif"), overwrite = T)
    
    message("EVI")
    country_EVI <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                                    crop(EVI, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                                    mask = T)
    
    template <- raster(extent(country_population), crs = crs(country_population), resolution = 0.1)
    pts <- as(country_population, "SpatialPoints")
    vals <- raster::extract(country_population, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_EVI_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
    
    writeRaster(country_EVI_resampled, paste0("output/global_raster/country_level/raster/", x, "_EVI_", year, ".tif"), overwrite = T)

    message("Elevation")
    country_elevation <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                                   crop(elevation, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                                   mask = T)
    
    template <- raster(extent(country_elevation), crs = crs(country_elevation), resolution = 0.1)
    pts <- as(country_elevation, "SpatialPoints")
    vals <- raster::extract(country_elevation, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_elevation_resampled <- resample(rasterize(pts2, template, field = "vals", fun = mean, na.rm = T),
                                            country_population_resampled, method = "ngb")
    
    writeRaster(country_elevation_resampled, paste0("output/global_raster/country_level/raster/", x, "_elevation_", year, ".tif"), overwrite = T)
    
    #Irrigation raster
    message("Irrigation")
    country_irrigation <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                                    crop(irrigation, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                                    mask = T)
    
    template <- raster(extent(country_elevation_resampled), crs = crs(country_elevation_resampled), resolution = 0.1)
    pts <- as(country_irrigation, "SpatialPoints")
    vals <- raster::extract(country_irrigation, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_irrigation_resampled <- rasterize(pts2, template, field = "vals", fun = median, na.rm = T)
    
    writeRaster(country_irrigation_resampled, paste0("output/global_raster/country_level/raster/", x, "_irrigation_", year, ".tif"), overwrite = T)
    
    message("Landcover")
    country_LC <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                            crop(landcover_data, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                            mask = T)
    
    template <- raster(extent(country_LC), crs = crs(country_LC), resolution = 0.1)
    pts <- as(country_LC, "SpatialPoints")
    vals <- raster::extract(country_LC, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    
    #Run for each LC type
    LC_all_done <- do.call(cbind, sapply(unique(unique(pts2$vals)), function(a){
      
      this <- as.data.frame(pts2)
      this[which(this$vals != a), ]$vals <- 0
      this[which(this$vals == a), ]$vals <- 1
      
      country_LC_resampled <- resample(rasterize(SpatialPointsDataFrame(pts2, data.frame(vals = this[, "vals"])), 
                                                        template, field = "vals", fun = mean, na.rm = T),
                                              country_population_resampled,
                                              method = "ngb")
      
      writeRaster(country_LC_resampled, paste0("output/global_raster/country_level/raster/", x, "_LC", a, "_", year, ".tif"), overwrite = T)
      
      dz <- data.frame(value = country_LC_resampled[])
      colnames(dz) <- paste0("LC_", a)
      dz
      
    }, simplify = FALSE))
    
    message("Worldclim")
    country_worldclim_data <- do.call(cbind, sapply(worldclim_data, function(y){
      
      country_worldclim <- rasterize(total_countries[which(total_countries$GID_0 == x), ], 
                                     crop(y, extent(total_countries[which(total_countries$GID_0 == x), ])), 
                                     mask = T)
      
      template <- raster(extent(country_worldclim), crs = crs(country_worldclim), resolution = 0.1)
      pts <- as(country_worldclim, "SpatialPoints")
      vals <- raster::extract(country_worldclim, pts)
      pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
      worldclim_resampled <- resample(rasterize(pts2, template, field = "vals", fun = mean, na.rm = T),
                                      country_population_resampled, method = "ngb")
      
      dfz <- data.frame(val = worldclim_resampled[])
      colnames(dfz) <- paste0("worldclim_", gsub("wc2.1_30s_bio_", "", names(y)))
      writeRaster(worldclim_resampled, paste0("output/global_raster/country_level/raster/", x, "_", colnames(dfz), "_", year, ".tif"), overwrite = T)
      
      dfz
      
    }, simplify = FALSE))
    
    #Plot
    message("Creating dataframe")
    data_go <- cbind(data.frame(ISO = x,
                                cell = 1:ncell(country_population_resampled),
                                elevation = country_elevation_resampled[],
                                population_density = country_population_resampled[],
                                EVI = country_EVI_resampled[],
                                cattle = country_cow_resampled[],
                                goat = country_cow_resampled[],
                                sheep = country_cow_resampled[],
                                irrigation = country_irrigation_resampled[]),
                     LC_all_done,
                     country_worldclim_data
    )
    
    data_go$nrow <- nrow(country_population_resampled)
    data_go$ncol <- ncol(country_population_resampled)
    data_go$lon <- coordinates(country_population_resampled)[, 1]
    data_go$lat <- coordinates(country_population_resampled)[, 2]
    data_go
    
  }, simplify = FALSE), fill = TRUE))
  
  message("Saving dataframe")
  if(year == 2020){
    fwrite(country_data_all[which(!is.na(country_data_all$population_density) & !is.na(country_data_all$elevation)), ], 
           paste0("output/global_raster/country_level/", ISO, "_environment_data.csv"), row.names = FALSE)
  }
  
}
