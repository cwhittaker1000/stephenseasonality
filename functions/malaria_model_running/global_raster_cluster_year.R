
global_raster_cluster_year <- function(ISO, what_do, year = "2020"){
  
  #Load shapefiles
  message("Load shapefiles")
  SSA <- readOGR("data/shp/SSA_shapefile_adm0.shp")
  total_countries <- readOGR("data/shp/gadm36_0.shp")
  wanted_countries <- read.csv("data/ISO_use.csv", stringsAsFactors = FALSE)[, 1]
  
  #Sort out extraction - fastest if we do it by country
  sapply(strsplit(year, ";")[[1]], function(x){
    
    message(x)
    
    dummy_data <- raster(list.files("data/modelling_data/population", full.names = T, pattern = "2020"))
    country_dummy <- rasterize(total_countries[which(total_countries$GID_0 == ISO), ], 
                                    crop(dummy_data, extent(total_countries[which(total_countries$GID_0 == ISO), ])), 
                                    mask = T)
    
    template <- raster(extent(country_dummy), crs = crs(country_dummy), resolution = 0.1)
    pts <- as(country_dummy, "SpatialPoints")
    vals <- raster::extract(country_dummy, pts)
    pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
    country_dummy_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
    
    #Now add in data on environment
    if(what_do == "falciparum_prevalence"){
      all_raster_falciparum <- raster(list.files("data/modelling_data/prevalence/PfPR/Raster Data/PfPR_rmean", pattern = paste0("PfPR_", as.character(x)), full.names = T))
      country_falciparum <- rasterize(total_countries[which(total_countries$GID_0 == ISO), ], 
                                      crop(all_raster_falciparum, extent(total_countries[which(total_countries$GID_0 == ISO), ])), 
                                      mask = T)
      
      template <- raster(extent(country_falciparum), crs = crs(country_falciparum), resolution = 0.1)
      pts <- as(country_falciparum, "SpatialPoints")
      vals <- raster::extract(country_falciparum, pts)
      pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
      country_falciparum_resampled <- rasterize(pts2, template, field = "vals", fun = mean, na.rm = T)
      
      writeRaster(country_falciparum_resampled, paste0("output/global_raster/country_level/raster/", ISO, "_falciparumprevalence_", x, ".tif"), overwrite = T)
    }
    
    if(what_do == "Population"){
      message("Population")
      population_data <- raster(list.files("data/modelling_data/population", full.names = T, pattern = x))
      country_population <- rasterize(total_countries[which(total_countries$GID_0 == ISO), ], 
                                      crop(population_data, extent(total_countries[which(total_countries$GID_0 == ISO), ])), 
                                      mask = T)
      
      template <- raster(extent(country_population), crs = crs(country_population), resolution = 0.1)
      pts <- as(country_population, "SpatialPoints")
      vals <- raster::extract(country_population, pts)
      pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
      country_population_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
      
      writeRaster(country_population_resampled, paste0("output/global_raster/country_level/raster/", ISO, "_population_", x, ".tif"), overwrite = T)
    }
    
    if(what_do == "EVI"){
      message("EVI")
      EVI <- sapply(list.files("data/modelling_data/EVI/annual/", full.names = T, pattern = x), function(x) raster(x), simplify = FALSE)[[1]]
      
      country_EVI <- rasterize(total_countries[which(total_countries$GID_0 == ISO), ], 
                               crop(EVI, extent(total_countries[which(total_countries$GID_0 == ISO), ])), 
                               mask = T)
      
      template <- raster(extent(country_EVI), crs = crs(country_EVI), resolution = 0.1)
      pts <- as(country_EVI, "SpatialPoints")
      vals <- raster::extract(country_EVI, pts)
      pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
      country_EVI_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
      
      writeRaster(country_EVI_resampled, paste0("output/global_raster/country_level/raster/", ISO, "_EVI_", x, ".tif"), overwrite = T)
    }
    
    if(what_do == "Landcover"){
      message("Landcover")
      landcover_data <- raster(list.files("data/modelling_data/landcover", full.names = T, pattern = x))
      
      country_LC <- rasterize(total_countries[which(total_countries$GID_0 == ISO), ], 
                              crop(landcover_data, extent(total_countries[which(total_countries$GID_0 == ISO), ])), 
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
                                         country_dummy_resampled,
                                         method = "ngb")
        
        writeRaster(country_LC_resampled, paste0("output/global_raster/country_level/raster/", ISO, "_LC", a, "_", x, ".tif"), overwrite = T)
        
        dz <- data.frame(value = country_LC_resampled[])
        colnames(dz) <- paste0("LC_", a)
        dz
        
      }, simplify = FALSE))
    }
    
  })
  
}