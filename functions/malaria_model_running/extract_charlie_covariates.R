
extract_charlie_covariates <- function(Country){
  
  #Load in Charlie locations
  charlie_locations <- as.data.frame(read_excel("data/mosquito/stephensi_locations.xlsx"))
  charlie_locations_subset <- charlie_locations[which(charlie_locations$Country %in% if(grepl(";", Country)) unlist(strsplit(Country, ";")) else Country), ]
  
  #Load shapefiles
  #Saudi has no adm2, so have to load in a different shapefile
  if(Country == "Saudi Arabia"){
    total_countries <- readOGR("data/shp/gadm36_SAU_1.shp")
    total_countries$GID_2 <- NA
  } else {
    globe <- readOGR("data/shp/gadm36_2.shp")
    total_countries <- subset(globe, NAME_0 %in% unique(charlie_locations$Country))
  }
  
  #Now add in data on environment
  landcover_data <- sapply(list.files("data/modelling_data/landcover", full.names = T, "LC-L4"), function(x) raster(x), simplify = FALSE)
  population_data <- sapply(list.files("data/modelling_data/population", full.names = T, ".tif"), function(x) raster(x), simplify = FALSE)
  worldclim_data <- sapply(list.files("data/modelling_data/worldclim", full.names = T, ".tif"), function(x) raster(x), simplify = FALSE)
  EVI <- sapply(list.files("data/modelling_data/EVI/annual/", full.names = T, pattern = ".tif"), function(x) raster(x), simplify = FALSE)
  elevation <- raster("data/modelling_data/landcover/alwdgg.tif")
  
  all_runs_run <- rbindlist(sapply(1:nrow(charlie_locations_subset), function(x){
    
    message(paste0("Row ", x, " of ", nrow(charlie_locations_subset)))
    
    #Find the admin unit
    this_row <- charlie_locations_subset[x, ]
    this_unit <- if(is.na(this_row$`Admin 2`)| this_row$`Admin 2` == "NA"){
      subset(total_countries, NAME_0 == this_row$Country & NAME_1 == this_row$`Admin 1`)
    } else {
      subset(total_countries, NAME_0 == this_row$Country & NAME_1 == this_row$`Admin 1` & NAME_2 == this_row$`Admin 2`)
    }
    
    #Extract year
    years_use <- this_row$`Year Start`:this_row$`Year End`
    
    #Run for each year
    all_year_run <- rbindlist(sapply(years_use, function(y){
      
      message(y)
      
      #Now we extract the data - these do not need to be looped because its a single value we want
      landcover_data_here <- rasterize(this_unit, 
                                       crop(landcover_data[[which(grepl(pmax(y, 1992), names(landcover_data)))]], extent(this_unit)), 
                                       mask = T)
      #Landcover we also want to look at the proportion of cells that are of each LC type
      LC_prop <- rbind(table(landcover_data_here[])/sum(table(landcover_data_here[]), na.rm = T))
      colnames(LC_prop) <- paste0("LC_", colnames(LC_prop))
      
      population_data_here <- rasterize(this_unit, 
                                        crop(population_data[[which(grepl(pmax(y, 2000), names(population_data)))]], extent(this_unit)), 
                                        mask = T)
      EVI_data_here <- rasterize(this_unit, 
                                 crop(EVI[[which(grepl(pmax(y, 2001), names(EVI)))]], extent(this_unit)), 
                                 mask = T)
      elevation_data_here <- rasterize(this_unit, 
                                       crop(elevation, extent(this_unit)), 
                                       mask = T)
      
      #This has to be looped because we have multiple difference worldclims we want data from
      worldclim_data_here <- do.call(cbind, sapply(worldclim_data, function(z){
        df <- data.frame(rasterize(this_unit, 
                                   crop(z, extent(this_unit)), 
                                   mask = T)[],
                         stringsAsFactors = FALSE)
        colnames(df) <- gsub("wc2.1_30s_bio_", "worldclim_", names(z))
        df
      }, simplify = FALSE))
      
      #Now we combine everything into a dataframe
      data.frame(this_row[, c("Time Series ID", "Country", "Admin 1", "Admin 2", "City?")],
                 year = y,
                 as.data.frame(this_unit)[, c("GID_1", "GID_2")],
                 population_per_1km = sum(population_data_here[], na.rm = T)/(area(this_unit)/1000^2), #This calculate is needed to work out the population per km2
                 EVI = median(EVI_data_here[], na.rm = T), #Median EVI
                 elevation = median(elevation_data_here[], na.rm = T), #Median elevation
                 LC_prop, # The preprepared landcover and worldclim data
                 rbind(apply(worldclim_data_here, 2, median, na.rm = T)), #This is to work out the median value per worldclim
                 stringsAsFactors = FALSE)
      
    }, simplify = FALSE), fill = T)
    
    all_year_run
    
  }, simplify = FALSE), fill = T)
  
  write.csv(all_runs_run, paste0("output/charlie_data/", Country, "_data.csv"), row.names = FALSE)
  
}