
create_abundance_environment_df_cluster <- function(id_of_presence = NA, id_of_pseudoabsence = NA){
  
  presence <- read.csv("data/processed_data/modelling_df/raw_data_going_into_presence.csv", stringsAsFactors = FALSE)
  absence <- read.csv("data/processed_data/modelling_df/raw_data_going_into_pseudoabsence.csv", stringsAsFactors = FALSE)
  
  run_these_presence <- if(grepl("@", id_of_presence)) strsplit(id_of_presence, "@")[[1]] else id_of_presence
  run_these_absence <- if(grepl("@", id_of_pseudoabsence)) strsplit(id_of_pseudoabsence, "@")[[1]] else id_of_pseudoabsence
  
  what_number_unique_presence <- which(unique(presence$unique_setting_ID) %in% run_these_presence)
  what_number_unique_pseudoabsence <- which(unique(absence$site_id) %in% run_these_absence)
  
  #Raster of the correct crs, chosen BHR because its tiny
  population_data <- raster("output/global_raster/country_level/raster/BHR_population_2020.tif")
  
  if(!is.na(id_of_presence)){
    malaria_environment_species <- as.data.frame(rbindlist(sapply(run_these_presence, function(x){
      
      message(paste0(which(unique(presence$unique_setting_ID) == x), " of ", length(unique(presence$unique_setting_ID))))
      
      these_locations <- presence[which(presence$unique_setting_ID == x), ]
      these_coord <- matrix(c(as.numeric(unique(strsplit(these_locations$lon, ";")[[1]])),
                              as.numeric(unique(strsplit(these_locations$lat, ";")[[1]]))), ncol = 2)
      
      coord_spdf <- SpatialPoints(these_coord)
      crs(coord_spdf) <- crs(population_data)
      
      coord_buffed <- gBuffer(coord_spdf, width = 0.05, byid = FALSE,
                              capStyle = "SQUARE")
      
      rast_loc <- "output/global_raster/country_level/raster"
      all_rast_files <- list.files(rast_loc, 
                                   pattern = paste0(paste(unique(these_locations$ISO), collapse = "|"), "_"), #expanded search than just ISO to make sure it doesnt pick up any combination of letters that may also be an variable
                                   full.names = T)
      
      LC_all <- all_rast_files[grepl("LC", all_rast_files)]
      population_all <- all_rast_files[grepl("population", all_rast_files)]
      EVI_all <- all_rast_files[grepl("EVI", all_rast_files)]
      
      worldclim <- all_rast_files[grepl("worldclim", all_rast_files)]
      elevation <- raster(all_rast_files[grepl("elevation", all_rast_files)])
      cattle <- raster(all_rast_files[grepl("cow", all_rast_files)])
      goat <- raster(all_rast_files[grepl("goat", all_rast_files)])
      sheep <- raster(all_rast_files[grepl("sheep", all_rast_files)])
      irrigation <- raster(all_rast_files[grepl("irrigation", all_rast_files)])
      
      go_through_all_years <- as.data.frame(rbindlist(sapply(unlist(unique(strsplit(as.character(these_locations$year), ";"))), function(y){
        
        this_LC <- LC_all[which(grepl(min(max(c(y, 1992)), 2019), (LC_all)))]
        this_population <- raster(population_all[[which(grepl(min(max(c(y, 2000)), 2020), (population_all)))]])
        this_EVI <- raster(EVI_all[[which(grepl(min(max(c(y, 2001)), 2020), (EVI_all)))]])
        
        # loc_extent_wide <- (gBuffer(coord_spdf, width = 0.2, byid = FALSE, capStyle = "SQUARE"))
        # LC_present <- rasterize(loc_extent_wide, crop(this_LC, loc_extent_wide), mask = T)
        # country_use <- globe[which(globe$GID_0 %in% these_locations$ISO), ]
        # LC_country_value <- rasterize(country_use, crop(this_LC, extent(country_use)), mask = T)
        # 
        # LC_all_table <- as.data.frame(rbind(t(table(rasterize(gBuffer(coord_spdf, width = 0.05, byid = FALSE,
        #                                                               capStyle = "SQUARE"), 
        #                                                       crop(this_LC, extent(gBuffer(coord_spdf, width = 0.05, byid = FALSE,
        #                                                                                    capStyle = "SQUARE"))), mask = T)[], 
        #                                             useNA = "no"))))
        # 
        # colnames(LC_all_table) <- paste0("LC_", colnames(LC_all_table))
        # LC_all_table <- LC_all_table/sum(LC_all_table)
        
        
        LC_all_2table <- sapply(this_LC, function(x){
          loaded_raster <- raster(x)
          df_go <- data.frame(value = rasterize(coord_buffed, 
                                                crop(loaded_raster, extent(coord_buffed)), mask = T)[])
          colnames(df_go) <- strsplit(x, "_")[[1]][(length(strsplit(x, "_")[[1]])-1)]
          df_go
          
        }, simplify = FALSE)
        nonNALC <- colSums(do.call(cbind, LC_all_2table), na.rm = T)/sum(colSums(do.call(cbind, LC_all_2table), na.rm = T))
        LC_all_table <- nonNALC[-which(nonNALC == 0)]
        
        #Population all
        pop_all <- sum(as.numeric(na.omit(rasterize(coord_buffed, 
                                                    crop(this_population, extent(coord_buffed)), mask = T)[])))
        
        #Elevation all
        elevation_all <- median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                             crop(elevation, extent(coord_buffed)), mask = T)[]))
        )
        
        #EVI all
        EVI_all <- median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                       crop(this_EVI, extent(coord_buffed)), mask = T)[]))
        )
        
        #Worldclim all
        worldclim_all <- do.call(cbind, sapply(worldclim, function(x){
          raster_this <- raster(x)
          worldclimvalue <- data.frame(value = median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                                                   crop(raster_this, extent(coord_buffed)), mask = T)[]))))
          name <- strsplit(x, "_")[[1]]
          colnames(worldclimvalue) <- paste(name[(length(name)-2):(length(name)-1)], collapse = "_")
          worldclimvalue
        }, simplify = FALSE))
        
        #Cattle
        cattle_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                   crop(cattle, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #goat
        goat_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                 crop(goat, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #sheep
        sheep_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                  crop(sheep, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #Irrigation
        irrigation_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                       crop(irrigation, extent(coord_buffed)), mask = T)[]))
        
        
        cbind(these_locations[which(these_locations$species == "stephensi" & grepl(y, these_locations$year)), ], 
              data.frame(studyID = x,
                         ISO = unique(these_locations$ISO),
                         year = y,
                         population_density = sum(pop_all),
                         elevation = median(elevation_all),
                         EVI = EVI_all,
                         cattle = sum(cattle_all),
                         goat = sum(goat_all),
                         sheep = sum(sheep_all),
                         irrigation = median(irrigation_all),
                         worldclim_all,
                         as.data.frame(t(LC_all_table))
              ))
        
      }, simplify = FALSE), fill = T))
      
      go_through_all_years
      
    }, simplify = FALSE), fill = T))
    
    if(nrow(malaria_environment_species) != 0) malaria_environment_species$stephensi_present <- 1
  }
  
  if(!is.na(id_of_pseudoabsence)){
    no_stephensi_environment <- as.data.frame(rbindlist(sapply(as.character(run_these_absence), function(x){
      
      message(paste0(which(unique(absence$site_id) == x), " of ", length(unique(absence$site_id))))
      
      these_locations <- absence[which(absence$site_id == x), ]
      these_coord <- matrix(c(as.numeric(unique(strsplit(as.character(these_locations$longitude), ";")[[1]])),
                              as.numeric(unique(strsplit(as.character(these_locations$latitude), ";")[[1]]))), ncol = 2)
      
      coord_spdf <- SpatialPoints(these_coord)
      crs(coord_spdf) <- crs(population_data)
      
      rast_loc <- "output/global_raster/country_level/raster"
      all_rast_files <- list.files(rast_loc, 
                                   pattern = paste0(paste(unique(these_locations$country_id), collapse = "|"), "_"), #expanded search than just ISO to make sure it doesnt pick up any combination of letters that may also be an variable
                                   full.names = T)
      
      LC_all <- all_rast_files[grepl("LC", all_rast_files)]
      population_all <- all_rast_files[grepl("population", all_rast_files)]
      EVI_all <- all_rast_files[grepl("EVI", all_rast_files)]
      
      worldclim <- all_rast_files[grepl("worldclim", all_rast_files)]
      elevation <- raster(all_rast_files[grepl("elevation", all_rast_files)])
      cattle <- raster(all_rast_files[grepl("cow", all_rast_files)])
      goat <- raster(all_rast_files[grepl("goat", all_rast_files)])
      sheep <- raster(all_rast_files[grepl("sheep", all_rast_files)])
      irrigation <- raster(all_rast_files[grepl("irrigation", all_rast_files)])
      
      go_through_all_years <- as.data.frame(rbindlist(sapply(((as.numeric(na.omit(unique(c(these_locations$year_start, these_locations$year_end)))))), function(y){
        
        this_LC <- LC_all[which(grepl(min(max(c(y, 1992)), 2019), (LC_all)))]
        this_population <- raster(population_all[[which(grepl(min(max(c(y, 2000)), 2020), (population_all)))]])
        this_EVI <- raster(EVI_all[[which(grepl(min(max(c(y, 2001)), 2020), (EVI_all)))]])
        
        coord_buffed <- gBuffer(coord_spdf, width = 0.05, byid = FALSE,
                                capStyle = "SQUARE")
        
        LC_all_2table <- sapply(this_LC, function(x){
          loaded_raster <- raster(x)
          df_go <- data.frame(value = rasterize(coord_buffed, 
                                                crop(loaded_raster, extent(coord_buffed)), mask = T)[])
          colnames(df_go) <- strsplit(x, "_")[[1]][(length(strsplit(x, "_")[[1]])-1)]
          df_go
          
        }, simplify = FALSE)
        nonNALC <- colSums(do.call(cbind, LC_all_2table), na.rm = T)/sum(colSums(do.call(cbind, LC_all_2table), na.rm = T))
        LC_all_table <- nonNALC[-which(nonNALC == 0)]
        
        #Population all
        pop_all <- sum(as.numeric(na.omit(rasterize(coord_buffed, 
                                                    crop(this_population, extent(coord_buffed)), mask = T)[])))
        
        #Elevation all
        elevation_all <- median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                             crop(elevation, extent(coord_buffed)), mask = T)[]))
        )
        
        #EVI all
        EVI_all <- median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                       crop(this_EVI, extent(coord_buffed)), mask = T)[]))
        )
        
        #Worldclim all
        worldclim_all <- do.call(cbind, sapply(worldclim, function(x){
          raster_this <- raster(x)
          worldclimvalue <- data.frame(value = median(as.numeric(na.omit(rasterize(coord_buffed, 
                                                                                   crop(raster_this, extent(coord_buffed)), mask = T)[]))))
          name <- strsplit(x, "_")[[1]]
          colnames(worldclimvalue) <- paste(name[(length(name)-2):(length(name)-1)], collapse = "_")
          worldclimvalue
        }, simplify = FALSE))
        
        #Cattle
        cattle_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                   crop(cattle, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #goat
        goat_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                 crop(goat, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #sheep
        sheep_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                  crop(sheep, extent(gBuffer(coord_spdf, width = 0.05, capStyle = "SQUARE"))), mask = T)[]))
        
        #Irrigation
        irrigation_all <- as.numeric(na.omit(rasterize(coord_buffed, 
                                                       crop(irrigation, extent(coord_buffed)), mask = T)[]))
        
        
        
        data.frame(studyID = x,
                   ISO = unique(these_locations$country_id),
                   lat = unique(these_locations$latitude),
                   lon = unique(these_locations$longitude),
                   year = y,
                   population_density = sum(pop_all),
                   EVI = EVI_all,
                   elevation = elevation_all,
                   cattle = sum(cattle_all),
                   goat = sum(goat_all),
                   sheep = sum(sheep_all),
                   irrigation = median(irrigation_all),
                   worldclim_all,
                   t(LC_all_table)
        )
        
      }, simplify = FALSE), fill = T))
      
      go_through_all_years
      
    }, simplify = FALSE), fill = T))
  }
  
  if(!is.na(id_of_presence)){
    write.csv(malaria_environment_species, file = paste0("data/processed_data/modelling_df/chunk/presence/row_", paste(what_number_unique_presence, collapse = ";"), "_of_", length(unique(presence$site_id)), ".csv"), row.names = FALSE)
  }
  if(!is.na(id_of_pseudoabsence)){
    write.csv(no_stephensi_environment, file = paste0("data/processed_data/modelling_df/chunk/pseudoabsence/row_", paste(what_number_unique_pseudoabsence, collapse = ";"), "_of_", length(unique(absence$site_id)), ".csv"), row.names = FALSE)
  }
  
}


