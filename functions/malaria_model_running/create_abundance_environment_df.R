
create_abundance_environment_df <- function(relative_abundance, redo_pseudoabsence = F, redo_all = F){
  
  #Load in previous to see what is new to be added
  if(file.exists("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")) old_file <- read.csv("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")
  if(redo_all != T){
    if(file.exists("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")){
      relative_abundance <- relative_abundance[-which(relative_abundance$unique_setting_ID %in% old_file$unique_setting_ID), ]
      relative_abundance <- relative_abundance[!duplicated(relative_abundance), ]
      if(nrow(relative_abundance) == 0){
        message("No new data, ending process")
        break()
      }
    }
  }
  
  #Malaria data for each point
  globe <- readOGR("data/shp/gadm36_0.shp")
  locations_use <- subset(globe, GID_0 %in% c("IND", "PAK", "IRN", "MMR", "THA", "SAU", 
                                              "AFG", "NPL", "BGD", "ETH", "DJI", "SDN",
                                              "KHM", "LAO", "TJK", "TKM", "UZB"))
  
  #Load in vector bionomic data
  vector_bionomic_data <- read.csv("data/mosquito/processed_species_bionomics.csv", stringsAsFactors = F)
  vector_bionomic_data <- vector_bionomic_data[-which(is.na(vector_bionomic_data$species)), ]
  
  #Load in contribution
  contribution_loaded <- read.csv("output/vector_competence/malariasimulation_stephensi_density_uncertainty.csv")
  
  # #Now add in data on environment
  landcover_data <- sapply(list.files("data/modelling_data/landcover", pattern = "LCCS", full.names = T), function(x) raster(x), simplify = FALSE)
  # population_data <- sapply(list.files("data/modelling_data/population", pattern = "tif", full.names = T), function(x) raster(x), simplify = FALSE)
  # worldclim <- sapply(list.files("data/modelling_data/worldclim", pattern = "tif", full.names = T), function(x) raster(x), simplify = FALSE)
  # elevation <- raster("data/modelling_data/landcover/alwdgg.tif")
  # EVI <- sapply(list.files("data/modelling_data/EVI/annual/", pattern = "tif", full.names = T), function(x) raster(x), simplify = FALSE)
  # cattle <- raster("data/modelling_data/livestock/5_Ct_2010_Da.tif")
  # goat <- raster("data/modelling_data/livestock/5_Gt_2010_Da.tif")
  # sheep <- raster("data/modelling_data/livestock/5_Sh_2010_Da.tif")
  # irrigation <- raster("data/modelling_data/gmia_v5_aei_pct.asc")
  
  #Where stephensi isnt
  these_ISO <- if(any(locations_use$GID_0 %in% c("ETH", "DJI", "SDN"))){
    locations_use$GID_0[-which(locations_use$GID_0 %in% c("ETH", "DJI", "SDN"))]
  } else {
    locations_use$GID_0
  }
  
  #This is current down because MAP hasnt renewed their security license
  # vector_occurence <- getVecOcc(ISO = these_ISO)
  # vector_occurence <- vector_occurence[which(vector_occurence$year_start >= 1990), ]
  vector_occurence <- read.csv("data/mosquito/Anopheline_Data.csv")
  vector_occurence <- vector_occurence[which(vector_occurence$country_id %in% these_ISO), ]
  
  mosquito_relative_abundance <- as.data.frame(read_excel("data/mosquito/mosquito_relative_abundance.xlsx"))
  no_stephensi_additional <- mosquito_relative_abundance[which(mosquito_relative_abundance$stephensi_absent == 1), ]
  no_stephensi_additional_format <- data.frame(site_id = no_stephensi_additional$ID,
                                               country_id = globe[which(globe$NAME_0 %in% no_stephensi_additional$country), ]$GID_0,
                                               year_start = sapply(strsplit(no_stephensi_additional$year_recorded, ";"), function(x) min(x)),
                                               latitude = no_stephensi_additional$lat,
                                               longitude = no_stephensi_additional$lon,
                                               stringsAsFactors = FALSE)
  
  where_is_stephensi <- vector_occurence[which(vector_occurence$species_plain == "Anopheles stephensi"), ]
  where_is_no_stephensi <- vector_occurence[-which(vector_occurence$site_id %in% where_is_stephensi$site_id), ]
  where_is_no_stephensi <- plyr::rbind.fill(where_is_no_stephensi, no_stephensi_additional_format)
  if(any(where_is_no_stephensi$longitude == "NA") | any(where_is_no_stephensi$latitude == "NA")) where_is_no_stephensi <- where_is_no_stephensi[-which((where_is_no_stephensi$longitude == "NA") | (where_is_no_stephensi$latitude == "NA")), ]
  
  #Here stephensi
  study_relative_contribution <- if(any(is.na(relative_abundance$lat) | is.na(relative_abundance$lon))) relative_abundance[-which(is.na(relative_abundance$lat) | is.na(relative_abundance$lon)), ] else relative_abundance
  
  here_stephensi_sp <- gBuffer(SpatialPoints(coordinates(matrix(c(c(where_is_stephensi$longitude, as.numeric(unlist(strsplit(subset(study_relative_contribution, species == "stephensi")$lon, ";")))),
                                                                  c(where_is_stephensi$latitude, as.numeric(unlist(strsplit(subset(study_relative_contribution, species == "stephensi")$lat, ";"))))), 
                                                                ncol = 2))), 
                               width = 0.1, byid = F)
  
  keep_these <- sapply(unique(where_is_no_stephensi$site_id), function(x){
    # print(x)
    this_site <- where_is_no_stephensi[which(where_is_no_stephensi$site_id == x), ]
    here_no_stephensi_sp <- gBuffer(SpatialPoints(coordinates(matrix(c(as.numeric(this_site$longitude), 
                                                                       as.numeric(this_site$latitude)), ncol = 2))), width = .1, byid = F)
    if(is.null(raster::intersect(here_no_stephensi_sp, here_stephensi_sp))) x else NA
  })
  
  where_is_no_stephensi_true <- where_is_no_stephensi[which(where_is_no_stephensi$site_id %in% keep_these), ]
  
  #Process these
  study_relative_contribution <- study_relative_contribution[which(study_relative_contribution$species == "stephensi"), ]
  
  write.csv(study_relative_contribution, "data/processed_data/modelling_df/raw_data_going_into_presence.csv", row.names = FALSE)
  write.csv(where_is_no_stephensi_true, "data/processed_data/modelling_df/raw_data_going_into_pseudoabsence.csv", row.names = FALSE)
  
  malaria_environment_species <- as.data.frame(rbindlist(sapply(unique(study_relative_contribution$unique_setting_ID), function(x){
    
    message(paste0(which(unique(study_relative_contribution$unique_setting_ID) == x), " of ", length(unique(study_relative_contribution$unique_setting_ID))))
    
    these_locations <- study_relative_contribution[which(study_relative_contribution$unique_setting_ID == x), ]
    these_coord <- matrix(c(as.numeric(unique(strsplit(these_locations$lon, ";")[[1]])),
                            as.numeric(unique(strsplit(these_locations$lat, ";")[[1]]))), ncol = 2)
    
    coord_spdf <- SpatialPoints(these_coord)
    crs(coord_spdf) <- crs(landcover_data[[1]])
    
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
  
  #Now repeat for absence points
  if(redo_pseudoabsence == T){
    no_stephensi_environment <- as.data.frame(rbindlist(sapply(unique(where_is_no_stephensi_true$site_id), function(x){
      
      message(paste0(which(unique(where_is_no_stephensi_true$site_id) == x), " of ", length(unique(where_is_no_stephensi_true$site_id))))
      
      these_locations <- where_is_no_stephensi_true[which(where_is_no_stephensi_true$site_id == x), ]
      these_coord <- matrix(c(as.numeric(unique(strsplit(as.character(these_locations$longitude), ";")[[1]])),
                              as.numeric(unique(strsplit(as.character(these_locations$latitude), ";")[[1]]))), ncol = 2)
      
      coord_spdf <- SpatialPoints(these_coord)
      crs(coord_spdf) <- crs(landcover_data[[1]])
      
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
    no_stephensi_environment$stephensi_present <- 0
    total_data_with_without_stephensi <- plyr::rbind.fill(malaria_environment_species, no_stephensi_environment)
  } else {
    total_data_with_without_stephensi <- malaria_environment_species
  }
  
  LC_value_NA <- as.matrix(total_data_with_without_stephensi[, which(grepl("LC_", colnames(total_data_with_without_stephensi)))])
  LC_value_NA[which(is.na(LC_value_NA))] <- 0
  total_data_with_without_stephensi[, which(grepl("LC_", colnames(total_data_with_without_stephensi)))] <- LC_value_NA
  
  #Save
  if(file.exists("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")){
    total_data_with_without_stephensi <- plyr::rbind.fill(old_file, total_data_with_without_stephensi)
  }
  
  
  write.csv(total_data_with_without_stephensi[!duplicated(total_data_with_without_stephensi), ], 
            "data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv", 
            row.names = FALSE)
  
}
