# covariate_formula = run_these$covariate_formula[1]
# fit_to = run_these$fit_to[1]
# model_type = run_these$model_type[1]
# row_take_samples = run_these$row_take_samples[1]
# ISO = "SDN"
# saveid = run_these$saveid[1]
# seed = 1
# prediction_type = "response"
# degree_take_samples <-  run_these$degree_take_samples[1]
# pseudoratio = run_these$pseudoratio[1]
# criteria_select = run_these$criteria_select[1]
# 


model_predict_fit <- function(covariate_formula, fit_to, model_type,
                              degree_take_samples,
                              row_take_samples, 
                              ISO = "", 
                              saveid, 
                              seed = 1,
                              pseudoratio = 1,
                              predict = T,
                              criteria_select = "RMSE"){
  
  #Set seed and remove scientific notation
  set.seed(seed)
  options(scipen = 999)
  
  message("Loading shapefiles")
  
  #Load shapefiles
  globe <- readOGR("data/shp/gadm36_0.shp")
  wanted_countries <- read.csv("data/ISO_use.csv", stringsAsFactors = FALSE)[, 1]
  total_countries <- subset(globe, GID_0 %in% wanted_countries)
  
  message("Loading raster data")
  
  #Load global raster
  global_raster_data <- as.data.frame(fread("output/global_raster/global_raster_wide_LC.csv"))
  LC_data <- as.matrix(global_raster_data[, grepl("LC_", colnames(global_raster_data))])
  LC_data[is.na(LC_data)] <- 0
  global_raster_data[, grepl("LC_", colnames(global_raster_data))] <- LC_data
  global_raster_data$year <- 2020
  for(i in 1:ncol(global_raster_data)){
    global_raster_data[which(is.infinite((global_raster_data)[, i]) | is.nan((global_raster_data)[, i]) | is.na((global_raster_data)[, i])), i] <- 0
  }
  
  message("Loading data to train model on")
  
  #Load stephensi environment data and make sure it all looks good
  data <- read.csv("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")
  colnames(data)[3] <- "id"
  data <- data[-which(data$ISO == "SAU"), ]
  data_fit_to <- read.csv("output/vector_competence/malariasimulation_stephensi_density_uncertainty.csv", stringsAsFactors = FALSE)
  
  data$id[-which(data$id %in% data_fit_to$id)]
  
  
  data <- left_join(data, data_fit_to, by = "id")
  
  data_stephensi_total <- data
  data_stephensi_total[which(is.na(data_stephensi_total$id)), ]$id <- paste0("nostephensi", 1:nrow(data_stephensi_total[which(is.na(data_stephensi_total$id)), ]))
  
  for(i in 1:ncol(data_stephensi_total)){
    data_stephensi_total[which(is.na(data_stephensi_total[, i])), i] <- 0
  }
  
  data_stephensi_total$stephensi_present <- 0
  data_stephensi_total$stephensi_present[which(data_stephensi_total$vectordensity_number_mid != 0)] <- 1
  
  
  #Work out model formula
  use_these_covariates_additional <- strsplit(covariate_formula, ";")[[1]]
  
  data_stephensi_total$year <- sapply(data_stephensi_total$year, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  data_stephensi_total$lat <- sapply(data_stephensi_total$lat, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  data_stephensi_total$lon <- sapply(data_stephensi_total$lon, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  
  data_stephensi_total$stephensi_present <- 0
  data_stephensi_total$stephensi_present[which(data_stephensi_total$EIR_number_high != 0)] <- 1
  
  run_cov_run <- if(any(use_these_covariates_additional %in% c("lat", "lon"))) use_these_covariates_additional[-which(use_these_covariates_additional %in% c("lat", "lon"))] else use_these_covariates_additional
  data_stephensi_total <- aggregate(data_stephensi_total[, run_cov_run],
                                    by = list(id = data_stephensi_total$id,
                                              ISO = data_stephensi_total$ISO,
                                              lat = data_stephensi_total$lat,
                                              lon = data_stephensi_total$lon),
                                    FUN = mean, na.rm = T)
  
  data_stephensi_total <- left_join(data_stephensi_total, data_fit_to, by = "id")
  data_stephensi_total <- data_stephensi_total[-which(data_stephensi_total$vectordensity_number_mid > 20), ]
  
  #Load in pre-determined which rows take
  which_rows_take <- read.csv(paste0("data/modelling_data/row_sample/which_row_take_sample_spatial_sampling_degree_", degree_take_samples, "_pseudoratio_", pseudoratio, "_refined.csv"))
  
  #Run through each spatial block sample of data
  run_row <- if(row_take_samples == "all") "all" else as.numeric(unlist(strsplit(as.character(row_take_samples), ";")))
  
  ek <- sapply(run_row, function(what_row_take){
    
    message(paste0("Sample row: ", what_row_take))
    
    #Work out which row to take
    if(what_row_take == "all"){
      data_stephensi <- data_stephensi_total
    } else {
      these_rows <- unlist(strsplit(which_rows_take$what_row_take[as.numeric(what_row_take)], ";"))
      data_stephensi <- data_stephensi_total[data_stephensi_total$id %in% these_rows, ]
    }

    #Set up training and fitting
    fitControl <- trainControl(method = "repeatedcv",
                               repeats = 10,
                               number = 3,
                               returnResamp = "final",
                               savePredictions = TRUE,
                               allowParallel = FALSE)
    
    #Scale
    for(i in use_these_covariates_additional){
      data_stephensi[, i] <- as.numeric(as.character(data_stephensi[, i]))
    }
    
    scale_this <- scale(data_stephensi[, use_these_covariates_additional], center = FALSE)
    data_stephensi[, use_these_covariates_additional] <- scale_this
    
    for(i in 1:ncol(data_stephensi)){
      data_stephensi[which(is.na(data_stephensi[, i])), i] <- 0
    }
    
    #Permutation will be ignored if the model is not ranger
    if(any(which(use_these_covariates_additional %in% fit_to))) use_these_covariates_additional <- use_these_covariates_additional[-which(use_these_covariates_additional %in% fit_to)]
    
    if(model_type == "ranger"){
      model_trained <- invisible(caret::train(form = as.formula(paste0(fit_to, " ~ ", paste(use_these_covariates_additional, collapse = " + "))),
                                              data = data_stephensi,
                                              trControl = fitControl,
                                              method = model_type,
                                              metric = criteria_select,
                                              importance  = "permutation"))
    } else {
      model_trained <- invisible(caret::train(form = as.formula(paste0(fit_to, " ~ ", paste(use_these_covariates_additional, collapse = " + "))),
                                              data = data_stephensi,
                                              trControl = fitControl,
                                              metric = criteria_select,
                                              method = model_type))
    }
    
    population_data <- raster(list.files("data/modelling_data/population", full.names = T, pattern = "2020"))
    
    message("Model fit and now predicting to countries")
    
    #Create folders
    save_here <- paste0("output/model_predictions/", saveid, "/", fit_to, "/", model_type, "/")
    if(!dir.exists(paste0(save_here, "trained_model_metric"))) dir.create(paste0(save_here, "trained_model_metric"), recursive = TRUE)
    
    #Save model
    model_trained[[4]]$criteria_select <- criteria_select
    fwrite(model_trained[[4]], paste0(paste0(save_here, "trained_model_metric/"), "best_tune_row_", what_row_take, ".csv"))
    varimp <- varImp(model_trained)[[1]]
    fwrite(data.frame(variable = row.names(varimp), importance = varimp[, 1]), paste0(paste0(save_here, "trained_model_metric/"), "variable_importance_row_", what_row_take, ".csv"))
    fwrite(data.frame(row = what_row_take, model_type = model_type, observed = data_stephensi[, fit_to], prediction = predict(model_trained)), paste0(paste0(save_here, "trained_model_metric/"), "model_prediction_vs_data_", what_row_take, ".csv"))
    
    #Run for each ISO 
    if(predict == T){
      prediction_versions <- sapply(unlist(strsplit(ISO, ";")), function(y){
        
        message(paste0("Sample row: ", y, " ", what_row_take))
        
        country_population_data <- rasterize(total_countries[which(total_countries$GID_0 == y), ],
                                             crop(population_data, extent(total_countries[which(total_countries$GID_0 == y), ])),
                                             mask = T)
        
        template <- raster(extent(country_population_data), crs = crs(country_population_data), resolution = 0.1)
        pts <- as(country_population_data, "SpatialPoints")
        vals <- raster::extract(country_population_data, pts)
        pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
        template_resampled <- rasterize(pts2, template, field = "vals", fun = mean, na.rm = T)
        
        country_dim <- global_raster_data[which(global_raster_data$ISO == y), ]
        colnames(country_dim)[grepl("LC_", colnames(country_dim))] <- gsub("_", "", colnames(country_dim)[grepl("LC_", colnames(country_dim))])
        
        country_dim_scaled <- as.data.frame(scale(country_dim[, use_these_covariates_additional], center = FALSE, scale = attr(scale_this, "scale")))
        for(i in 1:ncol(country_dim_scaled)){
          country_dim_scaled[which(is.infinite(country_dim_scaled[, i]) | is.nan(country_dim_scaled[, i]) | is.na(country_dim_scaled[, i])), i] <- 0
        }
        
        #Run models
        model_predict_ISO <- predict(model_trained, newdata = country_dim_scaled)#as.data.frame(country_dim_scaled[, use_these_covariates_additional]))
        
        country <- data.frame(cell = 1:(country_dim$ncol[1] * country_dim$nrow[1]),
                              value = rep(NA, country_dim$ncol[1] * country_dim$nrow[1]))
        country[which(country$cell %in% country_dim$cell), ]$value <- model_predict_ISO
        
        rast <- raster(matrix(country$value, ncol = country_dim$ncol[1], nrow = country_dim$nrow[1], byrow = T))
        rast_country <- rasterize(total_countries[which(total_countries$GID_0 == y), ], 
                                  crop(template_resampled, extent(total_countries[which(total_countries$GID_0 == y), ])), 
                                  mask = T)
        extent(rast) <- extent(total_countries[which(total_countries$GID_0 == y), ])
        crs(rast) <- crs(total_countries[which(total_countries$GID_0 == y), ])
        z <- resample(rast, rast_country, method = "ngb")
        if(!dir.exists(paste0(save_here, y, "/"))) dir.create(paste0(save_here, y, "/"), recursive = TRUE)
        writeRaster(z, paste0(save_here, y, "/country_raster_row_", what_row_take, ".tif"), overwrite = T)
        
      })
    }
    
  })  
  
}


