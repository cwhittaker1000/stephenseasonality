
cluster_load_combine_predictions <- function(fit_to, type, saveid, ISO, number_include = 10, criteria_select = "RMSE"){
  
  sapply(unlist(strsplit(ISO, ";")), function(a){
    
    message(a)
    
    all_these <- list.files(paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/", a),
                            pattern = "country_raster_row_",
                            full.names = T)
    
    all_these_metric <- list.files(paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/trained_model_metric"),
                                   pattern = "best_tune",
                                   full.names = T)
    
    all_these_metric_best <- as.data.frame(rbindlist(sapply(all_these_metric, function(x){
      these <- read.csv(x)
      data.frame(row = gsub(".csv", "", last(unlist(strsplit(x, "_")))), r2 = these[which.min(these$RMSE), ]$Rsquared,
                 RMSE = these[which.min(these$RMSE), ]$RMSE  )
    }, simplify = FALSE)))
    
    all_these_metric_best$RMSE_weight <- (max(all_these_metric_best$RMSE) - all_these_metric_best$RMSE)/sum((max(all_these_metric_best$RMSE) - all_these_metric_best$RMSE))
    all_these_metric_best$R2_weight <- all_these_metric_best$r2/sum(all_these_metric_best$r2)
    
    #Take the 10 best
    if(criteria_select == "RMSE"){
      best_10 <- all_these_metric_best[which(all_these_metric_best$RMSE_weight %in% rev(sort(all_these_metric_best$RMSE_weight))[1:number_include]), ]
    } else if(criteria_select == "r2"){
      best_10 <- all_these_metric_best[which(all_these_metric_best$r2 %in% (sort(all_these_metric_best$r2))[1:number_include]), ]
    }
    
    if(all(all_these_metric_best$row == "all")){
      best_10 <- all_these_metric_best
      best_10$weight_use <- 1
    } else {
      best_10$weight_use <- if(criteria_select == "RMSE") best_10$RMSE_weight/sum(best_10$RMSE_weight) else if(criteria_select == "r2") best_10$R2_weight/sum(best_10$R2_weight)
    }
    
    best_model_go <- paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/best_", number_include, "_models_metrics.csv")
    if(!file.exists(best_model_go)) write.csv(best_10, best_model_go, row.names = FALSE)
    
    all_these_subset <- all_these[grepl(paste(paste0("row_", best_10$row, ".tif"), collapse = "|"), all_these)]
    
    #Load now
    load_these <- sapply(all_these_subset, function(x) raster(x) * best_10[which(best_10$row == gsub(".tif", "", last(strsplit(x, "_")[[1]]))), ]$weight_use, simplify = FALSE)
    names(load_these) <- NULL
    combined <- Reduce("+", load_these)
    writeRaster(combined, filename = paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/", a, "/", a, "_mean_predictions.tif"), overwrite = T)
    
  })

}


