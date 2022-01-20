
process_end_model_predictions_to_one <- function(model_types, fit_to, id){
  
  #Load in all model predictions
  all_model_pred_loaded <- sapply(unlist(strsplit(model_types, ";")), function(mod){
    raster(list.files(paste0("U:/Arran/stephensi_density/output/model_predictions/", id, "/", fit_to, "/", mod),
                      pattern = "global_mean", 
                      full.names = T))
  }, simplify = T)
  
  #Load in all model metrics
  all_model_metric_loaded <- as.data.frame(rbindlist(sapply(unlist(strsplit(model_types, ";")), function(mod){
    best <- read.csv(list.files(paste0("U:/Arran/stephensi_density/output/model_predictions/", id, "/", fit_to, "/", mod), "best", full.names = T))
    df <- rbindlist(sapply(list.files(paste0("U:/Arran/stephensi_density/output/model_predictions/", id, "/", fit_to, "/", mod, "/trained_model_metric/"),
                                      pattern = "best_tune", 
                                      full.names = T), function(y){
                                        d <- read.csv(y)
                                        d$row <- last(gsub(".csv", "", unlist(strsplit(y, "_"))))
                                        d[which.min(d$RMSE), ]
                                      }, simplify = FALSE))
    df <- df[which(df$row %in% if(nrow(best) == 0) "all" else best$row), ]
    df$mod <- mod
    df$weight <- if(nrow(best) == 0) 1 else (max(df$RMSE) - df$RMSE)/sum((max(df$RMSE) - df$RMSE))#df$RMSE/sum(df$Rsquared)
    df$overall_r2 <-  if(nrow(best) == 0) 1 else sum(df$weight * df$Rsquared)
    df
  }, simplify = F), fill = T))
  
  weight_value <- sapply(unique(all_model_metric_loaded$mod), function(x){
    sum((max(all_model_metric_loaded$RMSE) - all_model_metric_loaded[which(all_model_metric_loaded$mod == x), ]$RMSE)/sum((max(all_model_metric_loaded$RMSE) - all_model_metric_loaded$RMSE)))
  })
  
  all_model_pred_loaded_weighted <- sapply(names(all_model_pred_loaded), function(x){
    all_model_pred_loaded[[x]] * weight_value[x]
  }, simplify = FALSE)
  
  pred_combo <- Reduce("+", all_model_pred_loaded_weighted)
  
  #Load variable importance
  all_variable_importance_loaded <- as.data.frame(rbindlist(sapply(unlist(strsplit(model_types, ";")), function(mod){
    print(mod)
    best <- read.csv(list.files(paste0("U:/Arran/stephensi_density/output/model_predictions/", id, "/", fit_to, "/", mod), "best", full.names = T))
    df <- as.data.frame(rbindlist(sapply(list.files(paste0("U:/Arran/stephensi_density/output/model_predictions/", id, "/", fit_to, "/", mod, "/trained_model_metric/"),
                                                    pattern = "importance", 
                                                    full.names = T), function(y){
                                                      d <- read.csv(y)
                                                      d$row <- last(gsub(".csv", "", unlist(strsplit(y, "_"))))
                                                      weight <- all_model_metric_loaded[which(all_model_metric_loaded$mod == mod & all_model_metric_loaded$row == unique(d$row)), ]$weight
                                                      d$weight <- ifelse(length(weight) == 0, 0, weight)
                                                      d
                                                    }, simplify = FALSE)))
    df$mod <- mod
    df <- df[which(df$row %in% best$row), ]
    df$importance_weighted <- df$importance * df$weight
    df
  }, simplify = F), fill = T))
  
  all_variable_importance_loaded_combined <- aggregate(x = list(importance_weighted = all_variable_importance_loaded$importance_weighted),
                                                       by = list(variable = all_variable_importance_loaded$variable,
                                                                 mod = all_variable_importance_loaded$mod),
                                                       FUN = sum)
  #Save
  writeRaster(pred_combo, paste0("output/model_predictions/", id, "/", fit_to, "/combined_prediction.tif"), overwrite = T)
  
  #Output
  list(prediction = pred_combo,
       variable_importance = all_variable_importance_loaded_combined,
       R2 = sum(sapply(unique(all_model_metric_loaded$mod), function(x) all_model_metric_loaded[which(all_model_metric_loaded$mod == x), ]$overall_r2[1] * weight_value[names(weight_value) == x])),
       individual_predictions = all_model_pred_loaded,
       individual_metrcs = all_model_metric_loaded)
  
}
