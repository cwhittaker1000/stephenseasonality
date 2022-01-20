
cluster_processed_the_processed <- function(fit_to, type, saveid, pattern = "mean_predictions.tif"){
  
  all_these <- list.files(paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/"),
                          pattern = pattern,
                          full.names = T,
                          recursive = T)
  
  if(any(grepl("overall_global", all_these))) all_these <- all_these[-which(grepl("overall_global", all_these))]
  
  load_these <- sapply(all_these, function(x) raster(x), simplify = FALSE)
  names(load_these) <- NULL
  load_these$fun <- mean
  load_these$na.rm <- TRUE
  load_these$tolerance <- 1
  
  combined <- if(length(all_these) == 1) load_these[[1]] else do.call(merge, load_these)
  writeRaster(combined, filename = paste0("output/model_predictions/", saveid, "/", fit_to, "/", type, "/overall_global_mean_predictions.tif"), overwrite = T)
  
}

