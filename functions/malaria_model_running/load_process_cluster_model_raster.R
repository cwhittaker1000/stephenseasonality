# fit_to = "vector_density_malaria_mid"
# model_type = "ranger"
# saveid = "GA1FZE"

load_process_cluster_model_raster <- function(fit_to,
                                              model_type,
                                              saveid,
                                              ISO = "",
                                              type = "median"){
  
  #Find all model runs
  all_run_loaded <- list.files(if(ISO == ""){
    paste0("output/model_predictions/", fit_to, "/", model_type, "/", saveid, "/")
  } else {
    paste0("output/model_predictions/", fit_to, "/", model_type, "/", saveid, "/", ISO, "/")
  },
  pattern = "country_raster_row_",
  full.names = T, recursive = T)
  
  all_run_loaded <- all_run_loaded[grepl("tif", all_run_loaded)]
  if(any(grepl("full_run", all_run_loaded))) all_run_loaded <- all_run_loaded[!grepl("full_run", all_run_loaded)]
    
  #Load in
  loaded_all <- sapply(all_run_loaded, function(x){raster(x)}, simplify = FALSE)
  
  #Process
  names(loaded_all) <- NULL
  loaded_all$fun <- if(type == "median") median else if(type == "sd") sd
  loaded_total <- do.call(mosaic, loaded_all)
  
  #Save
  save_name <- strsplit(all_run_loaded[1], "/")[[1]]
  writeRaster(loaded_total, 
              paste0(paste(save_name[1:(length(save_name)-2)], collapse = "/"), "/",
                     if(ISO == "") "all_countries_" else paste0(ISO, "_"),
                     fit_to, "_", model_type, "_", saveid, "_", type, "_full_run.tif"))
  
}

