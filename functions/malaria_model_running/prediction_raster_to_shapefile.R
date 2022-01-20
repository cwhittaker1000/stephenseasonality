
# raster_location <- "output/model_predictions/EIR_stephensi/hCG34o_combined_prediction.tif"
# ISO <- "ETH"
# savemod = "test"

prediction_raster_to_shapefile <- function(raster_location, admin = "1", ISO, savemod, falciparum_mask = F){
  
  globe <- readOGR(paste0("data/shp/gadm36_", admin, ".shp"))
  country_use <- globe[which(globe$GID_0 %in% unlist(strsplit(ISO, ";"))), ]
  pred_loaded <- raster(raster_location)
  population_data <- raster(list.files("data/modelling_data/population", full.names = T, pattern = "2020"))
  population_data_size <- rasterize(country_use, crop(population_data, extent(country_use)), mask = T)
  
  if(falciparum_mask == T){
    falciparum_here <- raster("data/mosquito/2010_Pf_Limits_Decompressed.geotiff")
    falciparum_here_crop <- rasterize(country_use, crop(falciparum_here, country_use), mask = T)
    value <- if(ISO == "DJI") c(1, 2) else 2

    falciparum_here_crop[!(falciparum_here_crop %in% value)] <- 0
    falciparum_here_crop[(falciparum_here_crop != 0)] <- 1

    template_falciparum <- raster(extent(population_data_size), crs = crs(population_data_size), resolution = res(population_data_size))
    pts_falciparum <- as(falciparum_here_crop, "SpatialPoints")
    vals_falciparum <- raster::extract(falciparum_here_crop, pts_falciparum)
    pts2_falciparum <- SpatialPointsDataFrame(pts_falciparum, data.frame(vals_falciparum))
    template_resampled_falciparum <- rasterize(pts2_falciparum, template_falciparum, field = "vals_falciparum", fun = max, na.rm = T)
    
    population_data_size <- population_data_size * template_resampled_falciparum
  }
  
  template <- raster(extent(pred_loaded), crs = crs(pred_loaded), resolution = res(pred_loaded))
  pts <- as(population_data_size, "SpatialPoints")
  vals <- raster::extract(population_data_size, pts)
  pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
  template_resampled <- rasterize(pts2, template, field = "vals", fun = sum, na.rm = T)
  
  df <- rbindlist(sapply(1:nrow(country_use), function(x){
    
    raster_combine <- rasterize(country_use[x, ], crop(pred_loaded, extent(country_use[x, ])), mask = T)
    pop_combine <- rasterize(country_use[x, ], crop(template_resampled, extent(country_use[x, ])), mask = T)
    
    pop_weighted <- raster_combine * pop_combine/sum(pop_combine[], na.rm = T)
    
    data.frame(country_use[x, grepl("NAME|GID", names(country_use))],
               population = sum(pop_combine[], na.rm = T),
               pred_pop_weighted = sum(pop_weighted[], na.rm = T))
    
  }, simplify = FALSE))
  
  save_loc <- paste0(gsub("combined_prediction.tif", "", raster_location), "admin_data/", ISO, "/")
  if(!dir.exists(save_loc)) dir.create(save_loc, recursive = T)
  
  fwrite(df, paste0(save_loc, "country_data_admin_", admin, "_", savemod, ".csv"))
  
}