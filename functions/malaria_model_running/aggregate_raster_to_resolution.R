
aggregate_raster_to_resolution <- function(raster, resolution, function_use){
  
  template <- raster(extent(raster), crs = crs(raster), resolution = resolution)
  pts <- as(raster, "SpatialPoints")
  vals <- raster::extract(raster, pts)
  pts2 <- SpatialPointsDataFrame(pts, data.frame(vals))
  raster_resampled <- rasterize(pts2, template, field = "vals", fun = function_use, na.rm = T)
  raster_resampled
  
}


