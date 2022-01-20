mean_all <- function(rasters, extent){
  re = lapply(rasters, function(r){extend(r, extent)})
  Reduce("+", re)/re
}