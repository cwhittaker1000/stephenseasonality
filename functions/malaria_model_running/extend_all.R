extend_all <- function(rasters){
  extent(Reduce(extend,rasters))
}