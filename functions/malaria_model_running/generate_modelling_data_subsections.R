
generate_modelling_data_subsections <- function(pseudoabsence_ratio, grid_degree){
  
  #Load stephensi environment data
  set.seed(1)
  
  #Load stephensi environment data and make sure it all looks good
  data <- read.csv("data/processed_data/modelling_df/vector_density_environmental_covariate_with_without_stephensi.csv")
  colnames(data)[1] <- "id"
  data <- data[, 1:which(colnames(data) == "stephensi_present")]
  
  data_fit_to <- read.csv("output/vector_competence/malariasimulation_stephensi_density_uncertainty.csv", stringsAsFactors = FALSE)
  
  data <- left_join(data, data_fit_to, by = "id")
  
  data_stephensi_total <- data
  if(any(is.na(data_stephensi_total$id))) data_stephensi_total[which(is.na(data_stephensi_total$id)), ]$id <- paste0("nostephensi", 1:nrow(data_stephensi_total[which(is.na(data_stephensi_total$id)), ]))
  
  for(i in 1:ncol(data_stephensi_total)){
    data_stephensi_total[which(is.na(data_stephensi_total[, i])), i] <- 0
  }
  
  data_stephensi_total$year <- sapply(data_stephensi_total$year, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  # data_stephensi_total$lat <- sapply(data_stephensi_total$lat, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  # data_stephensi_total$lon <- sapply(data_stephensi_total$lon, function(x) mean(as.numeric(strsplit(x, ";|:")[[1]]), na.rm = T))
  
  data_stephensi_total$stephensi_present <- 0
  data_stephensi_total$stephensi_present[which(data_stephensi_total$EIR_number_high != 0)] <- 1
  
  data_stephensi_total <- aggregate(data_stephensi_total[, which(colnames(data_stephensi_total) == "year"):which(colnames(data_stephensi_total) == "stephensi_present")],
                                    by = list(id = data_stephensi_total$id,
                                              lat = data_stephensi_total$lat,
                                              lon = data_stephensi_total$lon),
                                    FUN = mean)
  
  #Now work out which presence/absence points are within 0.5 degrees of other points, may be
  #oversampling from areas that are very similar
  data_stephensi_total$row_number <- 1:nrow(data_stephensi_total)
  
  data_stephensi_total$lat <- as.character(data_stephensi_total$lat)
  data_stephensi_total$lon <- as.character(data_stephensi_total$lon)
  
  data_stephensi_total <- as.data.frame(rbindlist(sapply(1:nrow(data_stephensi_total), function(x){
    # print(x)
    this_row <- data_stephensi_total[x, ]
    all_but_this_row <- data_stephensi_total[(1:nrow(data_stephensi_total))[-x], ]
    coords_all <- SpatialPoints(matrix(c(as.numeric(unlist(strsplit(all_but_this_row$lon, ";"))), 
                                         as.numeric(unlist(strsplit(all_but_this_row$lat, ";")))), ncol = 2))
    spdf_not_row <- SpatialPointsDataFrame(coords_all, data.frame(row = rep(all_but_this_row$row_number, sapply(all_but_this_row$lon, function(x) length(strsplit(x, ";")[[1]])))))
    gbuffer_ind <- rgeos::gBuffer(spdf_not_row, grid_degree, byid = T, id = 1:nrow(spdf_not_row))
    
    this_point_coord <- SpatialPoints(matrix(c(as.numeric(unlist(strsplit(this_row$lon, ";"))), 
                                               as.numeric(unlist(strsplit(this_row$lat, ";")))), ncol = 2))
    this_point_buff <- rgeos::gBuffer(SpatialPointsDataFrame(this_point_coord, 
                                                             data = data.frame(row = rep(x, length(this_point_coord)))), grid_degree, byid = F)
    which_intersect <- as.data.frame(intersect(gbuffer_ind, this_point_buff))
    
    this_row$intersects_with <- paste(sort(c(unique(x), unique(which_intersect$row))), collapse = ";")
    this_row
    
  }, simplify = FALSE)))
  
  
  #Account for geographic bias
  expand_df <- as.data.frame(data.table::rbindlist(sapply(1:nrow(data_stephensi_total), function(x){
    
    this_row <- data_stephensi_total[x, ]
    if(grepl(";", this_row$lat)){
      this_row <- this_row[rep(seq_len(nrow(this_row)), each = length(unlist(strsplit(this_row$lat, ";")))), ]
    }
    
    this_row$lat <- as.numeric(strsplit(this_row$lat, ";")[[1]])
    this_row$lon <- as.numeric(strsplit(this_row$lon, ";")[[1]])
    
    this_row$lat_round <- plyr::round_any(this_row$lat, grid_degree)
    this_row$lon_round <- plyr::round_any(this_row$lon, grid_degree)
    this_row$lon_lat_round <- paste0(this_row$lon_round, "_", this_row$lat_round)
    this_row
    
  }, simplify = FALSE)))
  
  #Now select by their unique grid reference
  unique_grid <- unique(expand_df$lon_lat_round)
  
  tz <- do.call(rbind, sapply(unique_grid, function(x){
    if(any(subset(expand_df, lon_lat_round == x)$stephensi_present == 1)) 1 else 0
  }, simplify = FALSE))
  
  all_100_row <- as.data.frame(rbindlist(sapply(1:100, function(a){
    # print(a)
    done <- do.call(rbind, sapply(unique_grid, function(y){
      these_values <- unique(expand_df[which(expand_df$lon_lat_round == y), ]$id)
      this_row <- if(length(these_values) == 1) these_values else sample(these_values, 1)
      data.frame(run = a,
                 unique_grid = y,
                 row_take = this_row,
                 stephensi_present = if(any(expand_df[which(expand_df$id == this_row), ]$stephensi_present == 1)) 1 else 0)
    }, simplify = FALSE))
    #Make presence and absence the same
    rbind(done[which(done$stephensi_present == 1), ],
          done[sample(which(done$stephensi_present == 0), round(nrow(done[which(done$stephensi_present == 1), ]) * pseudoabsence_ratio), replace = T), ])
    
  }, simplify = FALSE)))
  
  
  all_100_row$lon <- as.numeric(sapply(strsplit(all_100_row$unique_grid, "_"), function(g) g[1]))
  all_100_row$lat <- as.numeric(sapply(strsplit(all_100_row$unique_grid, "_"), function(g) g[2]))
  all_100_row$size <- NA
  for(i in unique(all_100_row$row_take)){
    all_100_row[which(all_100_row$row_take == i), ]$size <- nrow(all_100_row[which(all_100_row$row_take == i), ])
  }
  
  #Data frame
  run_100_samples_geography <- data.frame(run = 1:100, 
                                          grid_degree = grid_degree,
                                          pseudoabsence_ratio = pseudoabsence_ratio,
                                          what_row_take = sapply(1:max(all_100_row$run), function(x) paste0(all_100_row[which(all_100_row$run == x), ]$row_take, collapse = ";")),
                                          stephensi_present = sapply(1:max(all_100_row$run), function(x) sum(all_100_row[which(all_100_row$run == x), ]$stephensi_present)),
                                          stephensi_absent = sapply(1:max(all_100_row$run), function(x) nrow(all_100_row[which(all_100_row$run == x), ]) - sum(all_100_row[which(all_100_row$run == x), ]$stephensi_present)))
  
  write.csv(run_100_samples_geography, paste0("data/modelling_data/row_sample/which_row_take_sample_spatial_sampling_degree_", grid_degree, "_pseudoabsence_ratio_", pseudoabsence_ratio, "_refined.csv"), row.names = FALSE)
  
  
}