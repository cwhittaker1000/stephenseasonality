
what_didnt_run <- function(dataset){
  
  all_these <- as.data.frame(rbindlist(sapply(unique(dataset$saveid), function(a){
    rbindlist(sapply(unique(dataset$fit_to), function(y){
      rbindlist(sapply(unique(dataset$model_type), function(x){
        # print(x)
        z <- rbindlist(sapply(unique(unlist(strsplit(dataset$ISO, ";"))), function(ISO){
          # print(ISO)
          these_ran <- list.files(paste0("output/model_predictions/", a, "/", y, "/", x, "/", ISO), "tif")
          if(length(these_ran) == 0){
            these_ran <- list.files(paste0("output/model_predictions/", a, "/", y, "/", x, "/trained_model_metric/"), pattern = "model_prediction")
            number <- as.numeric(gsub("model_prediction_vs_data_|.csv", "", these_ran))
          } else {
            number <- as.numeric(gsub("country_raster_row_|.tif", "", list.files(paste0("output/model_predictions/", a, "/", y, "/", x, "/", ISO), "tif")))
          }
          
          number_need <- (1:100)[-number]
          
          if(length(number_need) != 0){
            these_go_go <- dataset[which(dataset$saveid == a & dataset$fit_to == y & dataset$model_type == x), ][1, ]
            these_go_go <- these_go_go[rep(1, length(number_need)), ]
            these_go_go$ISO <- ISO
            these_go_go$row_take_samples <- as.character(number_need)
            these_go_go 
          }
        }, simplify = FALSE), fill = F)
        z$ISO <- paste(unique(z$ISO), collapse = ";")
        # z$row_take_samples <- paste(as.character(unique(z$row_take_samples)), collapse = ";")
        if(all(z$ISO != "")) z[!duplicated(z), ]
      }, simplify = FALSE), fill = F)
    }, simplify = FALSE))
  }, simplify = FALSE)))
  
  all_these
  
}