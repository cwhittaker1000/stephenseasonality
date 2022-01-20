


vector_competence_cluster <- function(species, mv, itn_cov, irs_cov, savemod = "", savefolder = ""){
  
  options(scipen = 999)
  
  #Load in LHC data
  vector_LHC_data <- as.data.frame(fread("data/mosquito/LHC/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE))
  species_LHC_data <- vector_LHC_data[which(vector_LHC_data$species == species), ]
  
  #Run all LHC forumlations
  sapply(as.numeric(strsplit(itn_cov, ";")[[1]]), function(a){
    sapply(as.numeric(strsplit(irs_cov, ";")[[1]]), function(b){
      all_LHC_run_done <- as.data.frame(rbindlist(sapply(unique(species_LHC_data$LHC_sample), function(x){
        
        message(paste0(species, " - ", x))
        
        species_LHC_data_row <- species_LHC_data[which(species_LHC_data$LHC_sample == x), ]
        
        # creates the odin model
        density_vec <- rep(mv, 365 * 30)
        
        # out <- create_r_model(odin_model_path = "odin/odin_model.R",
        #                      #Set up age and heterogeneity brackets
        #                      het_brackets = 3,
        #                      age = c(0, 1, 2, 5, 10, 20, 80),
        #                      #No country or seasonality data
        #                      admin_unit = NULL,
        #                      country = NULL,
        #                      #Background EIR which determines end prevalence
        #                      init_EIR = 1,
        #                      #Treatment/ITN/IRS turn on and coverage
        #                      num_int = 4,
        #                      ITN_IRS_on = 2, #When it turns on
        #                      init_ft = 0,
        #                      itn_cov = a,
        #                      irs_cov = b,
        #                      #Mosquito life history parameters
        #                      Q0 = if(species_LHC_data_row$Q0 == 0) min(vector_LHC_data$Q0[which(vector_LHC_data$Q0 != 0)]) else species_LHC_data_row$Q0,
        #                      chi = species_LHC_data_row$chi,
        #                      bites_Bed = species_LHC_data_row$bites_Bed, 
        #                      bites_Indoors = species_LHC_data_row$bites_Indoors)
        
        wh <- create_r_model_epidemic(odin_model_path = "odin/odin_model_seasonality.R",
                                       het_brackets = 3,
                                       age = c(0, 1, 2, 5, 10, 20, 80),
                                       num_int = 4,
                                       init_ft = 0,
                                       init_EIR = 1,
                                       country = NULL,#country,
                                       admin_unit = NULL,#admin_unit,
                                       scalar = 1,
                                       # mu0 = mu0,
                                       Q0 = if(species_LHC_data_row$Q0 == 0) min(vector_LHC_data$Q0[which(vector_LHC_data$Q0 != 0)]) else species_LHC_data_row$Q0,
                                       chi = species_LHC_data_row$chi,
                                       bites_Bed = species_LHC_data_row$bites_Bed,
                                       bites_Indoors = species_LHC_data_row$bites_Indoors,
                                       custom_seasonality = NA,
                                       time_length = length(density_vec),
                                       density_vec = density_vec,
                                       delayMos = 10)
        
        did_it_work <- try({
          
          # generates model functions with initial state data
          mod <- wh$generator(user = wh$state, use_dde = TRUE)
          
          # Runs the model
          mod_run <- mod$run(t = 1:(365 * 30))
          out <- mod$transform_variables(mod_run)
          
          
          EIR_value <- sapply(1:length(density_vec), function(x) Reduce("+", out$EIR[x, , , ]))
          # EIR <- EIR_value#/(length(wh$state$age) * length(wh$state$het_wt) * length(wh$state$num_int))
          
          data.frame(prev = mean(out$prev[(365*27):(365*30)]),
                     prev_2_10 = mean(out$prev_2to10[(365*27):(365*30)]),
                     incidence = mean(out$Incidence[(365*27):(365*30)]),
                     annual_incidence = sum(out$Incidence[(365*27):(365*30)])/3,
                     vector_density = mean(out$mv[(365*27):(365*30)]),
                     EIR = mean(EIR_value[(365*27):(365*30)]),
                     sporozoite_rates = mean(out$Iv[(365*27):(365*30)] / out$mv[(365*27):(365*30)]))
          
        })
        
        if(class(did_it_work) == "try-error"){
          data.frame(prev = NA,
                     prev_2_10 = NA,
                     incidence = NA,
                     annual_incidence = NA,
                     vector_density = NA,
                     EIR = NA,
                     sporozoite_rates = NA)
        } else {
          did_it_work
        }
          
        
      }, simplify = FALSE)))
      
      species_processed <- data.frame(species = species,
                                      mv_base = mv,
                                      itn_cov = a,
                                      irs_cov = b,
                                      LHC_rows = length(unique(species_LHC_data$LHC_sample)),
                                      t(as.data.frame(apply(all_LHC_run_done, 2, median, na.rm = T))))
      
      save_this_words <- if(savemod == ""){
        paste0("output/vector_competence/", savefolder, "/vector_competence_", species, "_", mv, "_", a, "_", b ,".csv") 
      } else{
        paste0("output/vector_competence/", savefolder, "/vector_competence_", species, "_", mv, "_", a, "_", b , "_", savemod, ".csv")
      }
      
      create_dir <- paste(strsplit(save_this_words, "/")[[1]][1:(length(strsplit(save_this_words, "/")[[1]])-1)], collapse = "/")
      if(!dir.exists(create_dir)) dir.create(create_dir, recursive = T)
      
      write.csv(species_processed, 
                save_this_words, 
                row.names = FALSE)
    })
  })
  
}






