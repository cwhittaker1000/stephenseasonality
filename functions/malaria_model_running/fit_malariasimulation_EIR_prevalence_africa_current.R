

fit_malariasimulation_EIR_prevalence_africa_current <- function(row,
                                                                save_name_fill = ""){
  
  these_runs <- read.csv("output/africa_wide_prediction/baseline/run_specification_table.csv")
  
  sapply(unlist(strsplit(as.character(row), ";")), function(row_use){
    
    message("Row: ", row_use)
    
    id <- these_runs[row_use, "id"]
    human_population <- these_runs[row_use, "human_population"]
    malaria_prevalence_fit <- these_runs[row_use, "malaria_prevalence_fit"] 
    species <- these_runs[row_use, "species"]
    species_abundance <- these_runs[row_use, "species_abundance"]
    #Correct species
    these_species_present <- which(strsplit(species_abundance, ";")[[1]] != 0)
    species <- paste(strsplit(species, ";")[[1]][these_species_present], collapse = ";")
    species_abundance <- paste(strsplit(species_abundance, ";")[[1]][these_species_present], collapse = ";")
    itn_cov <- these_runs[row_use, "itn_cov"] 
    irs_cov <- these_runs[row_use, "irs_cov"] 
    runtime <- these_runs[row_use, "runtime"]
    total_runs <- these_runs[row_use, "total_runs"]
    n_tries <- these_runs[row_use, "n_tries"]
    fit_threshold <- these_runs[row_use, "fit_threshold"]
    
    EIR_use <- if(file.exists(save_name_fill)){
      old_df <- read.csv(save_name_fill)
      old_df$EIR_use
    } else {
      ifelse(malaria_prevalence_fit < 0.1, 10, 25)
    }
    
    #Initial guestimate fit
    guestimate_fit <- malariasimulation_EIR_prevalence_relationship_shoddy_fit(row = "mean",
                                                                               EIR = EIR_use,
                                                                               malaria_prevalence_fit = malaria_prevalence_fit,
                                                                               species = species, 
                                                                               species_abundance = species_abundance, 
                                                                               itn_cov = itn_cov, 
                                                                               irs_cov = irs_cov, 
                                                                               runtime = runtime,
                                                                               total_runs = total_runs,
                                                                               n_tries = n_tries,
                                                                               fit_threshold = fit_threshold)
    
    best_difference <- which.min(na.omit(abs(as.numeric(strsplit(guestimate_fit$differences, ";")[[1]]))))
    difference <- na.omit(as.numeric(strsplit(guestimate_fit$differences, ";")[[1]]))[best_difference]
    EIR_try <- na.omit(as.numeric(strsplit(guestimate_fit$EIR_tried, ";")[[1]]))[best_difference]
    
    output_df <- data.frame(total_runs = total_runs,
                            row = row_use,
                            runtime = runtime,
                            human_population = human_population,
                            itn_cov = itn_cov,
                            irs_cov = irs_cov,
                            species = species,
                            species_abundance = species_abundance,
                            malaria_prevalence_fit_to = malaria_prevalence_fit,
                            EIR_use = guestimate_fit$EIR,
                            rbind(guestimate_fit[, grepl("total_M_", colnames(guestimate_fit))]),
                            rbind(guestimate_fit[, grepl("EIR_", colnames(guestimate_fit))]),
                            malaria_prevalence_estimated = malaria_prevalence_fit - as.numeric(difference),
                            malaria_prevalence_difference = difference,
                            fit_threshold = fit_threshold)
    
    savename <- paste(row_use, human_population, itn_cov, irs_cov, id, runtime, total_runs, sep = "_")
    if(!dir.exists(paste0("output/africa_wide_prediction/baseline/", strsplit(id, "\\.")[[1]][1], "/"))) dir.create(paste0("output/africa_wide_prediction/baseline/", strsplit(id, "\\.")[[1]][1], "/"), recursive = T)
    save_name_fill <- paste0("output/africa_wide_prediction/baseline/", strsplit(id, "\\.")[[1]][1], "/", savename, ".csv")
    write.csv(output_df, save_name_fill, row.names = FALSE)
    
  })
  
}
