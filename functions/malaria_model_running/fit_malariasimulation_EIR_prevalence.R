# species <- "stephensi;dirus"
# species_abundance <- "0.5;0.5"
# EIR = 50
# itn_cov = 0.1
# irs_cov = 0
# runtime = 365 * 2
# run = 1
# total_runs = 5
# malaria_prevalence_fit = 0.02
# n_tries = 5
# human_population = 1000

# 
# #
# id = run_these[1, "id"]
# human_population = run_these[1, "human_population"]
# malaria_prevalence_fit = run_these[1, "malaria_prevalence_fit"]
# species = run_these[1, "species"]
# species_abundance = run_these[1, "species_abundance"]
# itn_cov = run_these[1, "itn_cov"]
# irs_cov = run_these[1, "irs_cov"]
# runtime = run_these[1, "runtime"]
# total_runs = run_these[1, "total_runs"]
# n_tries = run_these[1, "n_tries"]
# row = run_these[1, "row"]
# save_name_fill = run_these[1, "save_name_fill"]
# fit_threshold = run_these[1, "fit_threshold"]



fit_malariasimulation_EIR_prevalence <- function(id,
                                                 row,
                                                 human_population = 1000, 
                                                 malaria_prevalence_fit, 
                                                 species, 
                                                 species_abundance, 
                                                 itn_cov, 
                                                 irs_cov, 
                                                 runtime = 365 * 6, 
                                                 total_runs,
                                                 n_tries,
                                                 save_name_fill = "",
                                                 fit_threshold = 0.05){
  

  sapply(unlist(strsplit(as.character(row), ";")), function(row_use){
    
    message("Row: ", row_use)
    
    EIR_use <- if(file.exists(save_name_fill)){
      old_df <- read.csv(save_name_fill)
      old_df$EIR_use
    } else {
      25
    }
    
    #Initial guestimate fit
    guestimate_fit <- malariasimulation_EIR_prevalence_relationship_shoddy_fit(row = row_use,
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
    if(!dir.exists(paste0("output/vector_competence/updated_malariasimulation_runs/", row_use, "/"))) dir.create(paste0("output/vector_competence/updated_malariasimulation_runs/", row_use, "/"), recursive = T)
    save_name_fill <- paste0("output/vector_competence/updated_malariasimulation_runs/", row_use, "/", savename, ".csv")
    write.csv(output_df, save_name_fill, row.names = FALSE)
    
  })
  
}
