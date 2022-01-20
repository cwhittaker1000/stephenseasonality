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
# vectordensity_fit = run_these[1, "vectordensity_fit"]
# species = run_these[1, "species"]
# species_abundance = run_these[1, "species_abundance"]
# itn_cov = run_these[1, "itn_cov"]
# irs_cov = run_these[1, "irs_cov"]
# runtime = run_these[1, "runtime"]
# total_runs = run_these[1, "total_runs"]
# n_tries = run_these[1, "n_tries"]
# save_name_fill = ""

fit_malariasimulation_EIR_vectordensity <- function(id,
                                                 row,
                                                 human_population = 1000, 
                                                 vectordensity_fit, 
                                                 species, 
                                                 species_abundance, 
                                                 itn_cov, 
                                                 irs_cov, 
                                                 runtime = 365 * 6, 
                                                 total_runs,
                                                 n_tries,
                                                 save_name_fill = ""){
  
  
  sapply(unlist(strsplit(as.character(row), ";")), function(row_use){
    
    message("Row: ", row_use)
    
    EIR_use <- if(file.exists(save_name_fill)){
      old_df <- read.csv(save_name_fill)
      old_df$EIR_use
    } else {
      10
    }
    
    #Initial guestimate fit
    guestimate_fit <- malariasimulation_vectordensity_relationship_shoddy_fit(row = row_use,
                                                                               EIR = EIR_use,
                                                                               vectordensity_fit = vectordensity_fit,
                                                                               species = species, 
                                                                               species_abundance = species_abundance, 
                                                                               itn_cov = itn_cov, 
                                                                               irs_cov = irs_cov, 
                                                                               runtime = runtime,
                                                                               total_runs = total_runs,
                                                                               n_tries = n_tries)
    
    
    difference <- last(na.omit(as.numeric(strsplit(guestimate_fit$differences, ";")[[1]])))
    EIR_try <- last(na.omit(as.numeric(strsplit(guestimate_fit$EIR_tried, ";")[[1]])))
    
    output_df <- data.frame(id = id,
                            total_runs = total_runs,
                            row = row_use,
                            runtime = runtime,
                            human_population = human_population,
                            itn_cov = itn_cov,
                            irs_cov = irs_cov,
                            species = species,
                            species_abundance = species_abundance,
                            prevalence = guestimate_fit$mean_prevalence,
                            malaria_vectordensity_fit_to = vectordensity_fit,
                            total_M_guessed = guestimate_fit$total_M,
                            EIR_guessed = last(unlist(strsplit(guestimate_fit$EIR_tried, ";"))),
                            malaria_vectordensity_estimated = vectordensity_fit - as.numeric(difference),
                            malaria_vectordensity_difference = difference)
    
    savename <- paste(row_use, human_population, itn_cov, irs_cov, id, runtime, total_runs, sep = "_")
    if(!dir.exists(paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/"))) dir.create(paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/"), recursive = T)
    save_name_fill <- paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/", savename, ".csv")
    write.csv(output_df, save_name_fill, row.names = FALSE)
    
  })
  
}
