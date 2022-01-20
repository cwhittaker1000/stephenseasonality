# species <- "stephensi;dirus"
# species_abundance <- "0.5;0.5"
# EIR = 50
# itn_cov = 0.1
# irs_cov = 0
# runtime = 365 * 6
# run = 1
# total_runs = 3
# malaria_prevalence_fit = 0.02
# fit_or_results = "fit"

malariasimulation_EIR_prevalence_relationship <- function(EIR, malaria_prevalence_fit, species, 
                                                          species_abundance, itn_cov, irs_cov, runtime = 365 * 6, 
                                                          total_runs, fit_or_results = "fit"){
  
  #Load in vector bionomic data
  vector_bionomic_data <- read.csv("data/mosquito/processed_species_bionomics.csv", stringsAsFactors = F)
  vector_bionomic_data <- vector_bionomic_data[-which(is.na(vector_bionomic_data$species)), ]
  
  #Set up baseline dat 6
  year <- 365
  month <- 30
  sim_length <- 1 * year
  human_population <- 1000
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      model_seasonality = FALSE, 
      prevalence_rendering_min_ages = 2 * 365, #2-10 year old prevalence for map
      prevalence_rendering_max_ages = 10 * 365, #2-10 year old prevalence for map
      incidence_rendering_min_ages = 0,
      incidence_rendering_max_ages = 100 * 365,
      severe_enabled = 1
    )
  )
  
  simparams <- set_equilibrium(simparams, EIR)
  
  #Set up species
  #Set up bionomics
  all_bionomics <- sapply(strsplit(species, ";")[[1]], function(h){
    
    if(nrow(vector_bionomic_data[which(vector_bionomic_data$species == h), ]) == 0){
      vector_bionomic_data[which(vector_bionomic_data$species == h), ]
    }
    
    list(species = h,
         blood_meal_rates = 1/3, #species invariant due to lack of data
         foraging_time = 0.69, #species invariant due to lack of data
         Q0 = vector_bionomic_data[which(vector_bionomic_data$species == h), ]$anthropophagy,
         phi_bednets = vector_bionomic_data[which(vector_bionomic_data$species == h), ]$inbed_biting,
         phi_indoors = vector_bionomic_data[which(vector_bionomic_data$species == h), ]$indoor_biting,
         mum = 0.12) #species invariant due to lack of data
  }, simplify = FALSE)
  
  simparams <- set_species(simparams, all_bionomics, as.numeric(strsplit(species_abundance, ";")[[1]]))
  
  #Implement interventions if specified
  if(itn_cov != 0){
    itn_applications <- length(seq(1, runtime, by = 365 * 3))
    simparams <- set_bednets(
      parameters = simparams,
      timesteps = seq(1, runtime, by = 365 * 3),
      coverages = rep(itn_cov, length(seq(1, runtime, by = 365 * 3))),
      retention = 5 * year,
      dn0 = matrix(rep(.533, length(all_bionomics) * itn_applications), ncol = length(all_bionomics), byrow = T),
      rn = matrix(rep(.56, length(all_bionomics) * itn_applications), ncol = length(all_bionomics), byrow = T),
      rnm = matrix(rep(.24, length(all_bionomics) * itn_applications), ncol = length(all_bionomics), byrow = T),
      gamman = matrix(rep(2.64 * 365, itn_applications), nrow = 1, byrow = T)
    )
  }
  
  
  if(irs_cov != 0){
    irs_applications <- length(seq(1, runtime, by = 365 * 1))
    simparams <- set_spraying(
      simparams,
      timesteps = seq(1, runtime, by = 365 * 1),
      coverages = rep(irs_cov, length(seq(1, runtime, by = 365 * 1))),
      ls_theta = matrix(rep(2.025, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T),
      ls_gamma = matrix(rep(-0.009, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T),
      ks_theta = matrix(rep(-2.222, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T),
      ks_gamma = matrix(rep(0.008, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T),
      ms_theta = matrix(rep(-1.232, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T),
      ms_gamma = matrix(rep(-0.009, length(all_bionomics) * irs_applications), ncol = length(all_bionomics), byrow = T)
    )
  }
  
  all_run_done <- sapply(1:total_runs, function(a){
    
    output <- run_simulation(runtime, simparams)
    
    incidence_by_year <- matrix(output[, grepl("n_inc", colnames(output))]/output[, grepl("n_0", colnames(output))], ncol = 365, byrow = T)
    prevalence_by_year <- matrix(output[, grepl("n_inc", colnames(output))]/output[, grepl("n_0", colnames(output))], ncol = 365, byrow = T)
    
    dfz_save <- data.frame(total_runs = total_runs,
                           human_population = human_population,
                           EIR = EIR,
                           itn_cov = itn_cov,
                           irs_cov = irs_cov,
                           species = species,
                           species_abundance = species_abundance,
                           malaria_prevalence_fit = malaria_prevalence_fit,
                           last_year_annual_incidence = sum(last(incidence_by_year)),
                           last_year_mean_prevalence_2_10 = mean(prevalence_by_year),
                           stringsAsFactors = FALSE)
    
    # savename <- paste(human_population, EIR, itn_cov, irs_cov, species, species_abundance, run, total_runs, sep = "_")
    # write.csv(dfz_save, paste0("output/vector_competence/malariasimulation_runs/", savename, ".csv"), row.names = FALSE)
    
    list(likelihood =  -dnbinom(round(malaria_prevalence_fit * 1000, 0), 
                                1, 
                                mu = round(dfz_save$last_year_mean_prevalence_2_10 * 1000, 0), 
                                log = TRUE),
         data = dfz_save)
    
  }, simplify = FALSE)
  
  all_likelihood <- sapply(1:length(all_run_done), function(x) all_run_done[[x]][[1]])
  # print(all_likelihood)
  mean_likelihood <- mean(all_likelihood, na.rm = T)
  mean_outcome <- do.call(rbind, sapply(1:length(all_run_done), function(x) all_run_done[[x]][[2]], simplify = FALSE))
  mean_outcome_df <- mean_outcome[1, ]
  mean_outcome_df$last_year_annual_incidence <- mean(mean_outcome_df$last_year_annual_incidence)
  mean_outcome_df$last_year_mean_prevalence_2_10 <- mean(mean_outcome_df$last_year_mean_prevalence_2_10)
  this <- (malaria_prevalence_fit - mean_outcome_df$last_year_mean_prevalence_2_10)
  
  message(paste0("EIR: ", EIR, " - ", "Prevalence difference: ", this, " - Likelihood: ", mean_likelihood))
  
  if(fit_or_results == "fit"){
    mean_likelihood
  } else {
    mean_outcome_df
  }
  
  
}










