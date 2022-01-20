# species <- "stephensi;dirus"
# species_abundance <- "0.5;0.5"
# EIR = 50
# itn_cov = 0.1
# irs_cov = 0
# runtime = 365 * 2
# run = 1
# total_runs = 3
# vectordensity_fit = 0.02
# fit_or_results = "fit"


malariasimulation_vectordensity_relationship_shoddy_fit <- function(row, EIR, vectordensity_fit, species, 
                                                                     species_abundance, itn_cov, irs_cov, runtime = 365 * 6, 
                                                                     total_runs, n_tries){
  
  #Load in vector bionomic data
  vector_bionomic_data <- read.csv("data/mosquito/processed_species_bionomics_LHC.csv", stringsAsFactors = F)
  vector_bionomic_data <- vector_bionomic_data[which(vector_bionomic_data$run == row), ]
  vector_bionomic_data$run <- as.numeric(vector_bionomic_data$run)
  
  #Set up baseline dat 6
  year <- 365
  month <- 30
  sim_length <- 1 * year
  human_population <- 1000
  
  simparams <- get_parameters(
    list(
      individual_mosquitoes = FALSE,
      human_population = human_population,
      model_seasonality = FALSE, 
      prevalence_rendering_min_ages = 2 * 365, #2-10 year old prevalence for map
      prevalence_rendering_max_ages = 10 * 365, #2-10 year old prevalence for map
      incidence_rendering_min_ages = 0,
      incidence_rendering_max_ages = 100 * 365,
      severe_enabled = 1
    )
  )
  
  differences <- rep(NA, n_tries)
  EIR_tried <- rep(NA, n_tries)
  EIR_try <- EIR
  
  for(i in 1:n_tries){
    
    simparams <- set_equilibrium(simparams, EIR_try)
    
    #Set up species
    #Set up bionomics
    all_bionomics <- sapply(strsplit(species, ";")[[1]], function(h){
      
      if(nrow(vector_bionomic_data[which(vector_bionomic_data$species == h), ]) == 0){
        vector_bionomic_data <- rbind.fill(vector_bionomic_data,
                                           data.frame(species = h,
                                                      rbind(colMeans(vector_bionomic_data[, 2:ncol(vector_bionomic_data)]))))
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
      
      print(a)
      
      output <- run_simulation(runtime, simparams)
      
      # incidence_by_year <- matrix(output[, grepl("inc_", colnames(output))], ncol = 365, byrow = T)
      prevalence_by_year <- matrix(output[, grepl("n_infections", colnames(output))]/rowSums(output[, c("S_count", "D_count", "A_count", "Tr_count", "U_count")]), ncol = 365, byrow = T)
      
      total_M_all <- output[((nrow(output)-364):nrow(output)), grepl("total_M_", colnames(output))]
      EIR_all <- output[((nrow(output)-364):nrow(output)), grepl("EIR", colnames(output))]
      
      total_M_use <- if(is.null(nrow(total_M_all))) total_M = mean(total_M_all) else rbind(colMeans(total_M_all))
      EIR_use <- if(is.null(nrow(EIR_all))) EIR = mean(EIR_all) else rbind(colMeans(EIR_all))
      
      dfz_save <- data.frame(total_runs = total_runs,
                             human_population = human_population,
                             EIR = EIR_try,
                             itn_cov = itn_cov,
                             irs_cov = irs_cov,
                             species = species,
                             species_abundance = species_abundance,
                             vectordensity_fit = vectordensity_fit,
                             mean_prevalence = mean(prevalence_by_year),
                             total_M = total_M_use/human_population,
                             EIR_use = EIR_use,
                             stringsAsFactors = FALSE)
      
      list(likelihood =  -dnbinom(round(vectordensity_fit, 0), 
                                  1, 
                                  mu = round(as.numeric(as.character(dfz_save$total_M_stephensi)), 0), 
                                  log = TRUE),
           data = dfz_save)
      
    }, simplify = FALSE)
    
    mean_outcome <- do.call(rbind, sapply(1:length(all_run_done), function(x) all_run_done[[x]][[2]], simplify = FALSE))
    mean_outcome_df <- data.frame(species = mean_outcome$species[1],
                                  species_abundance = mean_outcome$species_abundance[1],
                                  rbind(colMeans(mean_outcome[, which(colnames(mean_outcome) == "mean_prevalence"):ncol(mean_outcome)])))
    
    this <- (vectordensity_fit - mean_outcome_df$total_M)
    
    message(paste0("EIR: ", EIR_try, " - ", "Vectordensity difference: ", this))
    EIR_tried[i] <- EIR_try
    
    differences[i] <- this
    
    if(abs(differences[i]) < 0.5){
      break
    } else {
      EIR_try <- pmax(EIR_try + rnorm(1, mean = pmax(5 * this/vectordensity_fit, 5) * sign(this)), 0.0000000000000001)
    }
  }
  
  mean_outcome_df$EIR_tried <- paste(EIR_tried, collapse = ";")
  mean_outcome_df$differences <- paste(differences, collapse = ";")
  
  mean_outcome_df
  
}










