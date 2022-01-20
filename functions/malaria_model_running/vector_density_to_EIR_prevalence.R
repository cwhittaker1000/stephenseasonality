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

vector_density_to_EIR_prevalence <- function(id,
                                             row,
                                             human_population = 1000, 
                                             vectordensity_fit, 
                                             species, 
                                             species_abundance, 
                                             itn_cov, 
                                             irs_cov, 
                                             runtime = 365 * 6, 
                                             save_name_fill = ""){
  
  
  sapply(unlist(strsplit(as.character(row), ";")), function(row_use){
    
    message("Row: ", row_use)
    
    #Load in vector bionomic data
    vector_bionomic_data <- read.csv("data/mosquito/processed_species_bionomics_LHC.csv", stringsAsFactors = F)
    vector_bionomic_data <- vector_bionomic_data[which(vector_bionomic_data$run == row_use), ]
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
    
    simparams <- set_equilibrium(simparams, 5)
    simparams <- parameterise_total_M(simparams, vectordensity_fit * human_population)
    
    output <- run_simulation(runtime, simparams)
    plot(output$n_infections)
    plot(output$EIR_stephensi)
    
    # incidence_by_year <- matrix(output[, grepl("inc_", colnames(output))], ncol = 365, byrow = T)
    prevalence_by_year <- matrix(rowSums(output[, c("D_count", "A_count", "Tr_count", "U_count")])/rowSums(output[, c("S_count", "D_count", "A_count", "Tr_count", "U_count")]), ncol = 365, byrow = T)
    
    total_M_all <- output[((nrow(output)- 365 * 3):nrow(output)), grepl("total_M_", colnames(output))]
    EIR_all <- output[((nrow(output)- 365 * 3):nrow(output)), grepl("EIR", colnames(output))]
    
    total_M_use <- if(is.null(nrow(total_M_all))) total_M = mean(total_M_all) else rbind(colMeans(total_M_all))
    EIR_use <- if(is.null(nrow(EIR_all))) EIR = mean(EIR_all) else rbind(colMeans(EIR_all))
    
    inc_use <- output$n_inc_0_36500
  
    output_df <- data.frame(id = id,
                            row = row_use,
                            human_population = human_population,
                            itn_cov = itn_cov,
                            irs_cov = irs_cov,
                            species = species,
                            species_abundance = species_abundance,
                            vectordensity_fit = vectordensity_fit,
                            mean_prevalence = mean(prevalence_by_year[(nrow(prevalence_by_year) - 2) : nrow(prevalence_by_year), ]),
                            annual_incidence_per_1000 = sum(inc_use[(length(inc_use) - (365 * 3) + 1):length(inc_use)])/3,
                            total_M = total_M_use/human_population,
                            EIR_use = EIR_use,
                            stringsAsFactors = FALSE)
    
    savename <- paste(row_use, human_population, itn_cov, irs_cov, id, runtime, sep = "_")
    if(!dir.exists(paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/"))) dir.create(paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/"), recursive = T)
    save_name_fill <- paste0("output/vector_competence/malariasimulation_vectordensity_runs/", row_use, "/", savename, ".csv")
    write.csv(output_df, save_name_fill, row.names = FALSE)
    
  })
  
}
