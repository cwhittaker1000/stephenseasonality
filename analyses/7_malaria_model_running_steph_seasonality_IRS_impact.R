# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(scales); library(patchwork); library(here); library(zoo);
library(tictoc)

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
source(here::here("functions", "time_series_characterisation_functions.R"))
overall <- readRDS(file = here::here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

# Loading and processing different IRS parameters
actellic <- read.csv(here("data", "IRS_parameters", "Actellic_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

bendiocarb <- read.csv(here("data", "IRS_parameters", "Bendiocarb_uncertainty.csv")) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

sumishield <- read.csv(here("data", "IRS_parameters", "Sumishield_uncertainty.csv")) %>%
  dplyr::rename(irs_decay_mort1 = Mort1, irs_decay_mort2 = Mort2, irs_decay_succ1 = Succ1, 
                irs_decay_succ2 = Succ2, irs_decay_det1 = Det1, irs_decay_det2 = Det2) %>%
  pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
  mutate(parameter = factor(parameter)) %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(median = median(value)) %>%
  pivot_wider(names_from = "parameter", values_from = "median")

# Extracting Mean Realisation for Each Time-Series
set.seed(10)
urban_rural <- overall$city
interpolating_points <- 2
mean_realisation <- matrix(nrow = dim(overall)[1], ncol = (12 * interpolating_points + 1))
prior <- "informative"
for (i in 1:length(overall$id)) {
  # Indexing 
  index <- overall$id[i]
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here::here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here::here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
  }  
  
  # Extracting the mean fitted time-series
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  ordered_timepoints <- all_timepoints[order(all_timepoints)]
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  negbinom_intensity_mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  mean_realisation[i, ] <- negbinom_intensity_mean
}
normalised_output <- t(apply(mean_realisation, 1, normalise_total))

# Converting stephensi profiles to be 365 days in length
scaled_year_output <- t(apply(normalised_output, 1, function(x) {
  temp <- approx(x, n = 365)$y
  temp <- 365 * temp/sum(temp)
  temp[temp < 0.1] <- 0.1
  temp }))
steph_seasonality_list <- vector(mode = "list", length = dim(scaled_year_output)[1])
for (i in 1:(dim(scaled_year_output)[1])) {
  steph_seasonality_list[[i]] <- scaled_year_output[i, ]
}
saveRDS(steph_seasonality_list, file = here::here("data/steph_seasonality_list.rds"))

# Setting Average Mosquito Density
density <- 20
years <- 20
density_vec <- rep(density, 365 * years)

# Generating Counterfactuals - I.e. No IRS
all_cores <- parallel::detectCores(logical = FALSE)
cl <- parallel::makeCluster(all_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(tidymodels); library(odin);
  source(here::here("functions", "time_series_characterisation_functions.R"))
  invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
  bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
  stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)
  steph_seasonality_list <- readRDS(here::here("data/steph_seasonality_list.rds"))
  density <- 20
  years <- 20
  density_vec <- rep(density, 365 * years)
  timings <- seq(0, 365, 15)
})

multi_outputs <- parallel::parLapply(cl = cl , 1:length(steph_seasonality_list), function(x){
  
  # Generate Model Formulation  
  set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                            
                                          #Model parameters to work out transmission
                                          init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                          
                                          # Mosquito Bionomics
                                          Q0 = stephensi_data$Q0, 
                                          chi = stephensi_data$chi, 
                                          bites_Bed = stephensi_data$bites_Bed,
                                          bites_Indoors = stephensi_data$bites_Indoors, 
                                          
                                          # Seasonal Variation In Density
                                          custom_seasonality = steph_seasonality_list[[x]], # NA for perennial
                                          time_length = length(density_vec),
                                          density_vec = density_vec,
                                          
                                          # Vector Control Interventions (Not On, Just Dummy Values)
                                          irs_decay_det1 = 1,
                                          irs_decay_det2 = 1, 
                                          irs_decay_succ1 = 1,
                                          irs_decay_succ2 = 1, 
                                          irs_decay_mort1 = 1, 
                                          irs_decay_mort2 = 1)
    
  # Running Model
  set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
  mod_run <- set_up_model$run(t = 1:length(density_vec))
  out <- set_up_model$transform_variables(mod_run)
  model_ran <- as.data.frame(out)
  
  # Storing Relevant Outputs
  output <- data.frame(seasonal_profile = x, t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
  
  # Removing and Garbage Collecting To Avoid Memory Issues
  rm(out)
  rm(mod_run)
  rm(model_ran)
  gc()
  
  # Returning the Output 
  return(output) 
})

saveRDS(multi_outputs, file = "outputs/malaria_model_IRS_running_counterfactuals.rds")
parallel::stopCluster(cl)

# Intervention Timings
timings <- seq(0, 365, 15)

# Running All Seasonal Profiles and All Intervention Timings - ACTELLIC
all_cores <- parallel::detectCores(logical = FALSE)
cl <- parallel::makeCluster(all_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(tidymodels); library(odin);
  source(here::here("functions", "time_series_characterisation_functions.R"))
  invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
  bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
  stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)
  actellic <- read.csv(here::here("data", "IRS_parameters", "Actellic_uncertainty.csv")) %>%
    pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
    mutate(parameter = factor(parameter)) %>%
    dplyr::group_by(parameter) %>%
    dplyr::summarise(median = median(value)) %>%
    pivot_wider(names_from = "parameter", values_from = "median")
  steph_seasonality_list <- readRDS(here::here("data/steph_seasonality_list.rds"))
  density <- 20
  years <- 20
  density_vec <- rep(density, 365 * years)
  timings <- seq(0, 365, 15)
})

multi_outputs <- parallel::parLapply(cl = cl , 1:length(steph_seasonality_list), function(x){

  output <- data.frame(id = -1, timing = -1, t = -1, incidence = -1, prevalence = -1)
  
  for (i in 1:length(timings)) {

    # Generate Model Formulation  
    set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                            
                                            #Model parameters to work out transmission
                                            init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                            
                                            # Mosquito Bionomics
                                            Q0 = stephensi_data$Q0, 
                                            chi = stephensi_data$chi, 
                                            bites_Bed = stephensi_data$bites_Bed,
                                            bites_Indoors = stephensi_data$bites_Indoors, 
                                            
                                            # Seasonal Variation In Density
                                            custom_seasonality = steph_seasonality_list[[x]], # NA for perennial
                                            time_length = length(density_vec),
                                            density_vec = density_vec,
                                            
                                            # Vector Control Interventions
                                            ITN_IRS_on = ((years - 2) * 365) + timings[i],
                                            num_int = 4,
                                            irs_decay_det1 = actellic$irs_decay_det1,
                                            irs_decay_det2 = actellic$irs_decay_det2, 
                                            irs_decay_succ1 = actellic$irs_decay_succ1,
                                            irs_decay_succ2 = actellic$irs_decay_succ2, 
                                            irs_decay_mort1 = actellic$irs_decay_mort1, 
                                            irs_decay_mort2 = actellic$irs_decay_mort2, 
                                            irs_cov = 0.8,
                                            IRS_interval = 10000,
                                            ITN_interval = 10000,
                                            itn_cov = 0)
    
    # Running Model
    set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
    mod_run <- set_up_model$run(t = 1:length(density_vec))
    out <- set_up_model$transform_variables(mod_run)
    model_ran <- as.data.frame(out)
    
    # Storing Relevant Outputs
    temp_output <- data.frame(id = i, timing = timings[i], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
    output <- rbind(output, temp_output)
    
    # Removing and Garbage Collecting To Avoid Memory Issues
    rm(out)
    rm(mod_run)
    rm(model_ran)
    gc()
    
  }
 
  # Returning the Output 
  return(output) 
})
saveRDS(multi_outputs, file = "outputs/malaria_model_IRS_running_actelic.rds")
parallel::stopCluster(cl)

# Running All Seasonal Profiles and All Intervention Timings - BENDIOCARB
all_cores <- parallel::detectCores(logical = FALSE)
cl <- parallel::makeCluster(all_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(tidymodels); library(odin);
  source(here::here("functions", "time_series_characterisation_functions.R"))
  invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
  bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
  stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)
  bendiocarb <- read.csv(here::here("data", "IRS_parameters", "Bendiocarb_uncertainty.csv")) %>%
    pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
    mutate(parameter = factor(parameter)) %>%
    dplyr::group_by(parameter) %>%
    dplyr::summarise(median = median(value)) %>%
    pivot_wider(names_from = "parameter", values_from = "median")
  steph_seasonality_list <- readRDS(here::here("data/steph_seasonality_list.rds"))
  density <- 20
  years <- 20
  density_vec <- rep(density, 365 * years)
  timings <- seq(0, 365, 15)
})

multi_outputs <- parallel::parLapply(cl = cl , 1:length(steph_seasonality_list), function(x){

  output <- data.frame(id = -1, timing = -1, t = -1, incidence = -1, prevalence = -1)
  
  for (i in 1:length(timings)) {

    # Generate Model Formulation  
    set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                            
                                            #Model parameters to work out transmission
                                            init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                            
                                            # Mosquito Bionomics
                                            Q0 = stephensi_data$Q0, 
                                            chi = stephensi_data$chi, 
                                            bites_Bed = stephensi_data$bites_Bed,
                                            bites_Indoors = stephensi_data$bites_Indoors, 
                                            
                                            # Seasonal Variation In Density
                                            custom_seasonality = steph_seasonality_list[[x]], # NA for perennial
                                            time_length = length(density_vec),
                                            density_vec = density_vec,
                                            
                                            # Vector Control Interventions
                                            ITN_IRS_on = ((years - 2) * 365) + timings[i],
                                            num_int = 4,
                                            irs_decay_det1 = bendiocarb$irs_decay_det1,
                                            irs_decay_det2 = bendiocarb$irs_decay_det2, 
                                            irs_decay_succ1 = bendiocarb$irs_decay_succ1,
                                            irs_decay_succ2 = bendiocarb$irs_decay_succ2, 
                                            irs_decay_mort1 = bendiocarb$irs_decay_mort1, 
                                            irs_decay_mort2 = bendiocarb$irs_decay_mort2, 
                                            irs_cov = 0.8,
                                            IRS_interval = 10000,
                                            ITN_interval = 10000,
                                            itn_cov = 0)
    
    # Running Model
    set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
    mod_run <- set_up_model$run(t = 1:length(density_vec))
    out <- set_up_model$transform_variables(mod_run)
    model_ran <- as.data.frame(out)
    
    # Storing Relevant Outputs
    temp_output <- data.frame(id = i, timing = timings[i], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
    output <- rbind(output, temp_output)
    
    # Removing and Garbage Collecting To Avoid Memory Issues
    rm(out)
    rm(mod_run)
    rm(model_ran)
    gc()
    
  }
  
  # Returning the Output 
  return(output) 
})
saveRDS(multi_outputs, file = "outputs/malaria_model_IRS_running_bendiocarb.rds")
parallel::stopCluster(cl)

# Running All Seasonal Profiles and All Intervention Timings - SUMISHIELD
all_cores <- parallel::detectCores(logical = FALSE)
cl <- parallel::makeCluster(all_cores)
doParallel::registerDoParallel(cl)
parallel::clusterEvalQ(cl, {
  library(tidymodels); library(odin);
  source(here::here("functions", "time_series_characterisation_functions.R"))
  invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
  bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
  stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)
  sumishield <- read.csv(here::here("data", "IRS_parameters", "Sumishield_uncertainty.csv")) %>%
    dplyr::rename(irs_decay_mort1 = Mort1, irs_decay_mort2 = Mort2, irs_decay_succ1 = Succ1, 
                  irs_decay_succ2 = Succ2, irs_decay_det1 = Det1, irs_decay_det2 = Det2) %>%
    pivot_longer(cols = irs_decay_mort1:irs_decay_det2, values_to = "value", names_to = "parameter") %>%
    mutate(parameter = factor(parameter)) %>%
    dplyr::group_by(parameter) %>%
    dplyr::summarise(median = median(value)) %>%
    pivot_wider(names_from = "parameter", values_from = "median")
  steph_seasonality_list <- readRDS(here::here("data/steph_seasonality_list.rds"))
  density <- 20
  years <- 20
  density_vec <- rep(density, 365 * years)
  timings <- seq(0, 365, 15)
})

multi_outputs <- parallel::parLapply(cl = cl , 1:length(steph_seasonality_list), function(x){

  output <- data.frame(id = -1, timing = -1, t = -1, incidence = -1, prevalence = -1)
  
  for (i in 1:length(timings)) {

    # Generate Model Formulation  
    set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                            
                                            #Model parameters to work out transmission
                                            init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                            
                                            # Mosquito Bionomics
                                            Q0 = stephensi_data$Q0, 
                                            chi = stephensi_data$chi, 
                                            bites_Bed = stephensi_data$bites_Bed,
                                            bites_Indoors = stephensi_data$bites_Indoors, 
                                            
                                            # Seasonal Variation In Density
                                            custom_seasonality = steph_seasonality_list[[x]], # NA for perennial
                                            time_length = length(density_vec),
                                            density_vec = density_vec,
                                            
                                            # Vector Control Interventions
                                            ITN_IRS_on = ((years - 2) * 365) + timings[i],
                                            num_int = 4,
                                            irs_decay_det1 = sumishield$irs_decay_det1,
                                            irs_decay_det2 = sumishield$irs_decay_det2, 
                                            irs_decay_succ1 = sumishield$irs_decay_succ1,
                                            irs_decay_succ2 = sumishield$irs_decay_succ2, 
                                            irs_decay_mort1 = sumishield$irs_decay_mort1, 
                                            irs_decay_mort2 = sumishield$irs_decay_mort2, 
                                            irs_cov = 0.8,
                                            IRS_interval = 10000,
                                            ITN_interval = 10000,
                                            itn_cov = 0)
    
    # Running Model
    set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
    mod_run <- set_up_model$run(t = 1:length(density_vec))
    out <- set_up_model$transform_variables(mod_run)
    model_ran <- as.data.frame(out)
    
    # Storing Relevant Outputs
    temp_output <- data.frame(id = i, timing = timings[i], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
    output <- rbind(output, temp_output)
    
    # Removing and Garbage Collecting To Avoid Memory Issues
    rm(out)
    rm(mod_run)
    rm(model_ran)
    gc()
    
  }
  
  # Returning the Output 
  return(output) 
})
saveRDS(multi_outputs, file = "outputs/malaria_model_IRS_running_sumishield.rds")

# list_of_lists <- vector(mode = "list", length = 2)
# tic()
# # multi_outputs <- parallel::parLapply(cl = cl , 1:length(steph_seasonality_list), function(x){
# for (i in 1:2) {
#   
#   parallel::clusterEvalQ(cl, {i <- i})
#   multi_outputs <- parallel::parLapply(cl = cl , 1:4, function(x){
#     
#     output <- data.frame(id = -1, timing = -1, t = -1, incidence = -1, prevalence = -1)
# 
#     # Generate Model Formulation  
#     set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
#                                               
#                                             #Model parameters to work out transmission
#                                             init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
#                                             
#                                             # Mosquito Bionomics
#                                             Q0 = stephensi_data$Q0, 
#                                             chi = stephensi_data$chi, 
#                                             bites_Bed = stephensi_data$bites_Bed,
#                                             bites_Indoors = stephensi_data$bites_Indoors, 
#                                             
#                                             # Seasonal Variation In Density
#                                             custom_seasonality = steph_seasonality_list[[i]], # NA for perennial
#                                             time_length = length(density_vec),
#                                             density_vec = density_vec,
#                                             
#                                             # Vector Control Interventions
#                                             ITN_IRS_on = ((years - 2) * 365) + timings[x],
#                                             num_int = 4,
#                                             irs_decay_det1 = actellic$irs_decay_det1,
#                                             irs_decay_det2 = actellic$irs_decay_det2, 
#                                             irs_decay_succ1 = actellic$irs_decay_succ1,
#                                             irs_decay_succ2 = actellic$irs_decay_succ2, 
#                                             irs_decay_mort1 = actellic$irs_decay_mort1, 
#                                             irs_decay_mort2 = actellic$irs_decay_mort2, 
#                                             irs_cov = 0.8,
#                                             IRS_interval = 10000,
#                                             ITN_interval = 10000,
#                                             itn_cov = 0)
#       
#     # Running Model
#     set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
#     mod_run <- set_up_model$run(t = 1:length(density_vec))
#     out <- set_up_model$transform_variables(mod_run)
#     model_ran <- as.data.frame(out)
#     
#     # Storing Relevant Outputs
#     temp_output <- data.frame(id = x, timing = timings[i], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
# 
#     # Removing and Garbage Collecting To Avoid Memory Issues
#     rm(out)
#     rm(mod_run)
#     rm(model_ran)
#     gc()
#     
#     print(x)
#     
#     # Returning the Output 
#     return(temp_output) 
#       
#   })
#   
#   list_of_lists[[i]] <- multi_outputs
# }
# toc()
# tic()
# multi_outputs <- parallel::parLapply(cl = cl , 1:length(timings), function(x){
#   
#   # Generate Model Formulation  
#   set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
#                                           
#                                           #Model parameters to work out transmission
#                                           init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
#                                           
#                                           # Mosquito Bionomics
#                                           Q0 = stephensi_data$Q0, 
#                                           chi = stephensi_data$chi, 
#                                           bites_Bed = stephensi_data$bites_Bed,
#                                           bites_Indoors = stephensi_data$bites_Indoors, 
#                                           
#                                           # Seasonal Variation In Density
#                                           custom_seasonality = steph_seasonality_list[[1]], # NA for perennial
#                                           time_length = length(density_vec),
#                                           density_vec = density_vec,
#                                           
#                                           # Vector Control Interventions
#                                           ITN_IRS_on = ((years - 2) * 365) + timings[x],
#                                           num_int = 4,
#                                           irs_decay_det1 = actellic$irs_decay_det1,
#                                           irs_decay_det2 = actellic$irs_decay_det2, 
#                                           irs_decay_succ1 = actellic$irs_decay_succ1,
#                                           irs_decay_succ2 = actellic$irs_decay_succ2, 
#                                           irs_decay_mort1 = actellic$irs_decay_mort1, 
#                                           irs_decay_mort2 = actellic$irs_decay_mort2, 
#                                           irs_cov = 0.8,
#                                           IRS_interval = 10000,
#                                           ITN_interval = 10000,
#                                           itn_cov = 0)
#   
#   # Running Model
#   set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
#   mod_run <- set_up_model$run(t = 1:length(density_vec))
#   out <- set_up_model$transform_variables(mod_run)
#   model_ran <- as.data.frame(out)
#   
#   # Storing Relevant Outputs
#   output <- data.frame(id = x, timing = timings[x], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
#   
#   # Removing and Garbage Collecting To Avoid Memory Issues
#   rm(out)
#   rm(mod_run)
#   rm(model_ran)
#   gc()
#   print(x)
#   
#   # Returning the Output 
#   return(output) 
# })
# toc()

# tic()
# multi_outputs <- lapply(1:length(timings[1:8]), function(x){
#     
#   # Generate Model Formulation  
#   set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
#                                           
#                                           #Model parameters to work out transmission
#                                           init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
#                                           
#                                           # Mosquito Bionomics
#                                           Q0 = stephensi_data$Q0, 
#                                           chi = stephensi_data$chi, 
#                                           bites_Bed = stephensi_data$bites_Bed,
#                                           bites_Indoors = stephensi_data$bites_Indoors, 
#                                           
#                                           # Seasonal Variation In Density
#                                           custom_seasonality = steph_seasonality_list[[1]], # NA for perennial
#                                           time_length = length(density_vec),
#                                           density_vec = density_vec,
#                                           
#                                           # Vector Control Interventions
#                                           ITN_IRS_on = ((years - 2) * 365) + timings[x],
#                                           num_int = 4,
#                                           irs_decay_det1 = actellic$irs_decay_det1,
#                                           irs_decay_det2 = actellic$irs_decay_det2, 
#                                           irs_decay_succ1 = actellic$irs_decay_succ1,
#                                           irs_decay_succ2 = actellic$irs_decay_succ2, 
#                                           irs_decay_mort1 = actellic$irs_decay_mort1, 
#                                           irs_decay_mort2 = actellic$irs_decay_mort2, 
#                                           irs_cov = 0.8,
#                                           IRS_interval = 10000,
#                                           ITN_interval = 10000,
#                                           itn_cov = 0)
#   
#   # Running Model
#   set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
#   mod_run <- set_up_model$run(t = 1:length(density_vec))
#   out <- set_up_model$transform_variables(mod_run)
#   model_ran <- as.data.frame(out)
#   
#   # Storing Relevant Outputs
#   output <- data.frame(id = x, timing = timings[x], t = model_ran$t, incidence = model_ran$Incidence, prevalence = model_ran$prev)
#   
#   # Removing and Garbage Collecting To Avoid Memory Issues
#   rm(out)
#   rm(mod_run)
#   rm(model_ran)
#   gc()
#   print(x)
#   
#   # Returning the Output 
#   return(output) 
# })
# toc()
