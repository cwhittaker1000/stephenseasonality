# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc)
library(ipred); library(ICDMM); #devtools::install_github("https://github.com/jhellewell14/ICDMM", dependencies = TRUE)

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

# Generating the vector of how mosquito density changes over time 
density_start <- 0.1
density_end <- 10
values <- sigmoid(seq(-10, 10, length.out = 365 * 10)) # How long you want introduction to last - here, 10 years
density_vec <- c(rep(density_start, 365 * 1), 
                 pmin((values * (density_end - density_start)) + density_start, density_end),
                 rep(density_end, 365 * 10))
plot(density_vec, bty = "l", ylab = "Vector density", xlab = "Time (days)")
abline(v = seq(365, length(density_vec), by = 365), lty = 2)

#Pre-loaded seasonal profiles
admin_units_seasonal <- readRDS("data/admin_units_seasonal.rds")
admin_matches <- admin_match(admin_unit = "Gujarat", country = "India", admin_units_seasonal = admin_units_seasonal)
loc_seasonality <- seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 1)
par(mfrow = c(1, 1))
plot(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 1), ylab = "Seasonality factor", xlab = "Time (days)", bty = "l", type = "l", col = "red", lwd = 2)
lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 0.5), col = "blue", lwd = 2) # scalar controls how flat
lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 0.1), col = "black", lwd = 2)

# Running the deterministic malaria model with increasing vector density over time and with/without seasonality 
two_outputs <- sapply(1:2, function(x){
  
  set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                          #Model parameters to work out transmission
                                          init_EIR = 0.000001, # initial EIR from which the endemic equilibria solution is created
                                          #These are the mosquito parameters (bionomics)
                                          Q0 = stephensi_data$Q0, 
                                          chi = stephensi_data$chi, 
                                          bites_Bed = stephensi_data$bites_Bed,
                                          bites_Indoors = stephensi_data$bites_Indoors, 
                                          #These are our seasonal variations
                                          scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                          custom_seasonality = if(x == 1) loc_seasonality else NA,
                                          #This sets up how long we want to run for and the density vec
                                          time_length = length(density_vec),
                                          density_vec = density_vec)
  
  #Process model formulation
  set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
  mod_run <- set_up_model$run(t = 1:length(density_vec))
  out <- set_up_model$transform_variables(mod_run)
  model_ran <- as.data.frame(out)
  model_ran
  
}, simplify = FALSE)

#Plot prev vs vector density
par(mfrow = c(3, 1))
plot(two_outputs[[1]]$prev, type = "l", bty = "l", ylab = "Prevalence", xlab = "Time (days)", col = "red", ylim = c(0, max(c(two_outputs[[1]]$prev, two_outputs[[2]]$prev))))
lines(two_outputs[[2]]$prev, col = "blue")

plot(two_outputs[[1]]$mv, type = "l", bty = "l", ylab = "Vector density", xlab = "Time (days)", col = "red", ylim = c(0, max(c(two_outputs[[1]]$mv, two_outputs[[2]]$mv))))
lines(two_outputs[[2]]$mv, col = "blue")

plot(colMeans(matrix(two_outputs[[1]]$mv, nrow = 365)), col = "red", lwd = 2, type = "l", bty = "l", ylab = "Mean annual vector density", xlab = "Time (years)")
lines(colMeans(matrix(two_outputs[[2]]$mv, nrow = 365)), col = "blue", lwd = 2)


#Run model to show extinction thresholds

intervention_seasonality <- sapply(seq(0, 1, by = 0.25), function(x){
  
  set_up_model <- create_r_model_epidemic(odin_model_path = "odin_model_seasonality.R",
                                          #Model parameters to work out transmission
                                          init_EIR = 0.001,
                                          larvicide_vector = "0.4;0.5;0.9", #This reduces adult emergence by the value specified
                                          t_larvicide_compact = "0;4000;5000",
                                          #These are the mosquito parameters
                                          Q0 = stephensi_data$Q0, 
                                          chi = stephensi_data$chi, 
                                          bites_Bed = stephensi_data$bites_Bed,
                                          bites_Indoors = stephensi_data$bites_Indoors, 
                                          #These are our seasonal variations
                                          scalar = 2, #This is the scalar shown above in the plots
                                          custom_seasonality = NA, #loc_seasonality,
                                          #This sets up how long we want to run for and the density vec
                                          time_length = length(density_vec),
                                          density_vec = density_vec,
                                          #Need to specify this if you want interventions
                                          num_int = 2)
  
  #Process model formulation
  set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
  mod_run <- set_up_model$run(t = 1:length(density_vec))
  out <- set_up_model$transform_variables(mod_run)
  model_ran <- as.data.frame(out)
  model_ran
  
}, simplify = FALSE)

par(mfrow = c(1, 1))
plot(intervention_seasonality[[1]]$mv, type = "l", bty = "l", ylab = "Prevalence", xlab = "Time (days)", col = "red",
     ylim = c(0, max(c(intervention_seasonality[[1]]$mv, intervention_seasonality[[2]]$mv))))
lines(intervention_seasonality[[2]]$mv, col = "blue")
lines(intervention_seasonality[[3]]$mv, col = "green")
lines(intervention_seasonality[[4]]$mv, col = "orange")
lines(intervention_seasonality[[5]]$mv, col = "purple")

