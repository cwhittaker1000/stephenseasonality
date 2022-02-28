# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(scales); library(patchwork); library(here); library(zoo)

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
source(here::here("functions", "time_series_characterisation_functions.R"))
overall <- readRDS(file = here::here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

## Arran says these are current stephensi bionomics - compare to above
# mu0 = 0.12345679 
# Q0 = 0.25
# chi = 0.5 
# bites_Bed = 0.4776
# bites_Indoors = 0.52186
# delayMos = 10

# Generating the vector of how mosquito density changes over time 
density_start <- 20
density_end <- 20
values <- sigmoid(seq(-10, 10, length.out = 365 * 4)) # How long you want introduction to last - here, 10 years
density_vec <- c(rep(density_start, 365 * 1), 
                 pmin((values * (density_end - density_start)) + density_start, density_end),
                 rep(density_end, 365 * 23))
dens_vec_df <- data.frame(vector_density = density_vec, time = seq(1:length(density_vec))/365)
plot(dens_vec_df$vector_density)

# Running Model - Arran's Old Version
set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                        #Model parameters to work out transmission
                                        init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                        #These are the mosquito parameters (bionomics)
                                        Q0 = stephensi_data$Q0, 
                                        chi = stephensi_data$chi, 
                                        bites_Bed = stephensi_data$bites_Bed,
                                        bites_Indoors = stephensi_data$bites_Indoors, 
                                        #These are our seasonal variations
                                        scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                        custom_seasonality = NA, # NA for perennial
                                        #This sets up how long we want to run for and the density vec
                                        time_length = length(density_vec),
                                        density_vec = density_vec,
                                        # Vector control interventions in
                                        ITN_IRS_on = 24 * 365,
                                        num_int = 4,
                                        IRS_interval = 10000,
                                        ITN_interval = 10000,
                                        irs_cov = 0.8,
                                        itn_cov = 0)
set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)

# Running Model With Ellie's Adapted IRS Formulation
set_up_model2 <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                        #Model parameters to work out transmission
                                        init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                        #These are the mosquito parameters (bionomics)
                                        Q0 = stephensi_data$Q0, 
                                        chi = stephensi_data$chi, 
                                        bites_Bed = stephensi_data$bites_Bed,
                                        bites_Indoors = stephensi_data$bites_Indoors, 
                                        #These are our seasonal variations
                                        scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                        custom_seasonality = NA, # NA for perennial
                                        #This sets up how long we want to run for and the density vec
                                        time_length = length(density_vec),
                                        density_vec = density_vec,
                                        # Vector control interventions in
                                        ITN_IRS_on = 24 * 365,
                                        num_int = 4,
                                        irs_decay_det1 = 1,
                                        irs_decay_det2 = 1, 
                                        irs_decay_succ1 = 1,
                                        irs_decay_succ2 = 1, 
                                        irs_decay_mort1 = 1, 
                                        irs_decay_mort2 = 1, 
                                        irs_cov = 0.8,
                                        IRS_interval = 10000,
                                        ITN_interval = 10000,
                                        itn_cov = 0)
set_up_model2 <- set_up_model2$generator(user = set_up_model2$state, use_dde = TRUE)

sum(!(names(set_up_model$contents()) %in% names(set_up_model2$contents())))
sum(!(names(set_up_model2$contents()) %in% names(set_up_model$contents())))
names(set_up_model2$contents())[!(names(set_up_model2$contents()) %in% names(set_up_model$contents()))]
names(set_up_model$contents())[!(names(set_up_model$contents()) %in% names(set_up_model2$contents()))]

x <- set_up_model$contents()
y <- set_up_model2$contents()
shared_names <- names(set_up_model$contents())[(names(set_up_model$contents()) %in% names(set_up_model2$contents()))]
for (i in 1:length(shared_names)) {
  temp1 <- x[[shared_names[i]]]
  temp2 <- y[[shared_names[i]]]
  identical <- identical(temp1, temp2)
  if (!identical) {
    print(i)
  }
  #print(paste0(i, " is ", identical))
}


statex <- set_up_model$state
statey <- set_up_model2$state
shared_names <- names(statex)[names(statex) %in% names(statey)]
for (i in 1:length(shared_names)) {
  temp1 <- statex[[shared_names[i]]]
  temp2 <- statey[[shared_names[i]]]
  identical <- identical(temp1, temp2)
  if (!identical) {
    print(i)
  }
  #print(paste0(i, " is ", identical))
}
names(set_up_model2)






# Running Model - Arran's Old Version
set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                        #Model parameters to work out transmission
                                        init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                        #These are the mosquito parameters (bionomics)
                                        Q0 = stephensi_data$Q0, 
                                        chi = stephensi_data$chi, 
                                        bites_Bed = stephensi_data$bites_Bed,
                                        bites_Indoors = stephensi_data$bites_Indoors, 
                                        #These are our seasonal variations
                                        scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                        custom_seasonality = NA, # NA for perennial
                                        #This sets up how long we want to run for and the density vec
                                        time_length = length(density_vec),
                                        density_vec = density_vec,
                                        # Vector control interventions in
                                        ITN_IRS_on = 24 * 365,
                                        num_int = 4,
                                        IRS_interval = 10000,
                                        ITN_interval = 10000,
                                        irs_cov = 0.2,
                                        itn_cov = 0)
set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
mod_run <- set_up_model$run(t = 1:length(density_vec))
out <- set_up_model$transform_variables(mod_run)
model_ran <- as.data.frame(out)

# Running Model With Ellie's Adapted IRS Formulation
set_up_model2 <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                         #Model parameters to work out transmission
                                         init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                         #These are the mosquito parameters (bionomics)
                                         Q0 = stephensi_data$Q0, 
                                         chi = stephensi_data$chi, 
                                         bites_Bed = stephensi_data$bites_Bed,
                                         bites_Indoors = stephensi_data$bites_Indoors, 
                                         #These are our seasonal variations
                                         scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                         custom_seasonality = NA, # NA for perennial
                                         #This sets up how long we want to run for and the density vec
                                         time_length = length(density_vec),
                                         density_vec = density_vec,
                                         # Vector control interventions in
                                         ITN_IRS_on = 24 * 365,
                                         num_int = 4,
                                         irs_decay_det1 = 1,
                                         irs_decay_det2 = 1, 
                                         irs_decay_succ1 = 1,
                                         irs_decay_succ2 = 1, 
                                         irs_decay_mort1 = 1, 
                                         irs_decay_mort2 = 1, 
                                         irs_cov = 0.2,
                                         IRS_interval = 10000,
                                         ITN_interval = 10000,
                                         itn_cov = 0)
set_up_model2 <- set_up_model2$generator(user = set_up_model2$state, use_dde = TRUE)
mod_run2 <- set_up_model2$run(t = 1:length(density_vec))
out2 <- set_up_model2$transform_variables(mod_run2)
model_ran2 <- as.data.frame(out2)

par(mfrow = c(2, 3))
plot(model_ran$prev, type = "l", ylim = c(0, max(c(model_ran$prev, model_ran2$prev))))
lines(model_ran2$prev, type = "l", col = "red")

plot(model_ran$Incidence, type = "l", ylim = c(0, max(c(model_ran$Incidence, model_ran2$Incidence))))
lines(model_ran2$Incidence, type = "l", col = "red")

plot(model_ran$Incidence[(23*365):(26*365)], type = "l", ylim = c(0, max(c(model_ran$Incidence[(23*365):(26*365)], model_ran2$Incidence[(23*365):(26*365)]))))
lines(model_ran2$Incidence[(23*365):(26*365)], type = "l", col = "red")

########
set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                        #Model parameters to work out transmission
                                        init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                        #These are the mosquito parameters (bionomics)
                                        Q0 = stephensi_data$Q0, 
                                        chi = stephensi_data$chi, 
                                        bites_Bed = stephensi_data$bites_Bed,
                                        bites_Indoors = stephensi_data$bites_Indoors, 
                                        #These are our seasonal variations
                                        scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                        custom_seasonality = NA, # NA for perennial
                                        #This sets up how long we want to run for and the density vec
                                        time_length = length(density_vec),
                                        density_vec = density_vec)

#Process model formulation
set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
mod_run <- set_up_model$run(t = 1:length(density_vec))
out <- set_up_model$transform_variables(mod_run)
model_ran <- as.data.frame(out)

# Running Model With Ellie's Adapted IRS Formulation
set_up_model2 <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
                                        #Model parameters to work out transmission
                                        init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
                                        #These are the mosquito parameters (bionomics)
                                        Q0 = stephensi_data$Q0, 
                                        chi = stephensi_data$chi, 
                                        bites_Bed = stephensi_data$bites_Bed,
                                        bites_Indoors = stephensi_data$bites_Indoors, 
                                        #These are our seasonal variations
                                        scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                        custom_seasonality = NA, # NA for perennial
                                        #This sets up how long we want to run for and the density vec
                                        time_length = length(density_vec),
                                        density_vec = density_vec,
                                        irs_decay_det1 = 1,
                                        irs_decay_det2 = 1, 
                                        irs_decay_succ1 = 1,
                                        irs_decay_succ2 = 1, 
                                        irs_decay_mort1 = 1, 
                                        irs_decay_mort2 = 1)

set_up_model2 <- set_up_model2$generator(user = set_up_model$state, use_dde = TRUE)
mod_run <- set_up_model$run(t = 1:length(density_vec))
out <- set_up_model$transform_variables(mod_run)
model_ran2 <- as.data.frame(out)

# par(mfrow = c(1, 3))
plot(model_ran$prev, type = "l", ylim = c(0, max(c(model_ran$prev, model_ran2$prev))))
lines(model_ran2$prev, type = "l", col = "red")

plot(model_ran$Incidence, type = "l", ylim = c(0, max(c(model_ran$Incidence, model_ran2$Incidence))))
lines(model_ran2$Incidence, type = "l", col = "red")

plot(model_ran$Incidence[(23*365):(26*365)], type = "l", ylim = c(0, max(c(model_ran$Incidence[(23*365):(26*365)], model_ran2$Incidence[(23*365):(26*365)]))))
lines(model_ran2$Incidence[(23*365):(26*365)], type = "l", col = "red")



#####
# x <- model_param_list_create(Q0 = stephensi_data$Q0,
#                              chi = stephensi_data$chi,
#                              bites_Bed = stephensi_data$bites_Bed,
#                              bites_Indoors = stephensi_data$bites_Indoors,
#                              scalar = 1,
#                              custom_seasonality = NA,
#                              time_length = length(density_vec),
#                              density_vec = density_vec,
#                              ITN_IRS_on = 24 * 365,
#                              num_int = 4,
#                              IRS_interval = 10000,
#                              ITN_interval = 10000,
#                              irs_cov = 0.8,
#                              itn_cov = 0)
# 
# y <- model_param_list_create(Q0 = stephensi_data$Q0,
#                              chi = stephensi_data$chi,
#                              bites_Bed = stephensi_data$bites_Bed,
#                              bites_Indoors = stephensi_data$bites_Indoors,
#                              scalar = 1,
#                              custom_seasonality = NA,
#                              time_length = length(density_vec),
#                              density_vec = density_vec,
#                              ITN_IRS_on = 24 * 365,
#                              num_int = 4,
#                              irs_decay_det1 = 1,
#                              irs_decay_det2 = 1,
#                              irs_decay_succ1 = 1,
#                              irs_decay_succ2 = 1,
#                              irs_decay_mort1 = 1,
#                              irs_decay_mort2 = 1,
#                              irs_cov = 0.8,
#                              IRS_interval = 10000,
#                              ITN_interval = 10000,
#                              itn_cov = 0)
# 
# sum(!(names(x) %in% names(y)))
# sum(!(names(y) %in% names(x)))
# names(y)[!(names(y) %in% names(x))]
# 
# shared_names <- names(x)
# for (i in 1:length(shared_names)) {
# 
#   temp1 <- x[[shared_names[i]]]
#   temp2 <- y[[shared_names[i]]]
#   identical <- identical(temp1, temp2)
#   print(paste0(i, " is ", identical))
# 
# }

# set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
#                                         #Model parameters to work out transmission
#                                         init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
#                                         #These are the mosquito parameters (bionomics)
#                                         Q0 = stephensi_data$Q0, 
#                                         chi = stephensi_data$chi, 
#                                         bites_Bed = stephensi_data$bites_Bed,
#                                         bites_Indoors = stephensi_data$bites_Indoors, 
#                                         #These are our seasonal variations
#                                         scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
#                                         custom_seasonality = NA, # NA for perennial
#                                         #This sets up how long we want to run for and the density vec
#                                         time_length = length(density_vec),
#                                         density_vec = density_vec,
#                                         # Vector control interventions in
#                                         ITN_IRS_on = 24 * 365,
#                                         num_int = 4,
#                                         IRS_interval = 10000,
#                                         ITN_interval = 10000,
#                                         irs_cov = 0.8,
#                                         itn_cov = 0)
# set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
# 
# # Running Model With Ellie's Adapted IRS Formulation
# set_up_model2 <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality_new_IRS.R", # model file
#                                          #Model parameters to work out transmission
#                                          init_EIR = 1.5, # initial EIR from which the endemic equilibria solution is created
#                                          #These are the mosquito parameters (bionomics)
#                                          Q0 = stephensi_data$Q0, 
#                                          chi = stephensi_data$chi, 
#                                          bites_Bed = stephensi_data$bites_Bed,
#                                          bites_Indoors = stephensi_data$bites_Indoors, 
#                                          #These are our seasonal variations
#                                          scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
#                                          custom_seasonality = NA, # NA for perennial
#                                          #This sets up how long we want to run for and the density vec
#                                          time_length = length(density_vec),
#                                          density_vec = density_vec,
#                                          # Vector control interventions in
#                                          ITN_IRS_on = 24 * 365,
#                                          num_int = 4,
#                                          irs_decay_det1 = 1,
#                                          irs_decay_det2 = 1, 
#                                          irs_decay_succ1 = 1,
#                                          irs_decay_succ2 = 1, 
#                                          irs_decay_mort1 = 1, 
#                                          irs_decay_mort2 = 1, 
#                                          irs_cov = 0.8,
#                                          IRS_interval = 10000,
#                                          ITN_interval = 10000,
#                                          itn_cov = 0)
# set_up_model2 <- set_up_model2$generator(user = set_up_model2$state, use_dde = TRUE)
# 
# sum(!(names(set_up_model$contents()) %in% names(set_up_model2$contents())))
# sum(!(names(set_up_model2$contents()) %in% names(set_up_model$contents())))
# names(set_up_model2$contents())[!(names(set_up_model2$contents()) %in% names(set_up_model$contents()))]
# names(set_up_model$contents())[!(names(set_up_model$contents()) %in% names(set_up_model2$contents()))]
# 
# shared_names <- names(set_up_model$contents())[(names(set_up_model$contents()) %in% names(set_up_model2$contents()))]
# for (i in 1:length(shared_names)) {
#   temp1 <- x[[shared_names[i]]]
#   temp2 <- y[[shared_names[i]]]
#   identical <- identical(temp1, temp2)
#   if (!identical) {
#     print(i)
#   }
#   #print(paste0(i, " is ", identical))
# }
# 
# mod_run <- set_up_model$run(t = 1:length(density_vec))
# out <- set_up_model$transform_variables(mod_run)
# model_ran <- as.data.frame(out)