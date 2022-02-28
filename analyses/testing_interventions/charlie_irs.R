library(data.table)
library(devtools)
library(ggplot2)
library(tidyverse)
library(odin)
library(ICDMM)
library(ggpubr)
library(viridis)
library(viridis)
library(tidyverse)
library(RColorBrewer)

options(scipen = 999)

invisible(sapply(list.files("R/functions/", full.names = T, recursive = T), function(x) source(x)))

#Read in the reference csv for matching EIR and density
EIR_prev_reference_df <- fread("//fi--didef3.dide.ic.ac.uk/Malaria/Arran/stephensi_east_africa/output/EIR_LHC_full_run_LHC_giant_v3_1_to_100.csv")

#Loading in pre-defined seasonality for illustration
admin_units_seasonal <- readRDS("//fi--didef3.dide.ic.ac.uk/Malaria/Arran/stephensi_east_africa/data/malaria_mosquito/admin_units_seasonal.rds")
djibouti_seasonality <- subset(admin_units_seasonal, admin1 == "Djibouti")
custom_seasonality <- seasonal_profile(djibouti_seasonality) #Using djibouti seasonality as a default to illustrate, put your vectors in here

match_irs_value <- charlie_irs_function(EIR_prev_reference_df = EIR_prev_reference_df,
                                        years_run = 25, #Number of years to run the simulation for
                                        year_irs_implemented = 15, #The year irs coverage is input
                                        day_irs_implemented = 180, #Accepts values 1:365
                                        irs_half_life = 0.5, #The irs half life
                                        IRS_interval = 1, #How often irs is implemented in years, default is 1
                                        irs_cov = 0.8, #irs coverage
                                        match_prev = 0.05, #what we want the initial prevalence to be - its very rough 
                                        custom_seasonality = custom_seasonality #your seasonal profile of length 365, set to NA if you want no seasonality
)


#You will want to edit out the burn in period from your graphs 
ggplot(data = subset(match_irs_value), 
       aes(x = t, y = value, color = type)) +
  geom_line() +
  facet_wrap(~type, scales = "free_y") +
  labs(x = "Time", y = "Value", color = "Metric") +
  geom_vline(aes(xintercept = (year_irs_implemented * 365) + day_irs_implemented), linetype = "dashed")











