#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools);
library(dismo); library(gbm); library(mltools); library(glmnet); library(caret); library(themis)
library(tidymodels); library(doParallel); library(vip); library(forcats); library(vip);
library(RColorBrewer); library(corrplot); library(DALEXtra)

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))
source(here("functions", "gp_fitting_functions.R"))

# Loading in extracted stephensi, processed but unsmoothed data and renaming variables
overall <- readRDS(file = here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
urban_rural <- overall$city

features_df <- readRDS(file = here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
features <- features_df[, 7:dim(features_df)[2]]

cluster_output <- readRDS(file = here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_output$cluster

# Exploring and Assessing Urban Association With Unimodal Dynamics

## Rural/urban associations with cluster membership
tab <- table(cluster_membership, urban_rural)
tab <- tab[, 2:3]
chisq.test(tab)

## Rural/urban association with number of peaks
tab <- table(features[, "peaks"], urban_rural)
chisq.test(tab[1:2, 2:3])

## Average Incidence Per 4 Months
table(features[, "peaks"], urban_rural)

mean(features$per_ind_4_months[urban_rural == "Urban"]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 1]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Urban" & features[, "peaks"] == 2]) # almost all urban time series are unimodal

mean(features$per_ind_4_months[urban_rural == "Rural"]) # mixture of unimodal and bimodal time series for rural
mean(features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 1]) # almost all urban time series are unimodal
mean(features$per_ind_4_months[urban_rural == "Rural" & features[, "peaks"] == 2]) # almost all urban time series are unimodal

# Need to explore the rural 1 vs 2 peaks in more depth - what's up with that and why are we seeing such
# divergent dynamics across those 2 groups.