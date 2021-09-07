#######################################################################################################
##                                                                                                   ##
##                            Loading Required Libraries and Functions                               ##
##                                                                                                   ##
#######################################################################################################
devtools::install_github("mrc-ide/umbrella")
umbrella::download(2017:2019, 1:366, "data/raw")

library(tidyverse)

admin1 <- readRDS(here::here("data/raw/stephensi_countries_admin1.RDS")) 
admin2 <- readRDS(here::here("data/raw/stephensi_countries_admin2.RDS")) 
