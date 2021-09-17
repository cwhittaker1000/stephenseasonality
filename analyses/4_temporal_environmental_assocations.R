#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools);
library(dismo); library(gbm); library(mltools); library(glmnet); library(caret);

source(here("functions", "time_series_characterisation_functions.R"))

#######################################################################################################
##                                                                                                   ##
##          Loading In Collated Environmental Data, and Time Series Temporal Properties              ##
##                                                                                                   ##
#######################################################################################################
ts_metadata <- readRDS(here("data", "processed", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "processed", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:worldclim_9, ~ mean(.x, na.rm = TRUE)))
  
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2")) 

subset <- overall %>%
  dplyr::select(id, entropy:weight, population_per_1km, EVI, worldclim_1:worldclim_9) %>%
  mutate(peaks = ifelse(peaks == 1, 1, 0))

ecological_variables <- subset %>%
  dplyr::select(population_per_1km:worldclim_9)
heatmap(cor(ecological_variables))
mean(cor(ecological_variables))

normalised_ecological_variables <- scale(ecological_variables)
summary(prcomp(normalised_ecological_variables))

cor(ecological_variables[, 1:4])

colnames(subset)

c("EVI", "worldclim_1", "population_per_1km", )


set.seed(42)
hit_elnet = train(
  prop_points ~ EVI + worldclim_9 + population_per_1km, 
  data = subset,
  preProcess = c("center", "scale"),
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 5),
  tuneLength = 3
)
hit_elnet$finalModel$beta
hit_elnet$finalModel$a0

coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)

coef(hit_elnet$finalModel, 0.001)

coef.glmnet(hit_elnet$finalModel)

?coef.glmnet

hit_elnet$bestTune$lambda
hit_elnet$bestTune$alpha

x[1, ]

, hit_elnet$bestTune$alpha)

coef(hit_elnet$finalModel)

tree_complexity <- 2
bag_fraction <- 0.95
max_trees <- 25000
learning_rate <- 0.001
brt <- gbm.step(data = subset, 
                gbm.x = 9:ncol(subset), 
                gbm.y = 6, 
                family = "bernoulli", # "gaussian" 
                tree.complexity = tree_complexity, 
                learning.rate = learning_rate, 
                bag.fraction = bag_fraction, 
                max.trees = max_trees, 
                n.folds = 10)
x <- summary(brt)
x

# prediction_table <- model_inputs %>%
#   filter_at(vars(-hospital_beds), all_vars(!is.na(.)))
# predicted <- predict.gbm(world_bank_brt, prediction_table[, c("maternal_mortality", "electricity", "prop_pop", "school_ratio", "rural", "domestic_GDP", "infant_mortality", "school_enrollment", "region_1", "region_2", "region_3", "region_4", "region_6", "region_7", "income_1", "income_2", "income_3", "income_4")],
#                          n.trees = world_bank_brt$gbm.call$best.trees, type = "response")
# prediction_table$pred_beds <- predicted
# observed_vs_predicted <- prediction_table %>%
#   mutate(income_group = factor(income_group, levels = c("Low income", "Lower middle income", "Upper middle income", "High income"))) %>%
#   filter(!is.na(hospital_beds)) %>%
#   select(income_group, hospital_beds, pred_beds)
# 
# a <- ggplot(observed_vs_predicted, aes(colour = income_group, x = hospital_beds, y = pred_beds)) +
#   scale_colour_manual(labels = c("Low Income", "Lower Middle Income", "Upper Middle Income", "High Income"),
#                       values = c("#B7C0EE", "#7067CF", "#362E91", "#241F60")) + 
#   guides(colour = guide_legend(override.aes = list(size = 4))) +
#   geom_point(stat = "identity", size = 2) +
#   geom_segment(x = 0, y = 0, xend = 13, yend = 13, colour = "black", size = 1, linetype = 2) +
#   xlim(c(0, 13)) +
#   ylim(c(0, 13)) +
#   theme_bw() +
#   theme(legend.position = c(0.3, 0.8), legend.title = element_blank(), 
#         axis.title.y = element_text(size = 13.5, vjust = -15), axis.title.x = element_text(size = 13.5, vjust = +0),
#         plot.margin = margin(0.2, 0.5, 0.5, 0.5, "cm"), legend.text = element_text(size = 12),
#         legend.background = element_rect(fill = NA), axis.text = element_text(size = 14, face = "bold")) +
#   labs(x = "Observed Hospital Beds /1000 People", y = "Predicted Hospital Beds\n /1000 People") 
# 
