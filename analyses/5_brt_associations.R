#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools);
library(dismo); library(gbm); library(mltools); library(glmnet); library(caret); library(finetune)
library(tidymodels); library(doParallel); library(vip); library(forcats); library(vip); library(xgboost)

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

# Prepping Data - Selecting Variables, Creating Recipe, Generating CV Folds
data <- overall %>%
  dplyr::select(period , population_per_1km:worldclim_9) %>%
  dplyr::select(-contains("LC"))
envt_recipe <- recipe(period  ~ ., data = data) %>% 
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) 
envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
juiced <- bake(envt_prepped, new_data = NULL)

cv_splits <- vfold_cv(data, v = 6) # v sets number of splits
perf_metrics <- metric_set(yardstick::rmse)

# Creating Test and Training Data
set.seed(123)
rf_split <- initial_split(data, prop = 5/6)
rf_train <- training(rf_split)
rf_test <- testing(rf_split)

# Boosted Regression Tree
xgb_spec <- boost_tree(trees = 1000, 
                       tree_depth = tune(), 
                       min_n = tune(),
                       mtry = tune(),
                       sample_size = tune(),
                       learn_rate = tune()) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

# Create Workflow and Grid of Hyperparameter Values to Try
xgb_wf <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(xgb_spec)
set.seed(123)
xgb_grid <- grid_max_entropy(
  tree_depth(c(5L, 10L)),
  min_n(c(10L, 40L)),
  mtry(c(5L, 10L)),
  sample_prop(c(0.5, 1.0)),
  learn_rate(c(-2, -1)),
  size = 150)

# Running All the Hyperparameter Combinations
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels); library(finetune)})
set.seed(234)
xgb_word_rs <- tune_race_anova(
  object = xgb_wf,
  resamples = cv_splits, 
  grid = xgb_grid,
  metrics = perf_metrics,
  control = control_race(verbose_elim = TRUE))
stopCluster(cl)

## add this in at somepoint I think
##    control = control_resamples(save_pred = TRUE)


xgb_word_rs$.metrics[[1]]$

# Selecting the Best Model and Exploring Quality of Fit
show_best(xgb_word_rs, "rmse")
best_rmse <- select_best(xgb_word_rs, "rmse")

xgb_final <- xgb_wf %>%
  finalize_workflow(best_rmse)

fit_call <- xgb_final %>% 
  fit(data = data) 

collect_predictions(xgb_final)
collect_predictions(xgb_wf)
collect_predictions(fit_call)
y <- fit_call %>%
  collect_predictions()

fitted(fit_call)

x <- fit_call %>% 
  predict(new_data = data) 
plot(data$period, x$.pred, pch = 20)

#plot(data$period, fit_call$fit$fit$fit$predictions, pch = 20, col = "orange")

fit_call$fit$actions$model$spec$args

random_forest_final <- rf_workflow %>%
  finalize_workflow(best_rmse)

fit_call <- random_forest_final %>% 
  fit(data = data)

x <- fit_call %>% 
  predict(data) 

plot(data$period, x$.pred, pch = 20)
sqrt(sum((data$period  - x$.pred)^2)/length(data$period ))

points(data$period , fit_call$fit$fit$fit$predictions, pch = 20, col = "orange")



plot_race(xgb_word_rs)

show_best(xgb_word_rs)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

xgb_last <-xgb_wf %>%
  finalize_workflow(select_best(xgb_word_rs, "rmse")) %>%
  last_fit(rf_split)

xgb_last %>%
  collect_metrics()

xgb_last %>%
  collect_predictions()


plot(xgb_last$.predictions[[1]]$period, xgb_last$.predictions[[1]]$.pred) 

extract_workflow(xgb_last) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point", num_features = 15)

