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
library(tidymodels); library(doParallel); library(vip); library(forcats); library(vip)

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
juiced <- juice(envt_prepped)

cv_splits <- vfold_cv(data, v = 6) # v sets number of splits
perf_metrics <- metric_set(yardstick::rmse)

# Creating Test and Training Data
set.seed(123)
rf_split <- initial_split(data, prop = 5/6)
rf_train <- training(rf_split)
rf_test <- testing(rf_split)

# Random Forest
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("regression") %>%
  set_engine("ranger")
rf_workflow <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(random_forest)
rf_grid <- grid_regular(mtry(range = c(5, 20)), min_n(range = c(2, 10)), levels = 10)
rf_grid %>% 
  ggplot(aes(x = mtry, y = min_n)) +
  geom_point()

# Running the Model in Parallel 
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels)})
set.seed(345)
tune_res <- tune_grid(object = rf_workflow,
                      resamples = cv_splits,
                      grid = rf_grid,
                      metrics = perf_metrics)
stopCluster(cl)

# Evaluating Performance
tune_res %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "rmse")

# Selecting the Best Model and Exploring Quality of Fit
show_best(tune_res, "rmse")
best_rmse <- select_best(tune_res, "rmse")

# Finalising Workflow and Fitting this Model to the Entire Data 

## WHAT IS THE DIFFERENCE BETWEEN FIT AND PREDICT HERE - WHAT IS GOING ON??
random_forest_final <- rf_workflow %>%
  finalize_workflow(best_rmse) %>% 
  fit(data = data)

x <- random_forest_final %>% 
  predict(data) 

## unclear to me currently why x and random_forest_final produce different predictions - need to work this out
plot(data$period, x$.pred, pch = 20)
points(data$period , random_forest_final$fit$fit$fit$predictions, pch = 20, col = "orange")
# think this is the only point I need to work out, and then the random forest workflow is done

random_forest_final %>% 
  predict(data) %>% 
  bind_cols(dplyr::select(data, period )) %>% 
  perf_metrics(truth = period , estimate = .pred)
sqrt(sum((data$period  - x$.pred)^2)/length(data$period ))
sqrt(sum((data$period  - random_forest_final$fit$fit$fit$predictions)^2)/length(data$period ))
# why is there a significant diffrence in the sum of squares between predictions from x and predictions from random_forest_final

final_rf <- finalize_model(random_forest, best_rmse)
y <- final_rf %>% 
  fit(period  ~ ., data = data) %>%
  predict(data) 
plot(data$period , y$.pred, pch = 20)
sqrt(sum((data$period  - y$.pred)^2)/length(data$period ))

# The below aren't actually different (you can use juiced or unjuiced data) 
# Any diff is because of stochasticity in the permutation used to calculate variable importance
# Become identical if you set same seed for each. 
set.seed(1)
final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(period  ~ ., data = data) %>%
  vip(geom = "point")
set.seed(1)
final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(period  ~ ., data = juiced) %>%
  vip(geom = "point")
unjuiced <- final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(period  ~ ., data = data)
juiced <- final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(period  ~ ., data = juiced) 
plot(unjuiced$fit$predictions, juiced$fit$predictions)

# Finalising workflow, fitting to training and testing data etc
final_wf <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(final_rf)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

final_res <- final_wf %>%
  last_fit(rf_split)
final_res %>%
  collect_metrics()
random_forest_final %>% 
  predict(data) %>% 
  bind_cols(dplyr::select(data, period )) %>% 
  perf_metrics(truth = period , estimate = .pred)

plot(final_res$.predictions[[1]]$period, final_res$.predictions[[1]]$.pred) 

