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

# Creating Test and Training Data
set.seed(123)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(period , population_per_1km:worldclim_9) %>%
  dplyr::select(-contains("LC"))
rf_split <- initial_split(data, prop = 5/6)
rf_train <- training(rf_split)
rf_test <- testing(rf_split)

# Prepping Data - Selecting Variables, Creating Recipe, Generating CV Folds
envt_recipe <- recipe(period  ~ ., data = rf_train) %>% 
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) 
envt_prepped <- prep(envt_recipe, training = rf_train, verbose = TRUE)
juiced <- juice(envt_prepped)

cv_splits <- vfold_cv(rf_train, v = 6) # v sets number of splits
perf_metrics <- metric_set(yardstick::rmse)

# Random Forest
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("regression") %>%
  set_engine("ranger", importance = "permutation")
rf_workflow <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(random_forest)
rf_grid <- grid_regular(mtry(range = c(3, 20)), min_n(range = c(2, 10)), levels = 15)
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
                      metrics = perf_metrics,
                      control = control_resamples(save_pred = TRUE))
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

# Extracting and plotting raw predictions
raw_predictions <- tune_res %>%
  collect_predictions()
raw_pred_best_hyp <- raw_predictions %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n)
order_raw_pred_best_hyp <- raw_pred_best_hyp %>%
  dplyr::arrange(., .row)

plot(order_raw_pred_best_hyp$.pred, rf_train$period)
sqrt(sum((rf_train$period  - order_raw_pred_best_hyp$.pred)^2)/length(rf_train$period ))

# Finalising Workflow and Fitting this Model to the Entire Data 
random_forest_final <- rf_workflow %>%
  finalize_workflow(best_rmse)

# checking that this gives us the OOB error model 
x <- random_forest_final %>%
  fit(data = juiced)
plot(x$fit$fit$fit$predictions, rf_train$period)
sqrt(sum((rf_train$period  - x$fit$fit$fit$predictions)^2)/length(rf_train$period ))

# calculating variable importance
set.seed(1)
random_forest_final %>%
  fit(data = juiced) %>%
  pull_workflow_fit() %>%
  vip(geom = "point")

# alternative way of calculating variable importance (confirm identical to above with same seed)
final_rf <- finalize_model(random_forest, best_rmse)
set.seed(1)
final_rf %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(period  ~ ., data = rf_train) %>%
  vip(geom = "point")

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

final_res <- random_forest_final %>%
  last_fit(rf_split)

set.seed(1)
extract_workflow(final_res) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")

final_res %>%
  collect_metrics()
plot(final_res$.predictions[[1]]$period, final_res$.predictions[[1]]$.pred, 
     xlim = c(5, 16), ylim = c(5, 16)) 
sqrt(sum((final_res$.predictions[[1]]$period - final_res$.predictions[[1]]$.pred)^2/length(final_res$.predictions[[1]]$.pred)))










# below is wrong because once you've finalised aim is to just evaluate once on the validation data set
fit_call <- random_forest_final %>% 
  fit(data = rf_train)
predictions <- fit_call %>% 
  predict(rf_train) 
sqrt(sum((rf_train$period  - predictions$.pred)^2)/length(data$period)) # predictions passing the data down every tree
sqrt(sum((rf_train$period  - fit_call$fit$fit$fit$predictions)^2)/length(data$period )) # oob predictions

fit_call %>% 
  collect_metrics()

  predict(rf_train) %>% 
  bind_cols(dplyr::select(rf_train, period )) %>% 
  perf_metrics(truth = period , estimate = .pred)

plot(rf_train$period, predictions$.pred, pch = 20)
points(rf_train$period , fit_call$fit$fit$fit$predictions, pch = 20, col = "orange")



sqrt(sum((data$period  - x$.pred)^2)/length(data$period ))
sqrt(sum((rf_train$period  - fit_call$fit$fit$fit$predictions)^2)/length(data$period ))

final_rf <- finalize_model(random_forest, best_rmse)
y <- final_rf %>% 
  fit(period  ~ ., data = data) %>%
  predict(data) 
plot(data$period , y$.pred, pch = 20)
points(data$period, x$.pred, pch = 20, col = "red")
points(data$period , fit_call$fit$fit$fit$predictions, pch = 20, col = "orange")
sqrt(sum((data$period  - y$.pred)^2)/length(data$period ))

# The below aren't actually different (you can use juiced or unjuiced data) 
# Any diff is because of stochasticity in the permutation used to calculate variable importance
# Become identical if you set same seed for each. 
# BUT THESE EVALUATIONS OF VARIABLE IMPORTANCE ARE USING THE WORSE MODEL AND I DON'T KNOW
# 1) WHY, AND 2) WHETHER THAT'S IMPORTANT.
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

plot(data$period , y$.pred, pch = 20)
points(data$period, x$.pred, pch = 20, col = "red")
sqrt(sum((data$period  - y$.pred)^2)/length(data$period ))
sqrt(sum((data$period  - x$.pred)^2)/length(data$period ))

plot(data$period , fit_call$fit$fit$fit$predictions, pch = 20, col = "orange")
points(data$period, unjuiced$fit$predictions, pch = 20, col = "blue")
points(data$period, juiced$fit$predictions, pch = 20, col = "purple")
sqrt(sum((data$period  - fit_call$fit$fit$fit$predictions)^2)/length(data$period ))
sqrt(sum((data$period  - unjuiced$fit$predictions)^2)/length(data$period ))
sqrt(sum((data$period  - juiced$fit$predictions)^2)/length(data$period ))

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
fit_call %>% 
  predict(data) %>% 
  bind_cols(dplyr::select(data, period )) %>% 
  perf_metrics(truth = period , estimate = .pred)

plot(final_res$.predictions[[1]]$period, final_res$.predictions[[1]]$.pred, 
     xlim = c(5, 16), ylim = c(5, 16)) 

# Scrap Code 
# for this to work need to do something like:
#    set_engine("ranger", importance = "impurity")
# ranger_obj <- extract_fit_parsnip(random_forest_final)$fit
# ranger_obj
# ranger_obj$variable.importance

# rf_mod <- # alt way of doing parallel apparently
#   rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>% 
#   set_engine("ranger", num.threads = cores) %>% 
#   set_mode("regression")

#rf_grid <- grid_max_entropy(mtry(range = c(3, 20)), min_n(range = c(2, 10)), size = 170)

## WHAT IS THE DIFFERENCE BETWEEN FIT AND PREDICT HERE - WHAT IS GOING ON??
## why is there a significant diffrence in the sum of squares between predictions from x and predictions from random_forest_final
## unclear to me currently why x and random_forest_final produce different predictions - need to work this out
## think this is the only point I need to work out, and then the random forest workflow is done

## so I *think* the one using fit using the exact model fitted using CV (and which has an approx RMSE of 2.28)
## and the other one is using the model fit to *all* the data? But unclear what's happening here then - is an
## entirely new model being fitted (with different trees but same hyperparameters)? Surely not. But if RMSE decreases
## when predicting on the same dataset (just with bit missing from)
