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

# Creating Test and Training Data
set.seed(123)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(period, country, population_per_1km:worldclim_9) %>%
  dplyr::select(-contains("LC")) %>%
  mutate(country = as.factor(country))
rf_split <- initial_split(data, prop = 5/6, strata = country)
rf_train <- training(rf_split)
rf_test <- testing(rf_split)

# Prepping Data - Selecting Variables, Creating Recipe, Generating CV Folds
envt_recipe <- recipe(period  ~ ., data = rf_train) %>% 
  step_center(all_predictors(), -country) %>% 
  step_scale(all_predictors(), -country) %>%
  step_dummy(country)
envt_prepped <- prep(envt_recipe, training = rf_train, verbose = TRUE)
juiced <- juice(envt_prepped)

cv_splits <- vfold_cv(rf_train, v = 6) # v sets number of splits
perf_metrics <- metric_set(yardstick::rmse)

# Boosted Regression Tree
xgb_spec <- boost_tree(trees = 1000, 
                       tree_depth = tune(), 
                       min_n = tune(),
                       mtry = tune(),
                       sample_size = tune(),
                       learn_rate = tune()) %>%
  set_engine("xgboost", importance = "permutation") %>%
  set_mode("regression")

# Create Workflow and Grid of Hyperparameter Values to Try
xgb_wf <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(xgb_spec)
set.seed(123)
xgb_grid <- grid_max_entropy(
  tree_depth(c(5L, 15L)), #10
  min_n(c(5L, 20L)),
  mtry(c(5L, 20)), #10
  sample_prop(c(0.5, 1.0)),
  learn_rate(c(-2, -1)),
  size = 450)

# CURRENT POINT OF CONFUSION: 
# -> When I train with rf_train (and do cv giving approx 44/9 test/train for hyperparameter optimisation),
#    I get 2.3ish best rmse. When I train that model again on the full rf_train, I get lower (about 1.7).
#    And I get something even lower for rf_test (not that much lower though, 1.5-1.6). 
#    I had rationalised this by saying "well, training on rf_train has given me more datapoints to learn with, 
#    so of course that's going to improve performance on rf_train. And I've avoided overfitting, hence rf_test
#    is ballpark the same as rf_train". 
# -> I did overall find this confusing a bit confusing though, so went back and instead trained with the full dataset 
#    (approx the size of rf_train for each hyperparam optimisation i.e. 53/12). My thinking was that because I'm now
#    training on about 53 samples, I should get an rmse of about 1.7, to match what I got when I took the tuned model
#    in the example above, and fit it to about the same number of data points at the final point in the above example.
#    I *DID* get a reduction in rmse, but not by much - i.e. 2.1ish (compared to 2.3ish for above example) - not as big a 
#    drop (to e.g. 1.7ish) as I was expecting.
# -> AM I MISSING SOMETHING? IS THERE SOME WEIRD OUT OF BAG/OUT OF SAMPLE FIT/PREDICT WEIRDNESS GOING ON FOR THE BRTS 
#    LIKE WITH THE RANDOM FORESTS? OR IS THIS STOCHASTICITY? OR IS THERE SOMETHING ELSE GOING ON? 

# Running All the Hyperparameter Combinations
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels); library(finetune)})
set.seed(234)
xgb_word_rs <- tune_grid(object = xgb_wf,
                         resamples = cv_splits, 
                         grid = xgb_grid,
                         metrics = perf_metrics,
                         control = control_resamples(save_pred = TRUE))
stopCluster(cl)

xgb_word_rs %>%
  collect_metrics() %>%
  filter(.metric == "rmse") %>%
  dplyr::select(mean, mtry:sample_size ) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter") %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rmse")

show_best(xgb_word_rs, "rmse")

best_rmse <- select_best(xgb_word_rs, "rmse")
best_rmse

best_sets <- show_best(xgb_word_rs, "rmse")

oos_predictions <- xgb_word_rs$.predictions %>%
  bind_rows %>%
  filter(mtry == best_rmse$mtry,
         min_n == best_rmse$min_n,
         tree_depth == best_rmse$tree_depth,
         learn_rate == best_rmse$learn_rate,
         sample_size == best_rmse$sample_size)

rmse_vec <- c()
sizes <- c()
for (i in 1:length(xgb_word_rs$.predictions)) {
  temp <- xgb_word_rs$.predictions[[i]] %>%
    filter(mtry == best_rmse$mtry,
           min_n == best_rmse$min_n,
           tree_depth == best_rmse$tree_depth,
           learn_rate == best_rmse$learn_rate,
           sample_size == best_rmse$sample_size)
  temp2 <- sqrt(sum((temp$.pred - temp$period)^2)/length(temp$period))
  rmse_vec <- c(rmse_vec, temp2)
  sizes <- c(sizes, length(temp$period))
}
mean(rmse_vec)
sizes

# this rmse calculation is NOT the same as reported in show_best. But as shown above, that's
# because this is calculating average of each of the individual predictions' errors, rather 
# average of the average error for each fold
plot(oos_predictions$.pred, oos_predictions$period, pch = 20)
sqrt(sum((oos_predictions$.pred - oos_predictions$period)^2)/length(oos_predictions$period))

# Selecting the Best Model and Exploring Quality of Fit
show_best(xgb_word_rs, "rmse")
best_rmse <- select_best(xgb_word_rs, "rmse")

xgb_final <- xgb_wf %>%
  finalize_workflow(best_rmse)

fit_call <- xgb_final %>% 
  fit(data = rf_train) 

x <- fit_call %>% 
  predict(new_data = rf_train) 
plot(rf_train$period, x$.pred, pch = 20)
sqrt(sum((rf_train$period - x$.pred)^2)/length(x$.pred))

x <- fit_call %>% 
  predict(new_data = rf_test) 
plot(rf_test$period, x$.pred, pch = 20)
sqrt(sum((rf_test$period  - x$.pred)^2)/length(rf_test$period ))

final_res <- xgb_final %>%
  last_fit(rf_split) 
set.seed(1)
extract_workflow(final_res) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")

final_res %>%
  collect_metrics()
plot(final_res$.predictions[[1]]$period, final_res$.predictions[[1]]$.pred, xlim = c(5, 16), ylim = c(5, 16)) 
sqrt(sum((final_res$.predictions[[1]]$period - final_res$.predictions[[1]]$.pred)^2/length(final_res$.predictions[[1]]$.pred)))

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

