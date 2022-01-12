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

source(here("functions", "time_series_characterisation_functions.R"))

# function to close connections made by doing things in parallel, which can cause errors elsewhere
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#######################################################################################################
##                                                                                                   ##
##          Loading In Collated Environmental Data, and Time Series Temporal Properties              ##
##                                                                                                   ##
#######################################################################################################
ts_metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE)))
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))

#######################################################################################################
##                                                                                                   ##
##                      Setting Up All The Components Required to Run the Model                      ##
##             I.e. Test/Train Split, CV Folds, Random Forest Engine, Performance Metrics            ##
##                                                                                                   ##
#######################################################################################################
# Subsetting Outcome and Variables for Analysis + Log_10'ing Population (Easier for Visualisation Later On)
set.seed(234)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(peaks, country, population_per_1km:mean_temperature_driest_quarter, -LC_190, -LC_210, -mean_temperature_driest_quarter) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country_peaks = paste0(peaks, "_", country))
data$peaks <- ifelse(data$peaks == 1, "one", "two")
data$peaks <- as.factor(data$peaks)
data$population_per_1km <- log(data$population_per_1km)

# Prepping Data - Selecting Variables, Creating Recipe - No Upsampling
set.seed(915)
envt_recipe <- recipe(peaks  ~ ., data = data) %>%
  update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
  step_center(all_numeric()) %>% 
  step_scale(all_numeric()) %>%
  step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
  step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
  step_dummy(country)
envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
juiced <- juice(envt_prepped)
colnames(juiced)
dim(juiced)
table(juiced$peaks)
new_envt_recipe <- recipe(peaks ~ ., data = juiced) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
  update_role(country_peaks, new_role = "ID")           # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
# Prepping Data - Selecting Variables, Creating Recipe - Upsampling to Balance Amount of Data for 1 and 2 Peaks
set.seed(915)
envt_recipe_ups <- recipe(peaks  ~ ., data = data) %>%
  update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
  step_center(all_numeric()) %>% 
  step_scale(all_numeric()) %>%
  step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
  step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
  step_dummy(country) %>%
  step_smote(peaks)
envt_prepped_ups <- prep(envt_recipe_ups, training = data, verbose = TRUE)
juiced_ups <- juice(envt_prepped_ups)
colnames(juiced_ups)
dim(juiced_ups)
table(juiced_ups$peaks)
new_envt_recipe_ups <- recipe(peaks ~ ., data = juiced_ups) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
  update_role(country_peaks, new_role = "ID")                   # when we do this, recipe is applied within each fold, and diff covariates are removed for each

# Checking Correlations - No Upsampling
rem <- c("country_India", "country_Iran", "peaks", "country_peaks")
rem_index <- which(colnames(juiced) %in% rem)
correlation_matrix <- cor(juiced[, -rem_index])
corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))
corr_vector <- as.vector(correlation_matrix[as.vector(correlation_matrix) != 1])
hist(unique(corr_vector))
mean(unique(corr_vector))

# Checking Correlations - Upsampling to Balance Peaks
rem <- c("country_India", "country_Iran", "peaks", "country_peaks")
rem_index <- which(colnames(juiced_ups) %in% rem)
correlation_matrix <- cor(juiced_ups[, -rem_index])
corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))
corr_vector <- as.vector(correlation_matrix[as.vector(correlation_matrix) != 1])
hist(unique(corr_vector))
mean(unique(corr_vector))

# Setting Up The Random Forest Framework, CV Folds and Performance Metrics to Evaluate
perf_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::sensitivity, yardstick::specificity)
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("classification") %>%
  set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  

## No upsampling
cv_splits <- vfold_cv(juiced, v = 6, strata = country_peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
  add_recipe(new_envt_recipe) %>%
  add_model(random_forest)
rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 20) # Creates a grid of hyperparameter values to try

## Upsampled
cv_splits_ups <- vfold_cv(juiced_ups, v = 6, strata = country_peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
rf_workflow_ups <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
  add_recipe(new_envt_recipe_ups) %>%
  add_model(random_forest)
rf_grid_ups <- grid_regular(mtry(range = c(2, dim(juiced_ups)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 20) # Creates a grid of hyperparameter values to try

#######################################################################################################
##                                                                                                   ##
##                      Running the Model In Parallel and Evaluating Performance                     ##
##                                                                                                   ##
#######################################################################################################

# Running the Model in Parallel 
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels)})

# Note - juicing is applied inside tune_grid i.e. to each training set associated with each cv fold
#        ==> means different predictors get removed by the recipe (because reducing sample size changes correlations)
#        ==> Solve by pre-juicing on full dataset, and then passing the juiced data in with a blank recipe.
# Checked by passing function to extract argument (as suggested here: https://tune.tidymodels.org/reference/tune_grid.html)
#                                                  add "extract = function (x) extract_recipe(x)))" to the control brackets

set.seed(345)
tune_res <- tune_grid(object = rf_workflow,
                      resamples = cv_splits,
                      grid = rf_grid,
                      metrics = perf_metrics,
                      control = control_resamples(save_pred = TRUE,
                                                  verbose = TRUE))
set.seed(345)
tune_res_ups <- tune_grid(object = rf_workflow_ups,
                          resamples = cv_splits_ups,
                          grid = rf_grid_ups,
                          metrics = perf_metrics,
                          control = control_resamples(save_pred = TRUE,
                                                      verbose = TRUE))

stopCluster(cl)

# Evaluating Performance
tune_res_perf <- tune_res %>%
  collect_metrics() 
tune_res_perf$type <- "no_ups"
tune_res_perf_ups <- tune_res_ups %>%
  collect_metrics() 
tune_res_perf_ups$type <- "ups"
tune_res_perf_overall <- rbind(tune_res_perf, tune_res_perf_ups)

tune_res_perf_overall %>%
  filter(.metric == "accuracy") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() + 
  facet_wrap(~type) +
  labs(y = "accuracy")

tune_res_perf_overall %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  facet_wrap(~type) +
  labs(y = "roc_auc")

tune_res_perf_overall %>%
  filter(.metric == "sens") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  facet_wrap(~type) +
  labs(y = "sensitivity")

tune_res_perf_overall %>%
  filter(.metric == "spec") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  facet_wrap(~type) +
  labs(y = "specificity")

#######################################################################################################
##                                                                                                   ##
##      Selecting the Best Fitting Set of Hyperparameters and Exploring Out-of-Bag Predictions       ##
##                                                                                                   ##
#######################################################################################################

# Exploring Quality of Fit and Selecting the Best Model 
show_best(tune_res, "roc_auc") # mtry = # predictors randomly sampled
show_best(tune_res_ups, "roc_auc") # mtry = # predictors randomly sampled
best_rmse <- select_best(tune_res, "roc_auc")
best_rmse_ups <- select_best(tune_res_ups, "roc_auc")
collect_predictions(tune_res) %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  group_by(id) %>%
  roc_curve(peaks, .pred_one) %>%
  autoplot()
collect_predictions(tune_res_ups) %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  group_by(id) %>%
  roc_curve(peaks, .pred_one) %>%
  autoplot()

# Checking performance on the different categories - VERY poor performance on 2 peaks for not upsampled data
predictions <- tune_res %>%
  collect_predictions() %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  dplyr::arrange(., .row)
one_peak <- predictions$peaks == "one"
two_peak <- predictions$peaks == "two"
sum(predictions$.pred_class[one_peak] == predictions$peaks[one_peak])/length(predictions$peaks[one_peak])
sum(predictions$.pred_class[two_peak] == predictions$peaks[two_peak])/length(predictions$peaks[two_peak])

predictions_ups <- tune_res_ups %>%
  collect_predictions() %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  dplyr::arrange(., .row)
one_peak_ups <- predictions_ups$peaks == "one"
two_peak_ups <- predictions_ups$peaks == "two"
sum(predictions_ups$.pred_class[one_peak] == predictions_ups$peaks[one_peak])/length(predictions_ups$peaks[one_peak])
sum(predictions_ups$.pred_class[two_peak] == predictions_ups$peaks[two_peak])/length(predictions_ups$peaks[two_peak])

#######################################################################################################
##                                                                                                   ##
##          Finalising Workflow and Re-Fitting to Entire Training Set, Evaluating Test Set           ##
##                                                                                                   ##
#######################################################################################################
# Finalising Workflows With Respective Best Hyperparameter Values  
random_forest_final <- rf_workflow %>%
  finalize_workflow(best_rmse)
random_forest_final_ups <- rf_workflow_ups %>%
  finalize_workflow(best_rmse_ups)

# Fitting a Random Forest With These Tuned Hyperparameters to the Entire Training Dataset
#  Note - random_forest_final has same recipe implicitly in there, so will process this data
#         according to that recipe. Important if you've got something in there that changes 
#         number of datapoints so as step_rose for oversampling. 
final_random_forest_fit <- random_forest_final %>%
  fit(data = juiced) 
fit_OOS_preds <- final_random_forest_fit$fit$fit$fit$predictions
fit_OOS_preds_df <- data.frame(truth = factor(juiced$peaks), estimate = fit_OOS_preds[, 1])
roc_auc(data = fit_OOS_preds_df, truth = truth, estimate)
peak_assign <- ifelse(fit_OOS_preds[, 1] > 0.50, "one", "two")
sum(peak_assign == juiced$peaks)/length(juiced$peaks)
sum(peak_assign[juiced$peaks == "one"] == juiced$peaks[juiced$peaks == "one"])/length(juiced$peaks[juiced$peaks == "one"])
sum(peak_assign[juiced$peaks == "two"] == juiced$peaks[juiced$peaks == "two"])/length(juiced$peaks[juiced$peaks == "two"])

final_random_forest_fit_ups <- random_forest_final_ups %>%
  fit(data = juiced_ups) 
fit_OOS_preds_ups <- final_random_forest_fit_ups$fit$fit$fit$predictions
fit_OOS_preds_df_ups <- data.frame(truth = factor(juiced_ups$peaks), estimate = fit_OOS_preds_ups[, 1])
roc_auc(data = fit_OOS_preds_df_ups, truth = truth, estimate)
peak_assign <- ifelse(fit_OOS_preds_ups[, 1] > 0.50, "one", "two")
sum(peak_assign == juiced_ups$peaks)/length(juiced_ups$peaks)
sum(peak_assign[juiced_ups$peaks == "one"] == juiced_ups$peaks[juiced_ups$peaks == "one"])/length(juiced_ups$peaks[juiced_ups$peaks == "one"])
sum(peak_assign[juiced_ups$peaks == "two"] == juiced_ups$peaks[juiced_ups$peaks == "two"])/length(juiced_ups$peaks[juiced_ups$peaks == "two"])

# Calculating variable importance
unregister_dopar()
extract_fit_parsnip(final_random_forest_fit) %>%
  vip(geom = "col", num_features = dim(juiced)[2])
extract_fit_parsnip(final_random_forest_fit_ups) %>%
  vip(geom = "col", num_features = dim(juiced_ups)[2])

# Collecting the evaluation metrics and assessing performance 
explainer <- explain_tidymodels(
  model = final_random_forest_fit,
  data = dplyr::select(juiced, -peaks),
  y = as.numeric(juiced$peaks),
  verbose = FALSE)
rem <- which(colnames(explainer$dat) == "country_peaks")
pdp <- model_profile(explainer, variables = colnames(explainer$data)[-rem], N = NULL)
plot(pdp)

explainer_ups <- explain_tidymodels(
  model = final_random_forest_fit_ups,
  data = dplyr::select(juiced_ups, -peaks),
  y = as.numeric(juiced_ups$peaks),
  verbose = FALSE)
rem_ups <- which(colnames(explainer_ups$dat) == "country_peaks")
pdp_ups <- model_profile(explainer_ups, variables = colnames(explainer_ups$data)[-rem_ups], N = NULL)
plot(pdp_ups)

# Scrap
# cv_splits <- vfold_cv(data, v = 6, strata = peaks) # v sets number of splits
# rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
#   add_recipe(envt_recipe) %>%
#   add_model(random_forest)

# df <- data.frame(log_pop = data$population_per_1km, city = overall$cit)
# boxplot(log_pop ~ city, df)
# How to access data used in a particular fold
#table(new_cv_splits$splits[[1]]$data$country_India[new_cv_splits$splits[[1]]$in_id])

# Checking to see which and if variables used differ between the folds (relevant for commented out tune_res above)
# tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # same removals within fold and recipe step
# tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # same removals within fold and recipe step
# tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # different removals between fold and within recipe step 3
# tune_res$.extracts[[2]]$.extracts[[1]]$steps[[3]]$removals # different removals between fold and within recipe step 3
# tune_res$.extracts[[1]]$.extracts[[1]]$steps[[4]]$removals # different removals between fold and within recipe step 4
# tune_res$.extracts[[2]]$.extracts[[1]]$steps[[4]]$removals # different removals between fold and within recipe step 4

# Below only applies to accuracy, but shows an important phenomenon
# Note: This accuracy calculation is NOT the same as reported in show_best exactly. But that's
#       because this is calculating average of each of the individual predictions' errors, rather 
#       average of the average error for each fold
# sum(order_raw_pred_best_hyp$.pred_class == order_raw_pred_best_hyp$peaks)/length(order_raw_pred_best_hyp$peaks)
# This DOES give us the same results as in show_best, because we're calculating the average of the 
# per fold performance average.
# fold_mean <- order_raw_pred_best_hyp %>%
#   group_by(id) %>%
#   summarise(prop = sum(.pred_class == peaks)/n())
# mean(fold_mean$prop)
