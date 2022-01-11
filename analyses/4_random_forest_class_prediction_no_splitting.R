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
# df <- data.frame(log_pop = data$population_per_1km, city = overall$cit)
# boxplot(log_pop ~ city, df)

# Prepping Data - Selecting Variables, Creating Recipe
set.seed(915)
envt_recipe <- recipe(peaks  ~ ., data = data) %>%
  update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
  step_center(all_numeric()) %>% 
  step_scale(all_numeric()) %>%
  step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
  step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
  step_dummy(country)# %>%
  #step_smote(peaks)
envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
juiced <- juice(envt_prepped)
colnames(juiced)
dim(juiced)
table(juiced$peaks)
new_envt_recipe <- recipe(peaks ~ ., data = juiced) %>%
  update_role(country_peaks, new_role = "ID")

# Checking Correlations
rem <- c("country_India", "country_Iran", "peaks", "country_peaks")
rem_index <- which(colnames(juiced) %in% rem)
correlation_matrix <- cor(juiced[, -rem_index])
corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))
corr_vector <- as.vector(correlation_matrix[as.vector(correlation_matrix) != 1])
hist(unique(corr_vector))
mean(unique(corr_vector))

# Setting Up The Random Forest Framework, CV Folds and Performance Metrics to Evaluate
perf_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::sensitivity, yardstick::specificity)
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("classification") %>%
  set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  

cv_splits <- vfold_cv(data, v = 6, strata = peaks) # v sets number of splits
rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
  add_recipe(envt_recipe) %>%
  add_model(random_forest)

new_cv_splits <- vfold_cv(juiced, v = 6, strata = country_peaks) # new cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
new_rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
  add_recipe(new_envt_recipe) %>%
  add_model(random_forest)

rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 15) # Creates a grid of hyperparameter values to try

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
#        ==> Unclear how much of a problem this is.
# Checked by passing function to extract argument (as suggested here: https://tune.tidymodels.org/reference/tune_grid.html)
set.seed(345)
# tune_res <- tune_grid(object = rf_workflow,
#                       resamples = cv_splits,
#                       grid = rf_grid,
#                       metrics = perf_metrics, 
#                       control = control_resamples(save_pred = TRUE,
#                                                   verbose = TRUE,
#                                                   extract = function (x) extract_recipe(x)))
tune_res <- tune_grid(object = new_rf_workflow,
                      resamples = new_cv_splits,
                      grid = rf_grid,
                      metrics = perf_metrics, 
                      control = control_resamples(save_pred = TRUE,
                                                  verbose = TRUE))
stopCluster(cl)

# How to access data used in a particular fold
table(new_cv_splits$splits[[1]]$data$country_India[new_cv_splits$splits[[1]]$in_id])

# Checking to see which and if variables used differ between the folds (relevant for commented out tune_res above)
tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # same removals within fold and recipe step
tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # same removals within fold and recipe step
tune_res$.extracts[[1]]$.extracts[[1]]$steps[[3]]$removals # different removals between fold and within recipe step 3
tune_res$.extracts[[2]]$.extracts[[1]]$steps[[3]]$removals # different removals between fold and within recipe step 3
tune_res$.extracts[[1]]$.extracts[[1]]$steps[[4]]$removals # different removals between fold and within recipe step 4
tune_res$.extracts[[2]]$.extracts[[1]]$steps[[4]]$removals # different removals between fold and within recipe step 4

# Evaluating Performance
tune_res %>%
  collect_metrics() %>%
  filter(.metric == "accuracy") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "accuracy")

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "roc_auc")

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "sens") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "sensitivity")

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "spec") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "specificity")

#######################################################################################################
##                                                                                                   ##
##      Selecting the Best Fitting Set of Hyperparameters and Exploring Out-of-Bag Predictions       ##
##                                                                                                   ##
#######################################################################################################

# Exploring Quality of Fit and Selecting the Best Model 
show_best(tune_res, "roc_auc") # mtry = # predictors randomly sampled
show_best(tune_res, "accuracy") # min_n = minmum number of data points required for node to be split further
best_rmse <- select_best(tune_res, "roc_auc")

collect_predictions(tune_res) %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  group_by(id) %>%
  roc_curve(peaks, .pred_one) %>%
  autoplot()

# Extracting and plotting raw predictions
raw_predictions <- tune_res %>%
  collect_predictions()
raw_pred_best_hyp <- raw_predictions %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n)
order_raw_pred_best_hyp <- raw_pred_best_hyp %>%
  dplyr::arrange(., .row)

# Below only applies to accuracy, but shows an important phenomenon
# Note: This accuracy calculation is NOT the same as reported in show_best exactly. But that's
#       because this is calculating average of each of the individual predictions' errors, rather 
#       average of the average error for each fold
sum(order_raw_pred_best_hyp$.pred_class == order_raw_pred_best_hyp$peaks)/length(order_raw_pred_best_hyp$peaks)
# This DOES give us the same results as in show_best, because we're calculating the average of the 
# per fold performance average.
fold_mean <- order_raw_pred_best_hyp %>%
  group_by(id) %>%
  summarise(prop = sum(.pred_class == peaks)/n())
mean(fold_mean$prop)

# Checking performance on the different categories - VERY poor performance on 2 peaks, good on 1 peak
one_peak <- order_raw_pred_best_hyp$peaks == "one"
two_peak <- order_raw_pred_best_hyp$peaks == "two"
sum(order_raw_pred_best_hyp$.pred_class[one_peak] == order_raw_pred_best_hyp$peaks[one_peak])/length(order_raw_pred_best_hyp$peaks[one_peak])
sum(order_raw_pred_best_hyp$.pred_class[two_peak] == order_raw_pred_best_hyp$peaks[two_peak])/length(order_raw_pred_best_hyp$peaks[two_peak])

#######################################################################################################
##                                                                                                   ##
##          Finalising Workflow and Re-Fitting to Entire Training Set, Evaluating Test Set           ##
##                                                                                                   ##
#######################################################################################################
# Finalising Workflow With Best Hyperparameter Values  
# random_forest_final <- rf_workflow %>%
#   finalize_workflow(best_rmse)
random_forest_final <- new_rf_workflow %>%
  finalize_workflow(best_rmse)

# Fitting a Random Forest With These Tuned Hyperparameters to the Entire Training Dataset
#  Note - random_forest_final has same recipe implicitly in there, so will process this data
#         according to that recipe. Important if you've got something in there that changes 
#         number of datapoints so as step_rose for oversampling. 
# final_random_forest_fit <- random_forest_final %>%
#   fit(data = data) 
final_random_forest_fit <- random_forest_final %>%
  fit(data = juiced) 
#x <- data.frame(truth = factor(data$peaks), estimate = final_random_forest_fit$fit$fit$fit$predictions[, 1])
x <- data.frame(truth = factor(juiced$peaks), estimate = final_random_forest_fit$fit$fit$fit$predictions[, 1])
roc_auc(data = x, truth = truth, estimate)
final_train_predictions <- final_random_forest_fit$fit$fit$fit$predictions
final_train_predictions <- ifelse(final_train_predictions[, 1] > 0.50, "one", "two")
#sum(final_train_predictions == data$peaks)/length(data$peaks)
sum(final_train_predictions == juiced$peaks)/length(juiced$peaks)

# Calculating variable importance
unregister_dopar()
extract_fit_parsnip(final_random_forest_fit) %>%
  vip(geom = "point")

# Collecting the evaluation metrics and assessing performance 
explainer <- explain_tidymodels(
  model = final_random_forest_fit,
  data = dplyr::select(data, -peaks),
  y = as.numeric(data$peaks),
  verbose = FALSE
)

explainer <- explain_tidymodels(
  model = final_random_forest_fit,
  data = dplyr::select(juiced, -peaks),
  y = as.numeric(juiced$peaks),
  verbose = FALSE
)

pdp_time <- model_profile(
  explainer,
  variables = colnames(explainer$data)[-15], 
  N = NULL)
plot(pdp_time)

pdp_time <- model_profile(
  explainer,
  variables = "precipitation_seasonality_cv", 
  N = NULL)
plot(pdp_time)

pdp_time <- model_profile(
  explainer,
  variables = colnames, 
  N = NULL)
plot(pdp_time)