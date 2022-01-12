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

# Storage Tibble for Results
iterations <- tibble(iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                     cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                     test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                     importance = list(1))
iterations_ups <- tibble(iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                         cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                         test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                         importance = list(1))

# Running Multiple Random Forest Models Varying Seed Each Time
seeds <- c(345, 234, 284, 102391, 19, 2948457, 294894, 189, 38484902, 284651, 83829, 72645,
           204758, 90, 7564, 5647, 8201, 3782, 284621, 98823, 762, 977503, 75759, 34428, 9, 3,
           10, 19, 19029, 8726, 23716, 17278, 92883, 827, 7162, 162, 1282, 8172, 128, 91, 981,
           456, 224, 8743, 362, 81, 9223, 753, 357, 99, 101)
number_iterations <- 20
for (i in 1:number_iterations) {
  
  # Set Seed
  set.seed(seeds[i])
  seed <- seeds[i]
  
  # Creating Test and Training Data
  set.seed(seed)
  rf_split <- initial_split(data, prop = 0.85, strata = country_peaks)
  rf_train <- training(rf_split)
  rf_test <- testing(rf_split)
  
  # Prepping Data - Selecting Variables, Creating Recipe - No Upsampling
  envt_recipe <- recipe(peaks  ~ ., data = rf_train) %>%
    update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    step_center(all_numeric()) %>% 
    step_scale(all_numeric()) %>%
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
    step_dummy(country)
  envt_prepped <- prep(envt_recipe, training = rf_train, verbose = TRUE)
  juiced <- juice(envt_prepped)
  new_envt_recipe <- recipe(peaks ~ ., data = juiced) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country_peaks, new_role = "ID")           # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Prepping Data - Selecting Variables, Creating Recipe - Upsampling to Balance Amount of Data for 1 and 2 Peaks
  set.seed(seed)
  envt_recipe_ups <- recipe(peaks  ~ ., data = rf_train) %>%
    update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    step_center(all_numeric()) %>% 
    step_scale(all_numeric()) %>%
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
    step_dummy(country) %>%
    step_smote(peaks)
  envt_prepped_ups <- prep(envt_recipe_ups, training = rf_train, verbose = TRUE)
  juiced_ups <- juice(envt_prepped_ups)
  new_envt_recipe_ups <- recipe(peaks ~ ., data = juiced_ups) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country_peaks, new_role = "ID")                   # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Setting Up The Random Forest Framework
  perf_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::sensitivity, yardstick::specificity)
  random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  
  random_forest$eng_args$seed <- seed
  
  ## No Upsampling
  cv_splits <- vfold_cv(juiced, v = 6, strata = peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
  rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
    add_recipe(new_envt_recipe) %>%
    add_model(random_forest)
  rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 15) # Creates a grid of hyperparameter values to try
  
  ## Upsampled
  cv_splits_ups <- vfold_cv(juiced_ups, v = 6, strata = peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
  rf_workflow_ups <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
    add_recipe(new_envt_recipe_ups) %>%
    add_model(random_forest)
  rf_grid_ups <- grid_regular(mtry(range = c(2, dim(juiced_ups)[2] - round(dim(juiced_ups)[2]/4))), min_n(range = c(2, dim(juiced_ups)[1] - (round(dim(juiced_ups)[1]/3)))), levels = 15) # Creates a grid of hyperparameter values to try
  
  # Running the Model in Parallel 
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {library(tidymodels)})
  
  ## No Upsampling
  set.seed(seeds[i])
  tune_res <- tune_grid(object = rf_workflow,
                        resamples = cv_splits,
                        grid = rf_grid,
                        metrics = perf_metrics,
                        control = control_resamples(save_pred = TRUE, verbose = TRUE))
  
  ## Upsampling
  set.seed(seeds[i])
  tune_res_ups <- tune_grid(object = rf_workflow_ups,
                            resamples = cv_splits_ups,
                            grid = rf_grid_ups,
                            metrics = perf_metrics,
                            control = control_resamples(save_pred = TRUE, verbose = TRUE))
  
  stopCluster(cl)

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
  
  # Extracting CV Predictions and Evaluating Accuracy
  
  ## No Upsampling
  best <- select_best(tune_res, "roc_auc")
  cv_predictions <- tune_res %>%
    collect_predictions() %>%
    dplyr::filter(mtry == best$mtry, min_n == best$min_n) %>%
    dplyr::arrange(., .row)
  cv_auc <- show_best(tune_res, "roc_auc")[1, "mean"]
  cv_accuracy <- sum(cv_predictions$.pred_class == cv_predictions$peaks)/length(cv_predictions$peaks)
  one_peak <- cv_predictions$peaks == "one"
  two_peak <- cv_predictions$peaks == "two"
  cv_one_peak_accuracy <- sum(cv_predictions$.pred_class[one_peak] == cv_predictions$peaks[one_peak])/length(cv_predictions$peaks[one_peak])
  cv_two_peak_accuracy <- sum(cv_predictions$.pred_class[two_peak] == cv_predictions$peaks[two_peak])/length(cv_predictions$peaks[two_peak])
  
  ## Upsampling
  best_ups <- select_best(tune_res_ups, "roc_auc")
  cv_predictions_ups <- tune_res_ups %>%
    collect_predictions() %>%
    dplyr::filter(mtry == best_ups$mtry, min_n == best_ups$min_n) %>%
    dplyr::arrange(., .row)
  cv_auc_ups <- show_best(tune_res_ups, "roc_auc")[1, "mean"]
  cv_accuracy_ups <- sum(cv_predictions_ups$.pred_class == cv_predictions_ups$peaks)/length(cv_predictions_ups$peaks)
  one_peak_ups <- cv_predictions_ups$peaks == "one"
  two_peak_ups <- cv_predictions_ups$peaks == "two"
  cv_one_peak_accuracy_ups <- sum(cv_predictions_ups$.pred_class[one_peak_ups] == cv_predictions_ups$peaks[one_peak_ups])/length(cv_predictions_ups$peaks[one_peak_ups])
  cv_two_peak_accuracy_ups <- sum(cv_predictions_ups$.pred_class[two_peak_ups] == cv_predictions_ups$peaks[two_peak_ups])/length(cv_predictions_ups$peaks[two_peak_ups])
  
  # Selecting the Best Model and Final Fitting
  unregister_dopar()
  
  ## No Upsampling
  random_forest_final <- rf_workflow %>%
    finalize_workflow(best)
  final_random_forest_fit <- random_forest_final %>%
    fit(data = juiced)
  baked_train <- bake(envt_prepped, rf_test)
  full_forest_pred <- predict(final_random_forest_fit, baked_train, "prob")
  full_forest_preds_df <- data.frame(truth = factor(baked_train$peaks), estimate = unname(full_forest_pred[, 1]))
  test_roc_auc <- roc_auc(data = full_forest_preds_df, truth = truth, estimate)
  peak_assign <- ifelse(full_forest_pred[, 1] > 0.50, "one", "two")
  test_accuracy <- unname(sum(peak_assign == rf_test$peaks)/length(rf_test$peaks))
  test_one_peak_accuracy <- sum(peak_assign[rf_test$peaks == "one"] == rf_test$peaks[rf_test$peaks == "one"])/length(rf_test$peaks[rf_test$peaks == "one"])
  test_two_peak_accuracy <- sum(peak_assign[rf_test$peaks == "two"] == rf_test$peaks[rf_test$peaks == "two"])/length(rf_test$peaks[rf_test$peaks == "two"])
  
  ## Upsampling
  random_forest_final_ups <- rf_workflow_ups %>%
    finalize_workflow(best_ups)
  final_random_forest_fit_ups <- random_forest_final_ups %>%
    fit(data = juiced_ups)
  baked_train_ups <- bake(envt_prepped_ups, rf_test)
  full_forest_pred_ups <- predict(final_random_forest_fit_ups, baked_train_ups, "prob")
  full_forest_preds_df_ups <- data.frame(truth = factor(baked_train_ups$peaks), estimate = unname(full_forest_pred_ups[, 1]))
  test_roc_auc_ups <- roc_auc(data = full_forest_preds_df_ups, truth = truth, estimate)
  peak_assign_ups <- ifelse(full_forest_pred_ups[, 1] > 0.50, "one", "two")
  test_accuracy_ups <- unname(sum(peak_assign_ups == rf_test$peaks)/length(rf_test$peaks))
  test_one_peak_accuracy_ups <- sum(peak_assign_ups[rf_test$peaks == "one"] == rf_test$peaks[rf_test$peaks == "one"])/length(rf_test$peaks[rf_test$peaks == "one"])
  test_two_peak_accuracy_ups <- sum(peak_assign_ups[rf_test$peaks == "two"] == rf_test$peaks[rf_test$peaks == "two"])/length(rf_test$peaks[rf_test$peaks == "two"])
  
  # Calculating variable importance
  
  ## No Upsampling
  var_imp <- extract_fit_parsnip(final_random_forest_fit) %>%
    vip(num_features = dim(juiced)[2])
  var_imp <- var_imp$data
  var_imp$id <- i
  var_imp$sampling <- "no_upsampling"
  
  ## Upsampling
  var_imp_ups <- extract_fit_parsnip(final_random_forest_fit_ups) %>%
    vip(num_features = dim(juiced_ups)[2])
  var_imp_ups <- var_imp_ups$data
  var_imp_ups$id <- i
  var_imp_ups$sampling <- "upsampling"
  
  # Assigning Outputs to Tibble
  
  # No Upsampling
  iterations[i, "iteration"] <- i
  iterations[i, "juiced"] <- list(list(juiced))
  iterations[i, "best_mtry"] <- best$mtry
  iterations[i, "best_min_n"] <- best$min_n
  iterations[i, "cv_roc_auc"] <- cv_auc
  iterations[i, "cv_accuracy"] <- cv_accuracy
  iterations[i, "cv_one_peak_accuracy"] <- cv_one_peak_accuracy
  iterations[i, "cv_two_peak_accuracy"] <- cv_two_peak_accuracy
  iterations[i, "test_roc_auc"] <- test_roc_auc$.estimate
  iterations[i, "test_accuracy"] <- test_accuracy
  iterations[i, "test_one_peak_accuracy"] <- test_one_peak_accuracy
  iterations[i, "test_two_peak_accuracy"] <- test_two_peak_accuracy
  iterations[i, "importance"] <- list(list(var_imp))
  
  ## Upsampling
  iterations_ups[i, "iteration"] <- i
  iterations_ups[i, "juiced"] <- list(list(juiced_ups))
  iterations_ups[i, "best_mtry"] <- best_ups$mtry
  iterations_ups[i, "best_min_n"] <- best_ups$min_n
  iterations_ups[i, "cv_roc_auc"] <- cv_auc_ups
  iterations_ups[i, "cv_accuracy"] <- cv_accuracy_ups
  iterations_ups[i, "cv_one_peak_accuracy"] <- cv_one_peak_accuracy_ups
  iterations_ups[i, "cv_two_peak_accuracy"] <- cv_two_peak_accuracy_ups
  iterations_ups[i, "test_roc_auc"] <- test_roc_auc_ups$.estimate
  iterations_ups[i, "test_accuracy"] <- test_accuracy_ups
  iterations_ups[i, "test_one_peak_accuracy"] <- test_one_peak_accuracy_ups
  iterations_ups[i, "test_two_peak_accuracy"] <- test_two_peak_accuracy_ups
  iterations_ups[i, "importance"] <- list(list(var_imp_ups))
  
  print(paste0("Iteration ", i, " - Seed is ", seed))
  
}

x <- bind_rows(iterations$importance) %>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x[rev(order(x$mean)), ]

x <- bind_rows(iterations_ups$importance) %>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x[rev(order(x$mean)), ]


x <- bind_rows(iterations_ups$importance, iterations$importance) %>%
  pivot_wider(names_from = sampling, values_from = Importance) %>%
  group_by(Variable) %>%
  summarise(mean_ups = mean(upsampling), 
            mean_no = mean(no_upsampling))


x <- bind_rows(iterations_ups$importance, iterations$importance) 
ggplot(x, aes(x= Variable, fill = sampling, y = Importance)) +
  geom_boxplot()


