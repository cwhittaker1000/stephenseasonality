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
ts_metadata <- readRDS(here("data", "processed", "metadata_and_time_series_features.rds"))
envt_variables <- read.csv(here("data", "processed", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:worldclim_9, ~ mean(.x, na.rm = TRUE)))
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2")) 

#######################################################################################################
##                                                                                                   ##
##                      Setting Up All The Components Required to Run the Model                      ##
##             I.e. Test/Train Split, CV Folds, Random Forest Engine, Performance Metrics            ##
##                                                                                                   ##
#######################################################################################################
# Subsetting Outcome and Variables for Analysis
data <- overall %>% 
  dplyr::select(peaks, country, population_per_1km:worldclim_9, -LC_190) %>%
  mutate(country = as.factor(country))
data$peaks <- ifelse(data$peaks == 1, "one", "two")
data$peaks <- as.factor(data$peaks)

# Storage Tibble for Results
iterations <- tibble(iteration = 1, 
                     training_data = list(1),
                     test_data = list(1),
                     juiced = list(1),
                     best_mtry = 1,  
                     best_min_n = 1,
                     cv_roc_auc = 1,
                     cv_accuracy = 1,
                     cv_one_peak_accuracy = 1, 
                     cv_two_peak_accuracy = 1,
                     test_roc_auc = 1,
                     test_accuracy = 1,
                     test_one_peak_accuracy = 1,
                     test_two_peak_accuracy = 1,
                     importance = list(1))

# Running Multiple Random Forest Models Varying Seed Each Time
seeds <- c(234, 284, 102391, 19, 2948457, 294894, 189, 38484902, 284651, 83829, 72645,
           204758, 90, 7564, 5647, 8201, 3782, 284621, 98823, 762, 977503, 75759, 34428, 9, 3,
           10, 19, 19029, 8726, 23716, 17278, 92883, 827, 7162, 162, 1282, 8172, 128, 91, 981,
           456, 224, 8743, 362, 81, 9223, 753, 357, 99, 101)
for (i in 1:length(seeds)) {
  
  # Set Seed
  set.seed(seeds[i])
  seed <- seeds[i]
  
  # Creating Test and Training Data
  rf_split <- initial_split(data, prop = 0.85, strata = country)
  rf_train <- training(rf_split)
  rf_test <- testing(rf_split)
  
  # Prepping Data - Selecting Variables, Creating Recipe
  envt_recipe <- recipe(peaks  ~ ., data = rf_train) %>%
    update_role(country, new_role = "ID") %>% 
    step_log(population_per_1km) %>%
    step_center(all_predictors()) %>% 
    step_scale(all_predictors()) %>%
    step_nzv(all_predictors(), freq_cut = 80/20) %>% 
    step_corr(all_predictors(), threshold = 0.65) 
  envt_prepped <- prep(envt_recipe, training = rf_train, verbose = FALSE)
  juiced <- juice(envt_prepped)
  
  # Generating CV Folds and Selecting the Performance Metric We'll Evaluate Each Time  
  cv_splits <- vfold_cv(rf_train, v = 6, strata = peaks) 
  
  # Setting Up The Random Forest Framework
  random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "permutation", seed = 111) # dummy seed
  random_forest$eng_args$seed <- seed
  rf_workflow <- workflow() %>% 
    add_recipe(envt_recipe) %>%
    add_model(random_forest)
  rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2]-3)), 
                          min_n(range = c(2, dim(juiced)[1] - 10)), levels = 25)
  
  # Running the Model in Parallel 
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {library(tidymodels)})
  set.seed(seeds[i])
  tune_res <- tune_grid(object = rf_workflow,
                        resamples = cv_splits,
                        grid = rf_grid,
                        control = control_resamples(save_pred = TRUE))
  stopCluster(cl)
  
  # Extracting CV Predictions and Evaluating Accuracy
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
  
  # Selecting the Best Model and Final Fitting
  random_forest_final <- rf_workflow %>%
    finalize_workflow(best)
  unregister_dopar()
  final_res <- random_forest_final %>%
    last_fit(rf_split) 
  test_metrics <- final_res %>%
    collect_metrics()
  test_auc <- test_metrics$.estimate[test_metrics$.metric == "roc_auc"]
  test_accuracy <- test_metrics$.estimate[test_metrics$.metric == "accuracy"]
  test_one_peak <- final_res$.predictions[[1]]$peaks == "one"
  test_two_peak <- final_res$.predictions[[1]]$peaks == "two"
  test_one_peak_accuracy <- sum(final_res$.predictions[[1]]$peaks[test_one_peak] == 
                                  final_res$.predictions[[1]]$.pred_class[test_one_peak])/length(final_res$.predictions[[1]]$.pred_class[test_one_peak])
  test_two_peak_accuracy <- sum(final_res$.predictions[[1]]$peaks[test_two_peak] == 
                                  final_res$.predictions[[1]]$.pred_class[test_two_peak])/length(final_res$.predictions[[1]]$.pred_class[test_two_peak])
  if(sum(test_two_peak) == 0) {
    test_two_peak_accuracy <- NA
  }
  
  # Calculating variable importance
  var_imp <- extract_workflow(final_res) %>%
    extract_fit_parsnip() %>%
    vip(num_features = dim(juiced)[2])
  var_imp <- var_imp$data
  var_imp$id <- i
  
  # Assigning Outputs to Tibble
  iterations[i, "iteration"] <- i
  iterations[i, "training_data"] <- list(list(rf_train))
  iterations[i, "test_data"] <- list(list(rf_test))
  iterations[i, "juiced"] <- list(list(juiced))
  iterations[i, "best_mtry"] <- best$mtry
  iterations[i, "best_min_n"] <- best$min_n
  iterations[i, "cv_roc_auc"] <- cv_auc
  iterations[i, "cv_accuracy"] <- cv_accuracy
  iterations[i, "cv_one_peak_accuracy"] <- cv_one_peak_accuracy
  iterations[i, "cv_two_peak_accuracy"] <- cv_two_peak_accuracy
  iterations[i, "test_roc_auc"] <- test_auc
  iterations[i, "test_accuracy"] <- test_accuracy
  iterations[i, "test_one_peak_accuracy"] <- test_one_peak_accuracy
  iterations[i, "test_two_peak_accuracy"] <- test_two_peak_accuracy
  iterations[i, "importance"] <- list(list(var_imp))
  
  print(paste0("Iteration ", i, " - Seed is ", seed))
  
}


x <- bind_rows(iterations$importance) %>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x[rev(order(x$mean)), ]
