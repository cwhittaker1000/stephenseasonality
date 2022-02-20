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
cluster_membership <- readRDS(here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_membership[, c("id", "cluster")]
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))
overall <- overall %>%
  left_join(cluster_membership, by = "id")

#######################################################################################################
##                                                                                                   ##
##                      Setting Up All The Components Required to Run the Model                      ##
##             I.e. Test/Train Split, CV Folds, Random Forest Engine, Performance Metrics            ##
##                                                                                                   ##
#######################################################################################################
# Subsetting Outcome and Variables for Analysis + Log_10'ing Population (Easier for Visualisation Later On)
set.seed(234)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(cluster, country, population_per_1km:mean_temperature_driest_quarter, -LC_190, -LC_210, -mean_temperature_driest_quarter) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country_peaks = paste0(cluster, "_", country))
data$cluster <- ifelse(data$cluster == 1, "one", "two")
data$cluster <- as.factor(data$cluster)
data$population_per_1km <- log(data$population_per_1km)

# Storage Tibble for Results
iterations <- tibble(seed = 1, model = list(1), iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                     cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                     test_predictions = list(1), test_roc_curve <- list(1),
                     test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                     importance = list(1))
iterations_ups <- tibble(seed = 1, model = list(1), iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                         cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                         test_predictions = list(1), test_roc_curve <- list(1),
                         test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                         importance = list(1))

# Running Multiple Random Forest Models Varying Seed Each Time
seeds <- c(234, 284, 102391, 19, 2948457, 294894, 189, 38484902, 284651, 83829, 72645,
           204758, 90, 7564, 5647, 8201, 3782, 284621, 98823, 762, 977503, 75759, 34428, 9, 3,
           10, 19, 19029, 8726, 23716, 17278, 92883, 827, 7162, 162, 1282, 8172, 128, 91, 981,
           456, 224, 8743, 362, 81, 9223, 753, 357, 99, 101)
number_iterations <- 25
for (i in 1:number_iterations) {
  
  # Set Seed
  set.seed(seeds[i])
  seed <- seeds[i]
  
  # Prepping Data - Selecting Variables, Creating Recipe - No Upsampling
  envt_recipe <- recipe(cluster  ~ ., data = data) %>%
    update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    step_center(all_numeric()) %>% 
    step_scale(all_numeric()) %>%
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
    step_dummy(country)
  envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
  juiced <- juice(envt_prepped)
  new_envt_recipe <- recipe(cluster ~ ., data = juiced) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country_peaks, new_role = "ID")           # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Prepping Data - Selecting Variables, Creating Recipe - Upsampling to Balance Amount of Data for 1 and 2 Peaks
  set.seed(915) ## does this need to be removed????
  envt_recipe_ups <- recipe(cluster  ~ ., data = data) %>%
    update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    step_center(all_numeric()) %>% 
    step_scale(all_numeric()) %>%
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv")), threshold = 0.50) %>%
    step_dummy(country) %>%
    step_smote(cluster)
  envt_prepped_ups <- prep(envt_recipe_ups, training = data, verbose = TRUE)
  juiced_ups <- juice(envt_prepped_ups)
  new_envt_recipe_ups <- recipe(cluster ~ ., data = juiced_ups) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country_peaks, new_role = "ID")                   # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Setting Up The Random Forest Framework
  perf_metrics <- metric_set(yardstick::roc_auc, yardstick::accuracy, yardstick::sensitivity, yardstick::specificity)
  random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  
  random_forest$eng_args$seed <- seed
  
  ## No Upsampling
  cv_splits <- vfold_cv(juiced, v = 6, strata = country_peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
  rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
    add_recipe(new_envt_recipe) %>%
    add_model(random_forest)
  rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 15) # Creates a grid of hyperparameter values to try
  
  ## Upsampled
  cv_splits_ups <- vfold_cv(juiced_ups, v = 6, strata = country_peaks) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
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
  
  # Extracting CV Predictions and Evaluating Accuracy
  
  ## No Upsampling
  best <- select_best(tune_res, "roc_auc")
  cv_predictions <- tune_res %>%
    collect_predictions() %>%
    dplyr::filter(mtry == best$mtry, min_n == best$min_n) %>%
    dplyr::arrange(., .row)
  cv_auc <- show_best(tune_res, "roc_auc")[1, "mean"]
  cv_accuracy <- sum(cv_predictions$.pred_class == cv_predictions$cluster)/length(cv_predictions$cluster)
  one_peak <- cv_predictions$cluster == "one"
  two_peak <- cv_predictions$cluster == "two"
  cv_one_peak_accuracy <- sum(cv_predictions$.pred_class[one_peak] == cv_predictions$cluster[one_peak])/length(cv_predictions$cluster[one_peak])
  cv_two_peak_accuracy <- sum(cv_predictions$.pred_class[two_peak] == cv_predictions$cluster[two_peak])/length(cv_predictions$cluster[two_peak])
  
  ## Upsampling
  best_ups <- select_best(tune_res_ups, "roc_auc")
  cv_predictions_ups <- tune_res_ups %>%
    collect_predictions() %>%
    dplyr::filter(mtry == best_ups$mtry, min_n == best_ups$min_n) %>%
    dplyr::arrange(., .row)
  cv_auc_ups <- show_best(tune_res_ups, "roc_auc")[1, "mean"]
  cv_accuracy_ups <- sum(cv_predictions_ups$.pred_class == cv_predictions_ups$cluster)/length(cv_predictions_ups$cluster)
  one_peak_ups <- cv_predictions_ups$cluster == "one"
  two_peak_ups <- cv_predictions_ups$cluster == "two"
  cv_one_peak_accuracy_ups <- sum(cv_predictions_ups$.pred_class[one_peak] == cv_predictions_ups$cluster[one_peak])/length(cv_predictions_ups$cluster[one_peak])
  cv_two_peak_accuracy_ups <- sum(cv_predictions_ups$.pred_class[two_peak] == cv_predictions_ups$cluster[two_peak])/length(cv_predictions_ups$cluster[two_peak])
  
  # Selecting the Best Model and Final Fitting
  unregister_dopar()
  
  ## No Upsampling
  random_forest_final <- rf_workflow %>%
    finalize_workflow(best)
  final_random_forest_fit <- random_forest_final %>%
    fit(data = juiced)
  fit_OOS_preds <- final_random_forest_fit$fit$fit$fit$predictions
  fit_OOS_preds_df <- data.frame(truth = factor(juiced$cluster), estimate = fit_OOS_preds[, 1])
  test_roc_auc <- roc_auc(data = fit_OOS_preds_df, truth = truth, estimate)
  test_preds <- fit_OOS_preds[, 1]
  test_roc_curve <- roc_curve(data = fit_OOS_preds_df, truth = "truth", estimate)
  peak_assign <- ifelse(fit_OOS_preds[, 1] > 0.50, "one", "two")
  test_accuracy <- sum(peak_assign == juiced$cluster)/length(juiced$cluster)
  test_one_peak_accuracy <- sum(peak_assign[juiced$cluster == "one"] == juiced$cluster[juiced$cluster == "one"])/length(juiced$cluster[juiced$cluster == "one"])
  test_two_peak_accuracy <- sum(peak_assign[juiced$cluster == "two"] == juiced$cluster[juiced$cluster == "two"])/length(juiced$cluster[juiced$cluster == "two"])
  
  ## Upsampling
  random_forest_final_ups <- rf_workflow_ups %>%
    finalize_workflow(best_ups)
  final_random_forest_fit_ups <- random_forest_final_ups %>%
    fit(data = juiced_ups)
  fit_OOS_preds_ups <- final_random_forest_fit_ups$fit$fit$fit$predictions
  fit_OOS_preds_df_ups <- data.frame(truth = factor(juiced_ups$cluster), estimate = fit_OOS_preds_ups[, 1])
  test_roc_auc_ups <- roc_auc(data = fit_OOS_preds_df_ups, truth = truth, estimate)
  test_preds_ups <- fit_OOS_preds_ups[, 1]
  test_roc_curve_ups <- roc_curve(data = fit_OOS_preds_df_ups, truth = "truth", estimate)
  peak_assign_ups <- ifelse(fit_OOS_preds_ups[, 1] > 0.50, "one", "two")
  test_accuracy_ups <- sum(peak_assign_ups == juiced_ups$cluster)/length(juiced_ups$cluster)
  test_one_peak_accuracy_ups <- sum(peak_assign_ups[juiced_ups$cluster == "one"] == juiced_ups$cluster[juiced_ups$cluster == "one"])/length(juiced_ups$cluster[juiced_ups$cluster == "one"])
  test_two_peak_accuracy_ups <- sum(peak_assign_ups[juiced_ups$cluster == "two"] == juiced_ups$cluster[juiced_ups$cluster == "two"])/length(juiced_ups$cluster[juiced_ups$cluster == "two"])
  
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
  iterations[i, "seed"] <- seed
  iterations[i, "model"] <- list(list(final_random_forest_fit))
  iterations[i, "iteration"] <- i
  iterations[i, "juiced"] <- list(list(juiced))
  iterations[i, "best_mtry"] <- best$mtry
  iterations[i, "best_min_n"] <- best$min_n
  iterations[i, "cv_roc_auc"] <- cv_auc
  iterations[i, "cv_accuracy"] <- cv_accuracy
  iterations[i, "cv_one_peak_accuracy"] <- cv_one_peak_accuracy
  iterations[i, "cv_two_peak_accuracy"] <- cv_two_peak_accuracy
  iterations[i, "test_predictions"] <- list(list(test_preds))
  iterations[i, "test_roc_curve"] <- list(list(test_roc_curve))
  iterations[i, "test_roc_auc"] <- test_roc_auc$.estimate
  iterations[i, "test_accuracy"] <- test_accuracy
  iterations[i, "test_one_peak_accuracy"] <- test_one_peak_accuracy
  iterations[i, "test_two_peak_accuracy"] <- test_two_peak_accuracy
  iterations[i, "importance"] <- list(list(var_imp))
  
  ## Upsampling
  iterations_ups[i, "seed"] <- seed
  iterations_ups[i, "model"] <- list(list(final_random_forest_fit_ups))
  iterations_ups[i, "iteration"] <- i
  iterations_ups[i, "juiced"] <- list(list(juiced_ups))
  iterations_ups[i, "best_mtry"] <- best_ups$mtry
  iterations_ups[i, "best_min_n"] <- best_ups$min_n
  iterations_ups[i, "cv_roc_auc"] <- cv_auc_ups
  iterations_ups[i, "cv_accuracy"] <- cv_accuracy_ups
  iterations_ups[i, "cv_one_peak_accuracy"] <- cv_one_peak_accuracy_ups
  iterations_ups[i, "cv_two_peak_accuracy"] <- cv_two_peak_accuracy_ups
  iterations_ups[i, "test_predictions"] <- list(list(test_preds_ups))
  iterations_ups[i, "test_roc_curve"] <- list(list(test_roc_curve_ups))
  iterations_ups[i, "test_roc_auc"] <- test_roc_auc_ups$.estimate
  iterations_ups[i, "test_accuracy"] <- test_accuracy_ups
  iterations_ups[i, "test_one_peak_accuracy"] <- test_one_peak_accuracy_ups
  iterations_ups[i, "test_two_peak_accuracy"] <- test_two_peak_accuracy_ups
  iterations_ups[i, "importance"] <- list(list(var_imp_ups))
  
  print(paste0("Iteration ", i, " - Seed is ", seed))
  
}

saveRDS(iterations, file = paste0(here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData.rds")))
saveRDS(iterations_ups, file = paste0(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds")))



ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData.rds"))
ups_vip_results <- tibble(id = 1, var = "bloop", x = 1, y = 1)
for (i in 1:dim(ups)[1]) {
  final_random_forest_fit_ups <- ups$model[[i]]
  explainer_ups <- explain_tidymodels(
    model = final_random_forest_fit_ups,
    data = dplyr::select(juiced_ups, -peaks),
    y = as.numeric(juiced_ups$peaks),
    verbose = FALSE)
  rem_ups <- which(colnames(explainer_ups$dat) == "country_peaks")
  pdp_ups <- model_profile(explainer_ups, variables = colnames(explainer_ups$data)[-rem_ups], N = NULL)
  prof_ups <- pdp_ups$agr_profiles
  df_ups <- data.frame(id = i, var = prof_ups$`_vname_`, x = prof_ups$`_x_`, y = prof_ups$`_yhat_`)
  ups_vip_results <- rbind(ups_vip_results, df_ups)
  print(i)
}
ups_vip_results <- ups_vip_results[-1, ] 
ups_summary_vip_results <- ups_vip_results %>%
  group_by(var, x) %>%
  summarise(mean = mean(y),
            lower = min(y),
            upper = max(y))
ups_profile_plots <- ggplot(ups_summary_vip_results, aes(x = x, y = mean, col = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = var), alpha = 0.2, col = NA) +
  facet_wrap(~var, scale = "free_x") +
  theme(legend.position = "none") +
  labs(x = "Covariate Value", y = "Average Prediction")
ggsave(filename = here("figures/Supp_Figure_Upsample_Covariate_Profiling.pdf"), plot = ups_profile_plots, width = 8, height = 8)

no_ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData.rds"))
no_ups_vip_results <- tibble(id = 1, var = "bloop", x = 1, y = 1)
for (i in 1:dim(no_ups)[1]) {
  final_random_forest_fit <- no_ups$model[[i]]
  explainer <- explain_tidymodels(
    model = final_random_forest_fit,
    data = dplyr::select(juiced, -peaks),
    y = as.numeric(juiced$peaks),
    verbose = FALSE)
  rem <- which(colnames(explainer$dat) == "country_peaks")
  pdp <- model_profile(explainer, variables = colnames(explainer$data)[-rem], N = NULL)
  prof <- pdp$agr_profiles
  df <- data.frame(id = i, var = prof$`_vname_`, x = prof$`_x_`, y = prof$`_yhat_`)
  no_ups_vip_results <- rbind(no_ups_vip_results, df)
  print(i)
}
no_ups_vip_results <- no_ups_vip_results[-1, ] 
no_ups_summary_vip_results <- no_ups_vip_results %>%
  group_by(var, x) %>%
  summarise(mean = mean(y),
            lower = min(y),
            upper = max(y))
no_ups_profile_plots <- ggplot(no_ups_summary_vip_results, aes(x = x, y = mean, col = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = var), alpha = 0.2, col = NA) +
  facet_wrap(~var, scale = "free_x") +
  theme(legend.position = "none") +
  labs(x = "Covariate Value", y = "Average Prediction")
ggsave(filename = here("figures/Supp_Figure_NoUpsample_Covariate_Profiling.pdf"), plot = no_ups_profile_plots, width = 8, height = 8)

# ups_profile_plots
# no_ups_profile_plots
# 
# no_ups_summary_vip_results$data <- "noUps"
# ups_summary_vip_results$data <- "Ups"
# overall <- rbind(no_ups_summary_vip_results, ups_summary_vip_results)
# ggplot(overall, aes(x = x, y = mean, col = data)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = data), alpha = 0.2, col = NA) +
#   facet_wrap(~var, scale = "free_x") +
#   theme(legend.position = "none") +
#   labs(x = "Covariate Value", y = "Average Prediction")
# 
# x <- bind_rows(iterations$importance) %>%
#   group_by(Variable) %>%
#   summarise(mean = mean(Importance),
#             sd = sd(Importance),
#             se = sd(Importance)/sqrt(n()))
# x[rev(order(x$mean)), ]
# x <- bind_rows(iterations_ups$importance) %>%
#   group_by(Variable) %>%
#   summarise(mean = mean(Importance),
#             sd = sd(Importance),
#             se = sd(Importance)/sqrt(n()))
# x[rev(order(x$mean)), ]
# x <- bind_rows(iterations_ups$importance, iterations$importance) %>%
#   pivot_wider(names_from = sampling, values_from = Importance) %>%
#   group_by(Variable) %>%
#   summarise(mean_ups = mean(upsampling), 
#             mean_no = mean(no_upsampling))
# x <- bind_rows(iterations_ups$importance, iterations$importance) 
# ggplot(x, aes(x= Variable, fill = sampling, y = Importance)) +
#   geom_boxplot()
