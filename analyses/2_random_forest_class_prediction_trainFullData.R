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

raw_envt_variables <- read.csv(here("data", "environmental_covariates", "location_ecological_data.csv")) %>%
  rename(id = Time.Series.ID, country = Country, admin1 = Admin.1, admin2 = Admin.2) %>%
  group_by(id, country, admin1, admin2) %>%
  summarise(across(population_per_1km:mean_temperature_driest_quarter, ~ mean(.x, na.rm = TRUE))) 

# aggregating LC variables into their top level values
envt_variables <- raw_envt_variables %>%
  mutate(LC_10 = LC_10 + LC_11 + LC_12) %>%
  dplyr::select(-LC_11, -LC_12) %>%
  mutate(LC_120 = LC_120 + LC_121 + LC_122) %>%
  dplyr::select(-LC_121, -LC_122) %>%
  mutate(LC_60 = LC_60 + LC_61 + LC_62) %>%
  dplyr::select(-LC_61, -LC_62) %>%
  mutate(LC_70 = LC_70 + LC_71) %>%
  dplyr::select(-LC_71) %>%
  mutate(LC_150 = LC_152 + LC_153) %>%
  dplyr::select(-LC_152, -LC_153) %>%
  mutate(LC_200 = LC_200 + LC_201 + LC_202) %>%
  dplyr::select(-LC_201, -LC_202)

colnames(raw_envt_variables)[grep("LC*", colnames(raw_envt_variables))]
colnames(envt_variables)[grep("LC*", colnames(envt_variables))]

cluster_membership <- readRDS(here("data", "systematic_review_results", "cluster_membership.rds"))
cluster_membership <- cluster_membership[, c("id", "cluster")]
overall <- ts_metadata %>%
  left_join(envt_variables, by = c("id", "country", "admin1", "admin2"))
overall <- overall %>%
  left_join(cluster_membership, by = "id")

# Adding monthly catch in
counts <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds")) %>%
  dplyr::select(id, Jan:Dec) %>%
  pivot_longer(cols = Jan:Dec, names_to = "month", values_to = "catch") %>%
  group_by(id) %>%
  summarise(n_months = sum(!is.na(catch)),
            monthly_catch = log(sum(catch, na.rm = TRUE)/n_months))
overall <- overall %>%
  left_join(counts, by = "id") %>%
  dplyr::select(-n_months)


#######################################################################################################
##                                                                                                   ##
##                      Setting Up All The Components Required to Run the Model                      ##
##             I.e. Test/Train Split, CV Folds, Random Forest Engine, Performance Metrics            ##
##                                                                                                   ##
#######################################################################################################
# Subsetting Outcome and Variables for Analysis + Log_10'ing Population (Easier for Visualisation Later On)
set.seed(234)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(cluster, country, rainfall_seas_3, monthly_catch, 
                population_per_1km:mean_temperature_driest_quarter,
                -LC_190, -elevation) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country_peaks = paste0(cluster, "_", country))
data$cluster <- ifelse(data$cluster == 1, "one", "two")
data$cluster <- as.factor(data$cluster)
data$population_per_1km <- log(data$population_per_1km + 1)

# Storage Tibble for Results
iterations <- tibble(seed = 1, model = list(1), iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                     cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                     test_predictions = list(1), test_roc_curve <- list(1),
                     test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                     importance = list(1), recipe = list(1))
iterations_ups <- tibble(seed = 1, model = list(1), iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                         cv_roc_auc = 1, cv_accuracy = 1, cv_one_peak_accuracy = 1, cv_two_peak_accuracy = 1,
                         test_predictions = list(1), test_roc_curve <- list(1),
                         test_roc_auc = 1, test_accuracy = 1, test_one_peak_accuracy = 1, test_two_peak_accuracy = 1,
                         importance = list(1), recipe = list(1))

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
    #step_center(all_numeric()) %>% # RFs don't need scaling and normalising, and you potentially lose info when doing it
    #step_scale(all_numeric()) %>% # RFs don't need scaling and normalising, and you potentially lose info when doing it
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "rainfall_seas_3", "monthly_catch")), threshold = 0.50) %>%
    step_dummy(country)
  envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
  juiced <- juice(envt_prepped)
  new_envt_recipe <- recipe(cluster ~ ., data = juiced) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country_peaks, new_role = "ID")           # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Prepping Data - Selecting Variables, Creating Recipe - Upsampling to Balance Amount of Data for 1 and 2 Peaks
  set.seed(915) ## does this need to be removed????
  envt_recipe_ups <- recipe(cluster  ~ ., data = data) %>%
    update_role(country_peaks, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    #step_center(all_numeric()) %>% # RFs don't need scaling and normalising, and you potentially lose info when doing it
    #step_scale(all_numeric()) %>% # RFs don't need scaling and normalising, and you potentially lose info when doing it
    step_nzv(all_numeric(), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "rainfall_seas_3", "monthly_catch")), threshold = 0.50) %>%
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
  iterations[i, "recipe"] <- list(list(envt_prepped))
  
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
  iterations_ups[i, "recipe"] <- list(list(envt_prepped_ups))
  
  print(paste0("Iteration ", i, " - Seed is ", seed))
  
}

saveRDS(iterations, file = paste0(here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData_ProperRainSeas.rds")))
saveRDS(iterations_ups, file = paste0(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData_ProperRainSeas.rds")))

# Loading in the fitted random forests and processing accuracy results 
no_ups <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData_ProperRainSeas.rds"))
ups <- readRDS(file = here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData_ProperRainSeas.rds"))
for (i in 1:length(no_ups$test_roc_curve)) {
  temp <- no_ups$test_roc_curve[[i]]
  temp_ups <- ups$test_roc_curve[[i]]
  if (i == 1) {
    df <- temp
    df$iteration <- 1
    df_ups <- temp_ups
    df_ups$iteration <- 1
  } else {
    temp$iteration <- i
    df <- rbind(df, temp)
    temp_ups$iteration <- i
    df_ups <- rbind(df_ups, temp_ups)
  }
}

# Plotting AUC and RF Results for Upsampled 
AUC_upsample_plot <- ggplot(df_ups, aes(x = 1-specificity, y = sensitivity, id = factor(iteration))) +
  geom_path(alpha = 0.5) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = "black", lty = 2) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("label", x = 0.67, y = 0.07,
           label = paste0("Mean AUC = ", round(mean(ups$test_roc_auc), 2)),
           label.padding = unit(0.35, "lines"), label.r = unit(0, "lines"),
           label.size = unit(0.35, "lines"), size = 4)

importance_upsample <- bind_rows(ups$importance)
importance_upsample <- importance_upsample %>%
  group_by(Variable) %>%
  summarise(mean_Importance = mean(Importance),
            stdev_Importance = sd(Importance),
            stder_Importance = sd(Importance)/sqrt(n()))
importance_upsample$lower <- pmax(rep(0, length(importance_upsample$mean_Importance)), 
                                  importance_upsample$mean_Importance - 1.96 * importance_upsample$stdev_Importance)
var_names_ups <- importance_upsample$Variable[order(importance_upsample$mean_Importance)]
new_names_ups <- var_names_ups
importance_upsample_plot <- ggplot(importance_upsample, aes(x = reorder(Variable, mean_Importance), y = mean_Importance, 
                                                            fill = mean_Importance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = pmax(0, mean_Importance - 1.96 * stdev_Importance),
                    ymax = mean_Importance + 1.96 * stdev_Importance),
                width = 0.5) +
  scale_x_discrete(labels = new_names_ups) +
  scale_fill_continuous(low = "grey", high = "#E14545") +
  xlab("") + ylab("Variable Importance") +
  lims(y = c(0, 0.086)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

rf_plot <- cowplot::plot_grid(AUC_upsample_plot, importance_upsample_plot, nrow = 1, ncol = 2, rel_widths = c(1, 2), align = "h", axis = "b")

# Supplementary Figure - Upsample PDP and Covariate Profiling Plots
ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_Upsampling_FullData_ProperRainSeas.rds"))
ups_vip_results <- tibble(id = 1, var = "bloop", x = 1, y = 1)
for (i in 1:dim(ups)[1]) {
  final_random_forest_fit_ups <- ups$model[[i]]
  explainer_ups <- explain_tidymodels(
    model = final_random_forest_fit_ups,
    data = dplyr::select(ups$juiced[[i]], -cluster),
    y = as.numeric(ups$juiced[[i]]$cluster),
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
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Covariate Value", y = "Average Prediction")
ggsave(filename = here("figures/Supp_Figure_UpsampleProperRainSeas_Covariate_Profiling.pdf"), plot = ups_profile_plots, width = 8, height = 8)

# Plotting No Upsampling AUC
no_ups_AUC <- ggplot(df, aes(x = 1-specificity, y = sensitivity, id = factor(iteration))) +
  geom_path(alpha = 0.5) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), col = "black", lty = 2) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  theme_bw() +
  theme(legend.position = "none") +
  annotate("label", x = 0.67, y = 0.07,
           label = paste0("Mean AUC = ", round(mean(no_ups$test_roc_auc), 2)),
           label.padding = unit(0.35, "lines"), label.r = unit(0, "lines"),
           label.size = unit(0.35, "lines"), size = 4)
no_ups_df <- bind_rows(no_ups$importance)
no_ups_df$data <- "full_data"
no_ups_imp <- no_ups_df %>%
  group_by(Variable) %>%
  summarise(mean_Importance = mean(Importance),
            stdev_Importance = sd(Importance),
            stder_Importance = sd(Importance)/sqrt(n()))
no_ups_imp$lower <- pmax(rep(0, length(no_ups_imp$mean_Importance)), no_ups_imp$mean_Importance - 1.96 * no_ups_imp$stdev_Importance)
no_ups_var_names <- no_ups_imp$Variable[order(no_ups_imp$mean_Importance)]
no_ups_new_names <- no_ups_var_names
importance_noUps_plot <- ggplot(no_ups_imp, aes(x = reorder(Variable, mean_Importance), y = mean_Importance, 
                                         fill = mean_Importance)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = pmax(0, mean_Importance - 1.96 * stdev_Importance),
                    ymax = mean_Importance + 1.96 * stdev_Importance),
                width = 0.5) +
  scale_x_discrete(labels = no_ups_new_names) +
  scale_fill_continuous(low = "grey", high = "#E14545") +
  xlab("") + ylab("Variable Importance") +
  lims(y = c(0, 0.065)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

noUps_Supp_Plot <- cowplot::plot_grid(no_ups_AUC, importance_noUps_plot, nrow = 1, ncol = 2, rel_widths = c(1, 2), align = "h", axis = "b")
noUps_Supp_Plot
ggsave(filename = here("figures/Supp_Figure_AUC_VIP_NoUps.pdf"), plot = noUps_Supp_Plot, width = 12, height = 5)


no_ups <- readRDS(here("outputs", "random_forest_outputs", "repeated_rf_noUpsampling_FullData_ProperRainSeas.rds"))
no_ups_vip_results <- tibble(id = 1, var = "bloop", x = 1, y = 1)
for (i in 1:dim(no_ups)[1]) {
  final_random_forest_fit_no_ups <- no_ups$model[[i]]
  explainer_no_ups <- explain_tidymodels(
    model = final_random_forest_fit_no_ups,
    data = dplyr::select(no_ups$juiced[[i]], -cluster),
    y = as.numeric(no_ups$juiced[[i]]$cluster),
    verbose = FALSE)
  rem_no_ups <- which(colnames(explainer_no_ups$dat) == "country_peaks")
  pdp_no_ups <- model_profile(explainer_no_ups, variables = colnames(explainer_no_ups$data)[-rem_no_ups], N = NULL)
  prof_no_ups <- pdp_no_ups$agr_profiles
  df_no_ups <- data.frame(id = i, var = prof_no_ups$`_vname_`, x = prof_no_ups$`_x_`, y = prof_no_ups$`_yhat_`)
  no_ups_vip_results <- rbind(no_ups_vip_results, df_no_ups)
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
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Covariate Value", y = "Average Prediction")
