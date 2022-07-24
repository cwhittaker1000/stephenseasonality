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

# colnames(raw_envt_variables)[grep("LC*", colnames(raw_envt_variables))]
# colnames(envt_variables)[grep("LC*", colnames(envt_variables))]

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
  dplyr::select(per_ind_3_months, country, rainfall_seas_3, monthly_catch, 
                population_per_1km:mean_temperature_driest_quarter, 
                -LC_190, -elevation) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) 
data$population_per_1km <- log(data$population_per_1km)
data$country2 <- data$country # new var for country so we can both dummy/hot-one "country" and retain "country2" for stratification purposes in CV

# Storage Tibble for Results
iterations <- tibble(seed = 1, model = list(1), iteration = 1, juiced = list(1), best_mtry = 1, best_min_n = 1,
                     cv_rsme = 1, test_predictions = list(1), test_rsme = 1, importance = list(1), recipe = list(1))

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
  envt_recipe <- recipe(per_ind_3_months  ~ ., data = data) %>%
    update_role(country2, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
    #step_center(all_numeric(), -per_ind_3_months) %>% 
    #step_scale(all_numeric(), -per_ind_3_months) %>%
    step_nzv(all_numeric(), -per_ind_3_months, freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
    step_corr(all_numeric(), -per_ind_3_months, -contains(c("population_per_1km", "temperature_seasonality", "rainfall_seas_3", "monthly_catch")), threshold = 0.50) %>%
    step_dummy(country)
  envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
  juiced <- juice(envt_prepped)
  new_envt_recipe <- recipe(per_ind_3_months ~ ., data = juiced) %>% # have to do this way to ensure recipe isn't applied *within* tune-grid
    update_role(country2, new_role = "ID")                           # when we do this, recipe is applied within each fold, and diff covariates are removed for each
  
  # Setting Up The Random Forest Framework
  perf_metrics <- metric_set(yardstick::rmse)
  random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
    set_mode("regression") %>%
    set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  
  random_forest$eng_args$seed <- seed
  
  # No Upsampling
  cv_splits <- vfold_cv(juiced, v = 6, strata = country2) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
  rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
    add_recipe(new_envt_recipe) %>%
    add_model(random_forest)
  rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 15) # Creates a grid of hyperparameter values to try

  # Running the Model in Parallel 
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {library(tidymodels)})
  
  set.seed(seeds[i])
  tune_res <- tune_grid(object = rf_workflow,
                        resamples = cv_splits,
                        grid = rf_grid,
                        metrics = perf_metrics,
                        control = control_resamples(save_pred = TRUE, verbose = TRUE))
  
  stopCluster(cl)
  
  # Extracting CV Predictions and Evaluating Accuracy
  best <- select_best(tune_res, "rmse")
  cv_predictions <- tune_res %>%
    collect_predictions() %>%
    dplyr::filter(mtry == best$mtry, min_n == best$min_n) %>%
    dplyr::arrange(., .row)
  cv_rmse <- show_best(tune_res, "rmse")[1, "mean"]
  # plot(cv_predictions$per_ind_3_months, cv_predictions$.pred, xlim = c(0, 1), ylim = c(0, 1))
  # lines(c(0, 1), c(0, 1), lty = "dashed")
 
  # Selecting the Best Model and Final Fitting
  unregister_dopar()
  random_forest_final <- rf_workflow %>%
    finalize_workflow(best)
  final_random_forest_fit <- random_forest_final %>%
    fit(data = juiced)
  fit_OOS_preds <- final_random_forest_fit$fit$fit$fit$predictions
  fit_OOS_preds_df <- data.frame(truth = juiced$per_ind_3_months, estimate = fit_OOS_preds)
  test_rmse <- sqrt(sum((fit_OOS_preds_df$truth - fit_OOS_preds_df$estimate)^2)/length(fit_OOS_preds_df$truth))
  #plot(fit_OOS_preds_df$truth, fit_OOS_preds_df$estimate, xlim = c(0, 1), ylim = c(0, 1)) 
  #lines(c(0, 1), c(0, 1), lty = "dashed")

  # Calculating variable importance
  var_imp <- extract_fit_parsnip(final_random_forest_fit) %>%
    vip(num_features = dim(juiced)[2])
  var_imp <- var_imp$data
  var_imp$id <- i
  var_imp$sampling <- "no_upsampling"

  # Assigning Outputs to Tibble
  iterations[i, "seed"] <- seed
  iterations[i, "model"] <- list(list(final_random_forest_fit))
  iterations[i, "iteration"] <- i
  iterations[i, "juiced"] <- list(list(juiced))
  iterations[i, "best_mtry"] <- best$mtry
  iterations[i, "best_min_n"] <- best$min_n
  iterations[i, "cv_rsme"] <- cv_rmse
  iterations[i, "test_predictions"] <- list(list(fit_OOS_preds))
  iterations[i, "test_rsme"] <- test_rmse
  iterations[i, "importance"] <- list(list(var_imp))
  iterations[i, "recipe"] <- list(list(envt_prepped))
  
  print(paste0("Iteration ", i, " - Seed is ", seed))
  
}

saveRDS(iterations, file = paste0(here("outputs", "random_forest_outputs", "regresstion_repeated_rf_FullData.rds")))

no_ups <- readRDS(here("outputs", "random_forest_outputs", "regresstion_repeated_rf_FullData.rds"))
no_ups_vip_results <- tibble(id = 1, var = "bloop", x = 1, y = 1)
for (i in 1:dim(no_ups)[1]) {
  final_random_forest_fit <- no_ups$model[[i]]
  explainer <- explain_tidymodels(
    model = final_random_forest_fit,
    data = juiced,
    y = as.numeric(juiced$per_ind_3_months),
    verbose = FALSE)
  rem <- which(colnames(explainer$dat) == "country2")
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

x <- bind_rows(no_ups$importance) %>%
  group_by(Variable) %>%
  summarise(mean = mean(Importance),
            sd = sd(Importance),
            se = sd(Importance)/sqrt(n()))
x <- x[rev(order(x$mean)), ]

labs <- gsub("C_", "C", x$Variable)
labs <- gsub("_", "\n", labs)
labs <- gsub("itation", ".", labs)
labs <- gsub("seasonality", "seas.", labs)
labs <- gsub("temperature", "temp.", labs)
labs <- gsub("population", "pop.", labs)

var_imp_plot <- ggplot(x, aes(x= reorder(Variable, mean), y = mean, fill = mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = pmax(0, mean - 1.96 * sd),
                    ymax = mean + 1.96 * sd),
                width = 0.5) +
  scale_fill_continuous(low = "grey", high = "#E14545") +
  theme_bw() +
  scale_x_discrete(labels = rev(labs)) +
  labs(x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7))

preds <- matrix(nrow = 25, ncol = 65)
for (i in 1:25) {
  temp <- no_ups$test_predictions[[i]]
  preds[i, ] <- temp
}

mean_pred <- apply(preds, 2, mean)
pred_df <- data.frame(Actual = juiced$per_ind_3_months, Prediction = mean_pred)
preds_plot <- ggplot(pred_df, aes(x = Prediction, y = Actual)) +
  geom_point() +
  lims(x = c(0, 1), y = c(0, 1)) +
  geom_abline(linetype = "dashed") +
  theme_bw()
cor(mean_pred, juiced$per_ind_3_months)

rf_reg_plot <- cowplot::plot_grid(preds_plot, var_imp_plot, nrow = 2, ncol = 1, rel_widths = c(0.66, 2), align = "h", axis = "b")
ggsave(filename = here("figures/Supp_Figure_RF_Regression.pdf"), plot = rf_reg_plot, width = 6.5, height = 6)

