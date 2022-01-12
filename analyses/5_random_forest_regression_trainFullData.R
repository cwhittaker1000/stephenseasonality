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
  dplyr::select(per_ind_4_months, country, population_per_1km:mean_temperature_driest_quarter, -LC_190, -LC_210, -mean_temperature_driest_quarter) %>% # LC190 is urban so correlated v strong with PopPer1km
  mutate(country = case_when((country == "Afghanistan" | country == "Djibouti" | 
                                country == "Myanmar" | country == "Pakistan") ~ "aOther",
                             TRUE ~ country)) %>%
  mutate(country = as.factor(country)) %>%
  mutate(country2 = as.factor(country)) # for stratifying later on
data$population_per_1km <- log(data$population_per_1km)

# Prepping Data - Selecting Variables, Creating Recipe 
set.seed(915)
envt_recipe <- recipe(per_ind_4_months ~ ., data = data) %>%
  update_role(country2, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
  step_center(all_numeric(), -contains("per_ind_4_months")) %>% 
  step_scale(all_numeric(), -contains("per_ind_4_months")) %>%
  step_nzv(all_numeric(), -contains("per_ind_4_months"), freq_cut = 8, unique_cut = 25) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
  step_corr(all_numeric(), -contains(c("population_per_1km", "temperature_seasonality", "precipitation_seasonality_cv", "per_ind_4_months")), threshold = 0.50) %>%
  step_dummy(country)
envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
juiced <- juice(envt_prepped)
colnames(juiced)
dim(juiced)
new_envt_recipe <- recipe(per_ind_4_months ~ ., data = juiced) # have to do this way to ensure recipe isn't applied *within* tune-grid
                                                               # when we do this, recipe is applied within each fold, and diff covariates are removed for each
# Checking Correlations 
rem <- c("country_India", "country_Iran", "peaks", "country_peaks", "country2")
rem_index <- which(colnames(juiced) %in% rem)
correlation_matrix <- cor(juiced[, -rem_index])
corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))
corr_vector <- as.vector(correlation_matrix[as.vector(correlation_matrix) != 1])
hist(unique(corr_vector))
mean(unique(corr_vector))

# Setting Up The Random Forest Framework, CV Folds and Performance Metrics to Evaluate
perf_metrics <- metric_set(yardstick::rmse)
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("regression") %>%
  set_engine("ranger", importance = "permutation", seed = 123) # Creates a model specification.  

## No upsampling
cv_splits <- vfold_cv(juiced, v = 6, strata = country2) # cv folds for pre-juiced data (to avoid issues later on with diff variables getting filtered out in each fold)
rf_workflow <- workflow() %>% # Creates a workflow and then adds our model spec and recipe to it. 
  add_recipe(new_envt_recipe) %>%
  add_model(random_forest)
rf_grid <- grid_regular(mtry(range = c(2, dim(juiced)[2] - round(dim(juiced)[2]/4))), min_n(range = c(2, dim(juiced)[1] - (round(dim(juiced)[1]/3)))), levels = 20) # Creates a grid of hyperparameter values to try

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
stopCluster(cl)

# Evaluating Performance
tune_res_perf <- tune_res %>%
  collect_metrics() 
tune_res_perf %>%
  filter(.metric == "rmse") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() + 
  labs(y = "accuracy")

#######################################################################################################
##                                                                                                   ##
##      Selecting the Best Fitting Set of Hyperparameters and Exploring Out-of-Bag Predictions       ##
##                                                                                                   ##
#######################################################################################################

# Exploring Quality of Fit and Selecting the Best Model 
show_best(tune_res, "rmse") # mtry = # predictors randomly sampled
best_rmse <- select_best(tune_res, "rmse")

# Checking performance on the different categories - VERY poor performance on 2 peaks for not upsampled data
predictions <- tune_res %>%
  collect_predictions() %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n) %>%
  dplyr::arrange(., .row)
sqrt(sum((predictions$.pred - predictions$per_ind_4_months)^2)/length(predictions$.pred))
plot(predictions$per_ind_4_months, predictions$.pred, pch = 20, xlim = c(0, 1), ylim = c(0, 1))
lines(0:1, 0:1, lty = 3)

#######################################################################################################
##                                                                                                   ##
##          Finalising Workflow and Re-Fitting to Entire Training Set, Evaluating Test Set           ##
##                                                                                                   ##
#######################################################################################################
# Finalising Workflows With Respective Best Hyperparameter Values  
random_forest_final <- rf_workflow %>%
  finalize_workflow(best_rmse)

# Fitting a Random Forest With These Tuned Hyperparameters to the Entire Training Dataset
final_random_forest_fit <- random_forest_final %>%
  fit(data = juiced) 
fit_OOS_preds <- final_random_forest_fit$fit$fit$fit$predictions
sqrt(sum((fit_OOS_preds - juiced$per_ind_4_months)^2)/length(juiced$per_ind_4_months))
plot(juiced$per_ind_4_months, fit_OOS_preds, pch = 20, xlim = c(0, 1), ylim = c(0, 1))
lines(0:1, 0:1, lty = 3)
cor(juiced$per_ind_4_months, fit_OOS_preds)

# Calculating variable importance
unregister_dopar()
extract_fit_parsnip(final_random_forest_fit) %>%
  vip(geom = "col", num_features = dim(juiced)[2])

# Collecting the evaluation metrics and assessing performance 
explainer <- explain_tidymodels(
  model = final_random_forest_fit,
  data = juiced,
  y = juiced$per_ind_4_months,
  verbose = FALSE)
rem <- which(colnames(explainer$dat) == "country2")
pdp <- model_profile(explainer, variables = colnames(explainer$data)[-rem], N = NULL)
plot(pdp)
