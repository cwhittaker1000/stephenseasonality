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
set.seed(234)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(peaks, country, population_per_1km:worldclim_9, -LC_190) %>%
  mutate(country = as.factor(country))
data$peaks <- ifelse(data$peaks == 1, "one", "two")
data$peaks <- as.factor(data$peaks)

# Exploring Correlations in the Data and Removing Highly Correlated Variables
rem <- c("peaks", "country", "LC_62", "LC_152", "LC_202", "LC_80", "LC_71", "LC_160") # variables step_nzv removes below + 
# outcomes we don't want to assess correlation of
rem_index <- which(colnames(data) %in% rem)
clim_index <- grep("worldclim", colnames(data))
correlation_matrix <- cor(data[, -rem_index])
correlation_matrix <- cor(data[, clim_index])

corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))

col <- colorRampPalette(c("blue", "white", "red"))(60)
heatmap(x =correlation_matrix, col = col, symm = TRUE)
heatmap(x =correlation_matrix, symm = TRUE)
col <- colorRampPalette(c("blue", "white", "blue"))(60)
heatmap(x =correlation_matrix, col = col, symm = TRUE,
        cexRow = 0.7, cexCol = 0.7) #, Rowv = NA, Colv = NA)

# Creating Test and Training Data
rf_split <- initial_split(data, prop = 0.85, strata = country)
rf_train <- training(rf_split)
rf_test <- testing(rf_split)

# Prepping Data - Selecting Variables, Creating Recipe
# Things to think about doing:
#   -> Up/downsampling to get more with 2 peaks - step_rose(peaks) for upsampling, need
#      to look up downsampling. 
envt_recipe <- recipe(peaks  ~ ., data = rf_train) %>%
  update_role(country, new_role = "ID") %>% # retained in the data, but not used for model fitting. Role = "predictor", "id" or "outcome" - used differently in model. 
  step_log(population_per_1km) %>%
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) %>%
  step_nzv(all_predictors(), freq_cut = 80/20) %>% # remove variables that are sparse and unbalanced - remove where >80% of variable values are same value
  step_corr(all_predictors(), threshold = 0.65)
  #step_rose(peaks) 
envt_prepped <- prep(envt_recipe, training = rf_train, verbose = TRUE)
juiced <- juice(envt_prepped)

dim(juiced)
rem <- c("country", "peaks")
rem_index <- which(colnames(juiced) %in% rem)
correlation_matrix <- cor(juiced[, -rem_index])
corrplot(correlation_matrix, type="upper", order="hclust", col = brewer.pal(n=8, name="RdYlBu"))

# Generating CV Folds and Selecting the Performance Metric We'll Evaluate Each Time  
cv_splits <- vfold_cv(rf_train, v = 6, strata = peaks) # v sets number of splits
#perf_metrics <- metric_set(yardstick::accuracy)
#perf_metrics <- metric_set(yardstick::roc_auc)


# Setting Up The Boosted Regression Tree Framework
# don't forget to add something to engine to ensure internal brt assessment
# and the external cross-validation are using the same accuracy metric 
#   set_engine("xgboost", objective = 'binary:logistic', eval_metric = 'logloss')
# https://xgboost.readthedocs.io/en/latest/parameter.html

## ACTUALLY I DON'T THINK EVAL_METRIC IS HAVING ANY EFFECT HERE (WHICH MAKES SENSE AS THE RESULTS DIDN'T CHANGE)
# ONLY HAS AN EFFECT IF WITHIN XGBOOST, YOU SET THE ARGUMENT "EVALS" TO CONTAIN SOME DATA. 
# EVALS IS INSTANCE OF TRAINING DATA WITHIN THE XGBOOST TRAINING (I.E. EVALUATED EACH TIME WHILST TRAINING THE SINGLE BRT)
# AND EVAL_METRIC SPECIFIES HOW TO EVALUATE THE MODEL FIT TO EVALS. 
# SEE EVALS DESCRIPTION HERE: https://xgboost.readthedocs.io/en/latest/python/python_api.html
# I THINK THE OBJECTIVE WILL CHANGE THINGS THOUGH, AS IT "PROVIDES GRADIENT INFORMATION": 
# https://xgboost.readthedocs.io/en/latest/tutorials/custom_metric_obj.html
# And this gradient information is what's used by XGBoost (rather than entropy or Gini indices):
# see https://inblog.in/XGBoost-as-Regressor-and-Classifier-4IVqLERTn8 for further information on how this
# is used. Hessian of Loss, so as long as Loss is differentiable, you can probably do something with it.
# But this leads me to believe setting objective *DOES* have an effect. 
xgb_spec <- boost_tree(trees = 1000, 
                       tree_depth = tune(), 
                       min_n = tune(),
                       mtry = tune(),
                       sample_size = tune(),
                       learn_rate = tune()) %>%
  set_mode("classification") %>%
  set_engine("xgboost", importance = "permutation", seed = 123,
             objective = 'binary:logistic', eval_metric = 'auc')
xgb_wf <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(xgb_spec)

set.seed(123)
xgb_grid <- grid_max_entropy(
  tree_depth(c(3L, 12L)), # tree depth
  min_n(c(1L, 10L)), # minimal tip size (I think)
  mtry(c(2L, 15)), # number of predictors used for each tree
  sample_prop(c(0.4, 0.8)), # proportion of data to be used for each tree fit
  learn_rate(c(-2, -1)), # shrinkage parameter determines contribution of each tree to growing model - low learning rate is regularising
  size = 600)

#######################################################################################################
##                                                                                                   ##
##                      Running the Model In Parallel and Evaluating Performance                     ##
##                                                                                                   ##
#######################################################################################################

# Running the Model in Parallel 
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels); library(finetune)})
set.seed(234)
xgb_word_rs <- tune_grid(object = xgb_wf,
                         resamples = cv_splits, 
                         grid = xgb_grid,
                         #metrics = perf_metrics,
                         control = control_resamples(save_pred = TRUE))
stopCluster(cl)

############

# Evaluating Performance
xgb_word_rs %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  dplyr::select(mean, mtry:sample_size ) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter") %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "rmse")

#######################################################################################################
##                                                                                                   ##
##      Selecting the Best Fitting Set of Hyperparameters and Exploring Out-of-Bag Predictions       ##
##                                                                                                   ##
#######################################################################################################

# Exploring Quality of Fit and Selecting the Best Model 
show_best(xgb_word_rs, "roc_auc")
show_best(xgb_word_rs, "accuracy")
best_rmse <- select_best(xgb_word_rs, "roc_auc")

collect_predictions(xgb_word_rs) %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n,
                tree_depth == best_rmse$tree_depth, learn_rate == best_rmse$learn_rate,
                sample_size == best_rmse$sample_size) %>%
  group_by(id) %>%
  roc_curve(peaks, .pred_one) %>%
  autoplot()

# Extracting and plotting raw predictions
raw_predictions <- xgb_word_rs %>% # can also access predictions through xgb_word_rs$.predictions
  collect_predictions()
raw_pred_best_hyp <- raw_predictions %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n,
                tree_depth == best_rmse$tree_depth, learn_rate == best_rmse$learn_rate,
                sample_size == best_rmse$sample_size)
order_raw_pred_best_hyp <- raw_pred_best_hyp %>%
  dplyr::arrange(., .row)

# Below only applies to accuracy, but shows an important phenomenon
# Note: This rmse calculation is NOT the same as reported in show_best exactly. But that's
#       because this is calculating average of each of the individual predictions' errors, rather 
#       average of the average error for each fold
sum(order_raw_pred_best_hyp$.pred_class == order_raw_pred_best_hyp$peaks)/length(order_raw_pred_best_hyp$peaks)
# This DOES give us the same results as in show_best, because we're calculating the average of the 
# per fold performance average.
fold_mean <- order_raw_pred_best_hyp %>%
  group_by(id) %>%
  summarise(prop = sum(.pred_class == peaks)/n())
mean(fold_mean$prop)

xgb_word_rs %>%
  collect_metrics() %>%
  filter(.metric == "accuracy") %>%
  dplyr::filter(mtry == best_rmse$mtry, min_n == best_rmse$min_n,
                tree_depth == best_rmse$tree_depth, learn_rate == best_rmse$learn_rate,
                sample_size == best_rmse$sample_size)

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
xgb_final <- xgb_wf %>%
  finalize_workflow(best_rmse)

# Fitting a Random Forest With These Tuned Hyperparameters to the Entire Training Dataset
#  Note - random_forest_final has same recipe implicitly in there, so will process this data
#         according to that recipe. Important if you've got something in there that changes 
#         number of datapoints so as step_rose for oversampling. 
final_xgb_fit <- xgb_final %>%
  fit(data = rf_train)

# Evaluating Accuracy of Model Fit on Full Training Data - These are OOS Predictions Using Pruned Trees
# CANT FIND THE EQUIVALENT OF THIS FOR XGB OBJECTS JUST YET - QUESTION ASKED AT RSTUDIO COMMUNITY AND STACK OVERFLOW
# WAITING TO HEAR RESPONSES BACK 
final_train_predictions <- final_xgb_fit$fit$fit$fit$predictions
final_train_predictions <- ifelse(final_train_predictions[, 1] > 0.50, "one", "two")
sum(final_train_predictions == rf_train$peaks)/length(rf_train$peaks)

# Predicting On Training Data Using Whole Forest of Trees (Not Good Idea to Do, Just for Comparison With Above)
train_full_forest_pred <- predict(final_xgb_fit, rf_train)
sum(rf_train$peaks == train_full_forest_pred$.pred_class)/length(train_full_forest_pred$.pred_class)
table(rf_train$peaks, train_full_forest_pred$.pred_class) # note good performance on 1 peak, poor performance on 2 peaks

# Predicting On Test Data Using Whole Forest of Trees
test_full_forest_pred <- predict(final_xgb_fit, rf_test)
sum(rf_test$peaks == test_full_forest_pred$.pred_class)/length(test_full_forest_pred$.pred_class)
table(rf_test$peaks, test_full_forest_pred$.pred_class) # note good performance on 1 peak, poor performance on 2 peaks

# ISSUE: model mainly predicts 1s, and because dataset is mainly 1s, does well. 
# CONSIDER: oversampling or undersampling to balance the classes more. 

# NOTE: Don't actually have to do the above of exploring predictions manually
#       as you'll get the same from the last_fit call below. Just another way
#       of exploring model predictions on training data. 

# Final model fit to the training data
#   Note: To get training data (OOS pruned tree) predictions which are identical to final_train_predictions
#         above, run extract_workflow() and then go to $fit$fit$fit$predictions
unregister_dopar()
set.seed(1345)
final_res <- xgb_final %>%
  last_fit(rf_split) 

# Calculating variable importance
extract_workflow(final_res) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")

# Collecting the evaluation metrics and assessing performance 
final_res %>%
  collect_metrics()

sum(final_res$.predictions[[1]]$peaks == final_res$.predictions[[1]]$.pred_class)/length(final_res$.predictions[[1]]$.pred_class)
table(final_res$.predictions[[1]]$peaks, final_res$.predictions[[1]]$.pred_class) # note that still getting poor performance on two peaks

explainer <- explain_tidymodels(
  model = final_xgb_fit,
  data = dplyr::select(rf_train, -peaks),
  y = as.numeric(rf_train$peaks),
  verbose = FALSE
)

pdp_time <- model_profile(
  explainer,
  variables = "LC_10", 
  N = NULL)
plot(pdp_time)


####################################################
####################################################
####################################################
####################################################
####################################################

# note when I use ROSE oversampling, this code doesn't run
# *because* x$fit$fit$fit$predictions has predictions from initial training,
# not from a new fit to the whole rf_train as far as I can make out. 
# Nah don't think so, think it's refitting to rf_train, but running rf_train through the recipe
# temp <- rf_train#[-c(1:10), ]  
# x <- random_forest_final %>%
#   fit(data = temp)
# dim(x$fit$fit$fit$predictions)[1] == dim(temp)[1] # always match, because temp is being run 
#                                                   # through the recipe associated with random_forest_final

# alternative way of calculating variable importance #2
# THIS IS A MODEL NOT A WORKFLOW, THEREFORE DOESN'T HAVE THE 
# RECIPE ASSOCIATED WITH IT; THEREFORE WILL NOT PROCESS
# THE COUNTRY COLUMN, THEREFORE DIFFERNET RESULTS
# final_rf <- finalize_model(random_forest, best_rmse)
# class(final_rf)
# final_rf %>%
#   set_engine("ranger", importance = "permutation", seed = 123) %>%
#   fit(peaks  ~ ., data = rf_train) %>%
#   vip(geom = "point")
# 
# a <- random_forest_final %>% # 1 "A workflow" - does the prepping step and one hot encoding etc
#   fit(data = rf_train) %>%
#   extract_fit_parsnip()
# 
# b <- final_rf %>% # 2 "A Model" - does not do the prepping step of one hot encoding etc
#   set_engine("ranger", importance = "permutation", seed = 123) %>%
#   fit(peaks  ~ ., data = rf_train) 
# 
# a$fit$predictions[1:3, ]
# b$fit$predictions[1:3, ]

# CURRENTLY THESE TWO ARE DIFFERENT TO EACH OTHER - because categorical is rolled into one
# factor for B (instead of one hot encoded), which I think probably explains why (maybe)
# a$fit$forest$independent.variable.names
# b$fit$forest$independent.variable.names
# length(a$fit$variable.importance)
# length(b$fit$variable.importance)

# Extracts the fitted model and test set predictions from the workflow - these are the OOS pruned tree
# predictions and so should match final_train_predictions above.
# test_train_oos_pred <- extract_workflow(final_res)
# test_train_oos_pred$fit$fit$fit$predictions
# pred <- ifelse(test_train_oos_pred$fit$fit$fit$predictions[, 1] > 0.50, "one", "two")
# sum(pred == rf_train$peaks)/length(rf_train$peaks)
# matches the fitting and OOS, pruning based prediction on lines 129 - 136 (approx)

# z <- predict(y, rf_train)
# sum(temp$peaks == z$.pred_class)/length(z$.pred_class)
# matches the fitting based prediction on the original call way on ~ lines 144 and 145 above. 
# both are predicting on the full random_forest
# d <- predict(y, rf_test)
# sum(as.numeric(rf_test$peaks) == d$.pred_class)/length(d$.pred_class)
# matches the above previous predict call on rf_test

