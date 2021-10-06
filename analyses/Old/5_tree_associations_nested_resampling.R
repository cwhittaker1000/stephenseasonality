#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(tidyverse); library(here); library(zoo); library(forecast); library(TSA); 
library(mgcv); library(GPfit); library(rstan); library(shinystan); library(reshape2); 
library(deSolve); library(parallel); library(matlib); library(matlab); library(pracma); 
library(rstan); library(ggplot2); library(invgamma); library(tictoc); library(DescTools);
library(dismo); library(gbm); library(mltools); library(glmnet); library(caret); library(scales)
library(tidymodels); library(doParallel); library(vip); library(forcats); library(vip)

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

# Creating Test and Training Data for Nested Cross-Validation
set.seed(123)
data <- overall %>% # need to figure out whether to do rf_train or data here 
  dplyr::select(period, country, population_per_1km:worldclim_9) %>%
  dplyr::select(-contains("LC")) %>%
  mutate(country = as.factor(country))
envt_recipe <- recipe(period  ~ ., data = data) %>% 
  step_center(all_predictors(), -country) %>% 
  step_scale(all_predictors(), -country) %>%
  step_dummy(country)
envt_prepped <- prep(envt_recipe, training = data, verbose = TRUE)
juiced <- juice(envt_prepped)
results <- nested_cv(data, 
                     outside = vfold_cv(v = 6), 
                     inside = vfold_cv(v = 5))
perf_metric <- metric_set(yardstick::rmse)

# Random Forest
random_forest <- rand_forest(mtry = tune(), trees = 1000, min_n = tune()) %>%
  set_mode("regression") %>%
  set_engine("ranger", importance = "permutation")
rf_workflow <- workflow() %>%
  add_recipe(envt_recipe) %>%
  add_model(random_forest)
rf_grid <- grid_regular(mtry(range = c(3, 20)), min_n(range = c(2, 10)), levels = 15)

# Running the nested resampling
all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {library(tidymodels)})

set.seed(345)

tuning_results <- map(results$inner_resamples, function(dat, mod_inner, params) {
  tune_grid(object = mod_inner,
            resamples = dat,
            grid = params,
            metrics = perf_metric)
},
rf_workflow, rf_grid)

stopCluster(cl)

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

pooled_inner <- tuning_results %>% 
  bind_rows
best_cost <- function(dat) dat[which.min(dat$.estimate),]

show_best(tuning_results[[1]], "rmse")
show_best(tuning_results[[2]], "rmse")
show_best(tuning_results[[3]], "rmse")
show_best(tuning_results[[4]], "rmse")

best_cost <- function(dat) {
  temp <- dat[which.min(dat$mean),] 
  return(temp)
}

best_cost(tuning_results[[1]])

tuning_results[[1]] %>%
  collect_metrics()

tuning_results[[2]] %>%
  collect_metrics()


tuning_results[[1]]

p <- 
  ggplot(pooled_inner, aes(x = cost, y = .estimate)) + 
  scale_x_continuous(trans = 'log2') +
  xlab("SVM Cost") + ylab("Inner RMSE")

tuning_results[[i]]
best_cost(tuning_results[[i]])

x <- tuning_results[[4]]
x$.metrics

for (i in 1:length(tuning_results))
  p <- p  +
  geom_line(data = tuning_results[[i]], alpha = .2) +
  geom_point(data = best_cost(tuning_results[[i]]), pch = 16, alpha = 3/4)

p <- p + geom_smooth(data = pooled_inner, se = FALSE)
p

