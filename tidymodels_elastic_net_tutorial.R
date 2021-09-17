# Tutorial from https://dnield.com/posts/tidymodels-intro/
# Extra bits from https://juliasilge.com/blog/lasso-the-office/ 

# Load required libraries
library(tidymodels); library(AmesHousing); library(doParallel); library(vip); library(forcats)

# Load data
ames <- make_ames()
ames$Sale_Price <- log(ames$Sale_Price)

## Setting seed
set.seed(100)

# Generate training and testing data 
ames_split <- initial_split(ames, strata = "Sale_Price")
ames_train <- training(ames_split)
ames_test <- testing(ames_split)

# Preprocessing using recipes package - creates list of preprocessing operations that should be performed
#   Skip only applies when we're baking/juicing - which I think happens under the hood
#   hood in predict - basically we still get the right answers out, but if this = FALSE
#   predict fails. See https://github.com/tidymodels/workflows/issues/31.
#   Potential issues re tuning etc, see https://github.com/tidymodels/workflows/issues/47 and 
#   https://github.com/tidymodels/hardhat/issues/129#issuecomment-626386670
#   -> Instead, just manually transform the outcome variable at the outset
ames_rec <- recipe(
  Sale_Price ~ ., 
  data = ames_train) %>% 
  #step_log(Sale_Price, base = 10, skip = FALSE) %>% 
  step_rm(matches("Qual"), matches("Cond")) %>% 
  step_dummy(all_nominal()) %>% 
  step_center(all_predictors()) %>% 
  step_scale(all_predictors()) %>% 
  step_pca(contains("SF"), contains("Area"), threshold = .75) %>% 
  step_nzv(all_predictors())

# Estimates the required parameters for the preprocessing steps (we'll apply these later)
#   During the process of preparing the recipe, each step is estimated via prep and then 
#   applied to the training set using bake before proceeding to the next step.
ames_rec_trained <- prep(ames_rec, training = ames_train, verbose = TRUE)

# Elastic Net Regression preparation
ames_lasso <- linear_reg(penalty = 0.001, mixture = 1) %>% # we'll later optimise these parameters
  set_engine("glmnet")

# Create workflow (container object that aggregates info required to fit and predict from a model i.e. a recipe, a model etc)
# Note: When a recipe is used in a workflow as in your example, the pre-processing functions are not required.
ames_lasso_wfl <- workflow() %>% 
  add_recipe(ames_rec) %>% 
  add_model(ames_lasso)

# Fit the model and predict the values
ames_lasso_fit <- fit(ames_lasso_wfl, ames_train)
predict(ames_lasso_fit, ames_train) 

x <- predict(ames_lasso_fit, ames_train) 
plot(ames_train$Sale_Price, x$.pred)

perf_metrics <- metric_set(rmse, rsq, ccc)

perf_lasso <- ames_lasso_fit %>% 
  predict(ames_train) %>% 
  bind_cols(juice(ames_rec_trained)) %>% 
  perf_metrics(truth = Sale_Price, estimate = .pred)

perf_lasso %>% 
  arrange(.metric)

cv_splits <- vfold_cv(ames_train) # v sets number of splits
cv_splits

cv_eval <- fit_resamples(ames_lasso_wfl, resamples = cv_splits, metrics = perf_metrics)
cv_eval

collect_metrics(cv_eval)
perf_lasso %>% 
  arrange(.metric)


# Model Tuning
ames_mixture <- linear_reg(penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet")
ames_mixture <- linear_reg(penalty = tune(), mixture = 0) %>%  
  set_engine("glmnet")

ames_mixture_wfl <- update_model(ames_lasso_wfl, ames_mixture)

mixture_param <- parameters(penalty(), mixture())

mixture_param <- penalty() 
mixture_param <- range_set(mixture_param, c(-5, 3))

regular_grid <- grid_regular(mixture_param, levels = c(50))

regular_grid %>% 
  ggplot(aes(x = mixture, y = penalty)) +
  geom_point() +
  scale_y_log10()

all_cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {library(tidymodels)})
ames_tune <- tune_grid(
  object = ames_mixture_wfl,
  #preprocessor = ames_rec,
  #model = ames_mixture,
  resamples = cv_splits,
  grid = regular_grid,
  metrics = perf_metrics
)

stopCluster(cl)

# Naive Lasso performance
collect_metrics(cv_eval)
show_best(ames_tune, "ccc")
show_best(ames_tune, "rmse")
show_best(ames_tune, "rsq")

collect_metrics(ames_tune) %>% 
  filter(.metric == "rmse") %>%
  #mutate(mixture = format(mixture)) %>% 
  ggplot(aes(x = penalty, y = mean)) +#, col = mixture)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  geom_vline(xintercept = 0.001, color = "purple", lty = "dotted")

best_mixture <- select_best(ames_tune, metric = "rmse", maximize = FALSE)
best_mixture

ames_mixture_final <- ames_mixture_wfl %>% 
  finalize_workflow(best_mixture) %>% 
  fit(data = ames_train)

tidy_coefs <- ames_mixture_final$fit$fit$fit %>% 
  broom::tidy() %>% 
  filter(term != "(Intercept)") %>% 
  select(-step, -dev.ratio)

delta <- abs(tidy_coefs$lambda - best_mixture$penalty)
lambda_opt <- tidy_coefs$lambda[which.min(delta)]

label_coefs <- tidy_coefs %>% 
  mutate(abs_estimate = abs(estimate)) %>% 
  filter(abs_estimate >= 0.01) %>% 
  distinct(term) %>% 
  inner_join(tidy_coefs, by = "term") %>% 
  filter(lambda == lambda_opt)

tidy_coefs %>% 
  ggplot(aes(x = lambda, y = estimate, group = term, col = term, label = term)) +
  geom_vline(xintercept = lambda_opt, lty = 3) +
  geom_line(alpha = .4) +
  theme(legend.position = "none") +
  scale_x_log10() +
  ggrepel::geom_text_repel(data = label_coefs)

label_coefs

ames_mixture_final %>% 
  predict(ames_test) %>% 
  bind_cols(select(ames_test, Sale_Price)) %>% 
  mutate(Sale_Price = Sale_Price) %>% 
  perf_metrics(truth = Sale_Price, estimate = .pred)
collect_metrics(cv_eval)
perf_lasso %>% 
  arrange(.metric)

ames_mixture_final %>%
  fit(ames_test) %>%
  extract_fit_parsnip() %>%
  vi(lambda = best_mixture$penalty) %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x = Importance, y = Variable, fill = Sign)) +
  geom_col() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y = NULL)

x <- ames_mixture_final %>% 
  predict(ames_test) 
plot(ames_test$Sale_Price, x$.pred)

### SCRAP ####

extract_fit_engine(ames_lasso_fit) # engine specific fit
extract_fit_parsnip(ames_lasso_fit) # parsnip object type version of the fit
extract_mold(ames_lasso_fit)  # preprocessed mold object
extract_spec_parsnip(ames_lasso_fit) # parsnip model specification
extract_preprocessor(ames_lasso_fit) # preprocessing information
extract_recipe(ames_lasso_fit) # returns the relevant recipe

x <- ames_rec_trained %>%
  #bake(new_data = NULL) %>%
  dplyr::slice(1:10) %>%
  dplyr::select(-Sale_Price)

x <- ames_train %>%
  #bake(new_data = NULL) %>%
  #dplyr::slice(1:10) %>%
  dplyr::select(-Sale_Price)

predict(ames_lasso_fit, x)



lm_model <-
  linear_reg() %>%
  set_engine("lm") %>%
  fit(mpg ~ ., data = mtcars %>% dplyr::slice(11:32))

pred_cars <-
  mtcars %>%
  dplyr::slice(1:10) %>%
  dplyr::select(-mpg)

predict(lm_model, pred_cars)


# sfd_grid <- grid_max_entropy(mixture_param, size = 25)
# 
# sfd_grid
# 
# sfd_grid %>% 
#   ggplot(aes(x = mixture, y = penalty)) +
#   geom_point() +
#   scale_y_log10()
