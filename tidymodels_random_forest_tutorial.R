# Loading required libraries 
library(tidyverse); library(tidymodels); library(doParallel)

# Read in, process and visualise data
sf_trees <- read_csv("https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2020/2020-01-28/sf_trees.csv")
trees_df <- sf_trees %>%
  mutate(legal_status = case_when(legal_status == "DPW Maintained" ~ legal_status,
                                  TRUE ~ "Other"),
         plot_size = parse_number(plot_size)) %>%
  select(-address) %>%
  na.omit() %>%
  mutate_if(is.character, factor)

trees_df %>%
  ggplot(aes(longitude, latitude, color = legal_status)) +
  geom_point(size = 0.5, alpha = 0.4) +
  labs(color = NULL)

# Exploring which types of trees are maintained by DPW workers or not
trees_df %>%
  count(legal_status, caretaker) %>%
  add_count(caretaker, wt = n, name = "caretaker_count") %>%
  filter(caretaker_count > 50) %>%
  group_by(legal_status) %>%
  mutate(percent_legal = n / sum(n)) %>%
  ggplot(aes(percent_legal, caretaker, fill = legal_status)) +
  geom_col(position = "dodge") +
  labs(fill = NULL, x = "% of trees in each category")

# Generating test and training data
set.seed(123)
trees_split <- initial_split(trees_df, strata = legal_status)
trees_train <- training(trees_split)
trees_test <- testing(trees_split)

# Create recipe
tree_rec <- recipe(legal_status ~ ., data = trees_train) %>%
  update_role(tree_id, new_role = "ID") %>%
  step_other(species, caretaker, threshold = 0.01) %>%
  step_other(site_info, threshold = 0.005) %>%
  step_dummy(all_nominal(), -all_outcomes()) %>%
  step_date(date, features = c("year")) %>%
  step_rm(date) %>%
  themis::step_downsample(legal_status)

# mtry is the number of predictors to sample at each split
# min_n is the number of observations needed to keep splitting nodes
tune_spec <- rand_forest(mtry = tune(), trees = 10, min_n = tune()) %>%
  set_mode("classification") %>%
  set_engine("ranger")

tune_wf <- workflow() %>%
  add_recipe(tree_rec) %>%
  add_model(tune_spec)

set.seed(234)
trees_folds <- vfold_cv(trees_train)

doParallel::registerDoParallel()

set.seed(345)
tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 20
)

tune_res

tune_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, min_n, mtry) %>%
  pivot_longer(min_n:mtry,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")

# 
# all_cores <- parallel::detectCores(logical = FALSE)
# cl <- makePSOCKcluster(all_cores)
# registerDoParallel(cl)
# 
# clusterEvalQ(cl, {library(tidymodels); library(ranger)})
# tune_res <- tune_grid(
#   tune_wf,
#   resamples = trees_folds,
#   grid = 20)
# 
# stopCluster(cl)


rf_grid <- grid_regular(
  mtry(range = c(10, 30)),
  min_n(range = c(2, 8)),
  levels = 5
)

rf_grid

set.seed(456)
regular_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = rf_grid
)

regular_res

regular_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")


doParallel::registerDoParallel()

set.seed(345)
tune_res <- tune_grid(
  tune_wf,
  resamples = trees_folds,
  grid = 20
)

tune_res






