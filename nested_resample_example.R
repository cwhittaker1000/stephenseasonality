library(mlbench)
sim_data <- function(n) {
  tmp <- mlbench.friedman1(n, sd = 1)
  tmp <- cbind(tmp$x, tmp$y)
  tmp <- as.data.frame(tmp)
  names(tmp)[ncol(tmp)] <- "y"
  tmp
}

set.seed(9815)
train_dat <- sim_data(100)
large_dat <- sim_data(10^5)

library(tidymodels)
results <- nested_cv(train_dat, 
                     outside = vfold_cv(repeats = 2), 
                     inside = bootstraps(times = 5))
results

library(kernlab)

svm_rmse <- function(object, cost = 1) {
  y_col <- ncol(object$data)
  mod <- 
    svm_rbf(mode = "regression", cost = cost) %>% 
    set_engine("kernlab") %>% 
    fit(y ~ ., data = analysis(object))
  
  holdout_pred <- 
    predict(mod, assessment(object) %>% dplyr::select(-y)) %>% 
    bind_cols(assessment(object) %>% dplyr::select(y))
  rmse(holdout_pred, truth = y, estimate = .pred)$.estimate
}

# In some case, we want to parameterize the function over the tuning parameter:
rmse_wrapper <- function(cost, object) svm_rmse(object, cost)

tune_over_cost <- function(object) {
  tibble(cost = 2 ^ seq(-2, 8, by = 1)) %>% 
    mutate(RMSE = map_dbl(cost, rmse_wrapper, object = object))
}

summarize_tune_results <- function(object) {
  # Return row-bound tibble that has the 25 bootstrap results
  map_df(object$splits, tune_over_cost) %>%
    # For each value of the tuning parameter, compute the 
    # average RMSE which is the inner bootstrap estimate. 
    group_by(cost) %>%
    summarize(mean_RMSE = mean(RMSE, na.rm = TRUE),
              n = length(RMSE),
              .groups = "drop")
}

tuning_results <- map(results$inner_resamples, summarize_tune_results) 

library(scales)

pooled_inner <- tuning_results %>% bind_rows

best_cost <- function(dat) dat[which.min(dat$mean_RMSE),]

p <- 
  ggplot(pooled_inner, aes(x = cost, y = mean_RMSE)) + 
  scale_x_continuous(trans = 'log2') +
  xlab("SVM Cost") + ylab("Inner RMSE")

for (i in 1:length(tuning_results))
  p <- p  +
  geom_line(data = tuning_results[[i]], alpha = .2) +
  geom_point(data = best_cost(tuning_results[[i]]), pch = 16, alpha = 3/4)

p <- p + geom_smooth(data = pooled_inner, se = FALSE)
p

cost_vals <- 
  tuning_results %>% 
  map_df(best_cost) %>% 
  select(cost)

results <- 
  bind_cols(results, cost_vals) %>% 
  mutate(cost = factor(cost, levels = paste(2 ^ seq(-2, 8, by = 1))))

ggplot(results, aes(x = cost)) + 
  geom_bar() + 
  xlab("SVM Cost") + 
  scale_x_discrete(drop = FALSE)

# here, we're basically taking the best set of hyperparameters from each set 
# of hyperparameter tuning cv, then evaluating this on the holdout data in splits
# (in this case 90/10 where 10 is what we evaluate on)
results <- 
  results %>% 
  mutate(RMSE = map2_dbl(splits, cost, svm_rmse))

set.seed(1)
x <- c()
for (i in 1:20) {
  results <- 
    results %>% 
    mutate(RMSE = map2_dbl(splits, cost, svm_rmse))
  x <- c(x, results$RMSE[1])
  print(i)
}

summary(results$RMSE)


x <- c()
for (i in 1:20) {
  object <- results$splits[[i]]
  cost <- results$cost[i]
  
  set.seed(1)
  z <- svm_rmse(object, cost)
  x <- c(x, z)
  
}

plot(x, results$RMSE)

summary(results$RMSE)


object <- results$splits[[1]]
cost <- 2

set.seed(1)
svm_rmse(results$splits[[3]], 2)

svm_rmse <- function(object, cost = 1, seed = 1) {
  set.seed(seed)
  y_col <- ncol(object$data)
  mod <- 
    svm_rbf(mode = "regression", cost = cost) %>% 
    set_engine("kernlab") %>% 
    fit(y ~ ., data = analysis(object))
  
  holdout_pred <- 
    predict(mod, assessment(object) %>% dplyr::select(-y)) %>% 
    bind_cols(assessment(object) %>% dplyr::select(y))
  rmse(holdout_pred, truth = y, estimate = .pred)$.estimate
}


