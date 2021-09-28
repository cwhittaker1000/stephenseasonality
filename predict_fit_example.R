# Load required libraries 
library(tidymodels); library(modeldata) 

# Set seed 
set.seed(123)

# Split up data into training and test
data(cells, package = "modeldata")

# Define Model
rf_mod <- rand_forest(trees = 1000) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

# Fit the model to training data and then predict on same training data
rf_fit <- rf_mod %>% 
  fit(class ~ ., data = cells)
rf_training_pred <- rf_fit %>%
  predict(cells, type = "prob")

# Evaluate accuracy 
data.frame(rf_fit$fit$predictions) %>%
  bind_cols(cells %>% select(class)) %>%
  roc_auc(truth = class, PS)

rf_training_pred %>%   
  bind_cols(cells %>% select(class)) %>%
  roc_auc(truth = class, .pred_PS)


# 
# 
# rf_fit$fit$predictions[1, ]
# rf_training_pred[1, ]
# 
# plot(rf_fit$fit$predictions[, "PS"], rf_training_pred$.pred_PS)
# plot(rf_fit$fit$predictions[, "WS"], rf_training_pred$.pred_WS)
# 
# hist(rf_fit$fit$predictions[, "PS"])
# hist(rf_training_pred$.pred_PS)
# 
# sum(abs(rf_fit$fit$predictions[, "PS"] - rf_training_pred$.pred_PS))/length(rf_training_pred$.pred_PS)







