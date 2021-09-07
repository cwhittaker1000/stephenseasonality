# Initial Function to Generate Starting Values for STAN
initial_function <- function(){
  return(list(length_scale_mean = 2.5, length_scale_sd = 1, alpha_mean = 5, alpha_sd = 5, 
              period_mean = 12, period_sd = 8, overdispersion_mean = 2, overdispersion_sd = 3))
}

# Extract Mean Realisations (note the output is ordered such that it's the training points' mean values first, and then the prediction points' values)
extract_mean <- function(chain_output) {
  f <- chain_output[, grepl("f", colnames(chain_output))]
  f_mean <- apply(f, 2, mean)
  negbinom_intensity_mean <- exp(f_mean)
  return(negbinom_intensity_mean)
}

# Function to Run a NegativeBinomial GP Fitting and Store All the Relevant Things
periodic_gp_fit <- function(time_series, interpolating_points, time_series_index, GP_model, 
                            iter, chains, prior_choice, in_progress_output) {
  
  # Processing the Time Series to Remove NAs Either Side and Ensure Timepoints for NAs in Middle Aren't Included
  observed_time_series <- na.trim(as.numeric(time_series))
  observed_time_series <- round(observed_time_series)
  index_for_NAs <- which(is.na(observed_time_series))
  timepoints <- seq(1, length(observed_time_series))
  if(length(index_for_NAs) > 0) {
    timepoints <- timepoints[-index_for_NAs]
    observed_time_series <- observed_time_series[-index_for_NAs]
  } 
  
  # Creating a Series of Timepoints to Interpolate At
  interpolation_timepoints <- round(seq(0, 12, len = 12 * interpolating_points + 1), digits = 2)
  obs_index <- which(interpolation_timepoints %in% timepoints)
  interpolation_timepoints <- interpolation_timepoints[-obs_index]
  
  # Creating a Dataframe of the Data for Input Into STAN 
  STAN_data <- list()
  STAN_data$y <- observed_time_series
  STAN_data$x_train <- timepoints
  STAN_data$N_train <- length(timepoints)
  STAN_data$x_pred <- interpolation_timepoints
  STAN_data$N_pred <- length(interpolation_timepoints)
  
  if (prior_choice == "informative") {
    STAN_data$length_scale_mean <- 2 
    STAN_data$length_scale_sd <- 1
    STAN_data$alpha_mean <- mean(STAN_data$y) 
    STAN_data$alpha_sd <- sd(STAN_data$y)
    STAN_data$period_mean <- 12 
    STAN_data$period_sd <- 6
    STAN_data$overdispersion_mean <- 0
    STAN_data$overdispersion_sd <- 3
  } else if (prior_choice == "uninformative") {
    STAN_data$length_scale_mean <- 0 
    STAN_data$length_scale_sd <- 2 
    STAN_data$alpha_mean <- 0 
    STAN_data$alpha_sd <- sd(STAN_data$y) 
    STAN_data$period_mean <- 12 
    STAN_data$period_sd <- 8 
    STAN_data$overdispersion_mean <- 0
    STAN_data$overdispersion_sd <- 1
  }
  
  
  # Running the Negative Binomial Gaussian Process Fitting Mechanism in STAN - 10,000 Iterations 
  fit <- sampling(GP_model, data = STAN_data, iter = iter, chains = chains, init = initial_function, refresh = iter) #, control = list(adapt_delta = 0.90,  max_treedepth = 15))
  chain_output <- as.matrix(fit)
  
  if (in_progress_output == TRUE) {
    print(paste0("This is Time Series ", time_series_index))
    print(check_hmc_diagnostics(fit))
    Rhat <- summary(fit)
    Rhat <- Rhat$summary
    Rhat <- Rhat[, "Rhat"]
    print(paste0("Time Series ", time_series_index, " has ", sum(Rhat < 0.99 | Rhat > 1.1), " parameters where Rhat is of a concern"))
    
    ordered <- c(timepoints, interpolation_timepoints)
    ordering_indices <- order(ordered)
    ordered <- ordered[ordering_indices]
    mean_output <- extract_mean(chain_output)[ordering_indices]
    
    plot(STAN_data$x_train, STAN_data$y, col = "black", pch = 20, 
         xlim = c(min(c(STAN_data$x_pred, STAN_data$x_train)), max(c(STAN_data$x_pred, STAN_data$x_train))), ylim = c(0, max(c(mean_output, STAN_data$y))))
    lines(ordered, mean_output, pch = 20, type = "l", col = "red", lwd = 2)
  }
  
  output <- list(chain_output = chain_output, timepoints = timepoints, Rhat = Rhat, 
                 all_timepoints = c(timepoints, interpolation_timepoints))
  
  if (prior_choice == "informative") {
    saveRDS(output, file = paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", time_series_index, ".rds"))
  } 
  else if (prior_choice == "uninformative") {
    saveRDS(output, file = paste0(here("outputs/neg_binom_gp_fitting/uninformative"), "/uninf_periodic_fit_", time_series_index, ".rds"))
  }
  
  return(output)
}


mean_realisation_extract <- function(time_series_index, time_series_matrix, prior, print_output) {
  
  observed_time_series <- time_series_matrix[time_series_index, ]
  observed_time_series <- observed_time_series[!is.na(observed_time_series)]
  
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", time_series_index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/uninf_periodic_fit_", time_series_index, ".rds"))
  }
  
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  ordered_timepoints <- all_timepoints[order(all_timepoints)]
  
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  negbinom_intensity_mean <- exp(f_mean)[order(all_timepoints)]
  
  if (print_output == TRUE) {
    plot(timepoints, observed_time_series, col = "black", pch = 20, ylim = c(0, max(c(negbinom_intensity_mean, observed_time_series))))
    lines(ordered_timepoints, negbinom_intensity_mean, pch = 20, type = "l", col = "red", lwd = 2)
    legend("topleft", paste0("Time Series ", time_series_index), bty="n")
  }
  
  return(list(mean = negbinom_intensity_mean, timepoints = ordered_timepoints))
  
}
