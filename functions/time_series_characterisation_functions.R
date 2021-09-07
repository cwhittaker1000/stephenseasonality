# Entropic measure from Nguyen et al- note, requires time series normalised by total
entropic_measure <- function(fitting_output) {
  entropy <- 0
  for (i in 1:length(fitting_output)) {
    entropy <- entropy + fitting_output[i] * log2(fitting_output[i]/(1/length(fitting_output)))
  }
  return(entropy)
}

# Calculates the Peak Distance from January On a Circle
calculate_peak_distance_from_jan <- function(timepoints, fitting_output) {
  peak_index <- which(fitting_output == max(fitting_output))
  peak_time <- timepoints[peak_index]
  back_to_jan <- peak_time - 1
  forwards_to_jan <- 13 - peak_time
  distance_from_jan <- min(back_to_jan, forwards_to_jan)
  return(distance_from_jan)
}

# Calculates the Proportion of Points that Are Above a Specified Threshold (as a Function of the Mean)
points_greater_than_mean <- function(multiple, fitting_output) {
  mean <- mean(fitting_output)
  greater_than_mean <- sum(fitting_output > (multiple * mean))/length(fitting_output)
  return(greater_than_mean)
} 

# Normalise Time Series by Total Density 
normalise_total <- function(x) {
  temp <- na.trim(x)
  total <- sum(temp)
  y <- temp/total
  return(y)
}

# Fits a 1 Component Von Mises Distribution to the Data 
calc_squared_distance_one_component <- function(parameters, normalised_output) {
  x <- seq(1, length(normalised_output), 1)
  num_points <- length(normalised_output)
  mu <- parameters[1]
  k <- parameters [2]
  f <- exp(k * cos((2*pi*x/num_points) - mu))/(2*pi*besselI(k, 0))
  squared_distance <- sum((normalised_output - f)^2)
  return(squared_distance)
}

# Fits a 2 Component Von Mises Distribution to the Data 
calc_squared_distance_two_components <- function(parameters, normalised_output) {
  x <- seq(1, length(normalised_output), 1)
  num_points <- length(normalised_output)
  mu_1 <- parameters[1]
  k_1 <- parameters [2]
  mu_2 <- parameters[3]
  k_2 <- parameters [2] # assume same k for both components
  w <- parameters[4]
  f <- w * (exp(k_1 * cos((2*pi*x/num_points) - mu_1))/(2*pi*besselI(k_1, 0))) + (1 - w) * (exp(k_2 * cos((2*pi*x/num_points) - mu_2))/(2*pi*besselI(k_2, 0)))
  squared_distance <- sum((normalised_output - f)^2)
  return(squared_distance)
}


# Fits 1 and 2 Component Von Mises Distribution & Extracts Respective Governing Parameters
von_mises_fitting <- function(normalised_output, plotting) {
  
  x <- seq(1, length(normalised_output), 1)
  num_points <- length(x)

  # Fitting One Component Von Mises Distribution 
  one_comp_parameters <- c(3, 6) # initial values for Means, Dispersion 
  #calc_squared_distance_one_component(one_comp_parameters, normalised_output)
  one_comp_optim <- optim(one_comp_parameters, calc_squared_distance_one_component, normalised_output = normalised_output, lower = c(0, 0.1), upper = c(2*pi, 20), method = "L-BFGS-B")
  one_f <- exp(one_comp_optim$par[2] * cos((2*pi*x/num_points) - one_comp_optim$par[1]))/(2*pi*besselI(one_comp_optim$par[2], 0))
  
  # Fitting Two Component Von Mises Distribution
  two_comp_parameters <- c(2, 5, 4, 0.33) # initial Values for Means, Dispersions (same for both components) and Weight of Each Component
  #calc_squared_distance_two_components(two_comp_parameters, normalised_output)
  two_comp_optim <- optim(two_comp_parameters, calc_squared_distance_two_components, normalised_output = normalised_output, lower = c(0, 0.1, 0, 0.1, 0), upper = c(2*pi, 10, 2*pi, 1), method = "L-BFGS-B")
  two_f <- two_comp_optim$par[4] * (exp(two_comp_optim$par[2] * cos((2*pi*x/num_points) - two_comp_optim$par[1]))/(2*pi*besselI(two_comp_optim$par[2], 0))) + (1-two_comp_optim$par[4]) * (exp(two_comp_optim$par[2] * cos((2*pi*x/num_points) - two_comp_optim$par[3]))/(2*pi*besselI(two_comp_optim$par[2], 0)))
  
  # Plotting the Output
  if (plotting == TRUE) {
    plot(2*pi*x/36, normalised_output, ylim = c(0, max(c(one_f, two_f, normalised_output))), pch = 20)
    lines(2*pi*x/36, one_f, col = "blue", lwd = 2)
    lines(2*pi*x/36, two_f, col = "red", lwd = 2)
    legend("topleft", legend = c("One Component", "Two Component"), col = c("blue", "red"), lty = 1, cex=0.8)  
  }
  
  # Storing the Optimised Parameters  
  optimised_parameters <- c(one_comp_optim$par, two_comp_optim$par)
  names(optimised_parameters) <- c("1_Mean", "1_K", "2_1_Mean", "2_K", "2_2_Mean", "2_W")
  
  return(optimised_parameters)
}




