# Loading required libraries
library(data.table); library(devtools); library(ggplot2); library(tidyverse)
library(lubridate); library(rgdal); library(rgeos); library(raster); library(viridis)
library(ggpolypath); library(maptools); library(tidyverse); library(plyr); library(e1071)
library(odin); library(ggpubr); library(viridis); library(Hmisc); library(cowplot)
library(ipred); library(ICDMM); #devtools::install_github("https://github.com/jhellewell14/ICDMM", dependencies = TRUE)
library(scales); library(patchwork)

# Loading functions and mosquito bionomics data
options(scipen = 999)
invisible(sapply(list.files("functions/malaria_model_running/", full.names = TRUE, recursive = TRUE), function(x) source(x)))
bionomics_data <- read.csv("data/bionomic_species_all_LHC_100.csv", stringsAsFactors = FALSE)
stephensi_data <- round(as.data.frame(rbind(colMeans(subset(bionomics_data, species == "stephensi")[, c("Q0", "chi", "bites_Bed", "bites_Indoors")]))), 2)

# Generating the vector of how mosquito density changes over time 
density_start <- 0.1
density_end <- 10
values <- sigmoid(seq(-10, 10, length.out = 365 * 4)) # How long you want introduction to last - here, 10 years
density_vec <- c(rep(density_start, 365 * 1), 
                 pmin((values * (density_end - density_start)) + density_start, density_end),
                 rep(density_end, 365 * 23))
dens_vec_df <- data.frame(vector_density = density_vec, time = seq(1:length(density_vec))/365)

#Pre-loaded seasonal profiles
admin_units_seasonal <- readRDS("data/admin_units_seasonal.rds")
admin_matches <- admin_match(admin_unit = "Fatick", country = "Senegal", admin_units_seasonal = admin_units_seasonal)
loc_seasonality <- seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 1)
scalar_values <- seq(1, 0.05, -0.05)
seasonality_list <- vector(mode = "list", length = length(scalar_values))
for (i in 1:length(scalar_values)) {
  temp <- seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = scalar_values[i])
  temp_df <- data.frame(scalar = scalar_values[i], density = temp, time = 1:length(temp))
  seasonality_list[[i]] <- temp_df
}
seasonal_profiles <- bind_rows(seasonality_list)
seasonal_profiles <- rbind(seasonal_profiles, data.frame(scalar = 0, density = 1, time = 1:length(temp)))

# Running the deterministic malaria model with increasing vector density over time and with/without seasonality
loc_seasonality_list <- vector(mode = "list", length = length(scalar_values) + 1)
for (i in 1:length(scalar_values)) {
  loc_seasonality_list[[i]] <- seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = scalar_values[i])
}
loc_seasonality_list[i + 1] <- NA 

multi_outputs <- sapply(1:(length(scalar_values) + 1), function(x){
  
  set_up_model <- create_r_model_epidemic(odin_model_path = "models/odin_model_seasonality.R", # model file
                                          #Model parameters to work out transmission
                                          init_EIR = 0.000001, # initial EIR from which the endemic equilibria solution is created
                                          #These are the mosquito parameters (bionomics)
                                          Q0 = stephensi_data$Q0, 
                                          chi = stephensi_data$chi, 
                                          bites_Bed = stephensi_data$bites_Bed,
                                          bites_Indoors = stephensi_data$bites_Indoors, 
                                          #These are our seasonal variations
                                          scalar = 1, #This is the scalar shown above in the plots - DOESN'T SEEM TO DO ANYTHING - ASK ARRAN!!!
                                          #custom_seasonality = if(x == 1) loc_seasonality else NA,
                                          custom_seasonality = loc_seasonality_list[[x]], # NA for perennial
                                          #This sets up how long we want to run for and the density vec
                                          time_length = length(density_vec),
                                          density_vec = density_vec)
  
  #Process model formulation
  set_up_model <- set_up_model$generator(user = set_up_model$state, use_dde = TRUE)
  mod_run <- set_up_model$run(t = 1:length(density_vec))
  out <- set_up_model$transform_variables(mod_run)
  model_ran <- as.data.frame(out)
  model_ran <- data.frame(id = x, model_ran)
  
}, simplify = FALSE)

saveRDS(multi_outputs, file = "outputs/malaria_model_running.rds")

# Loading In Saved Runs and Combining Into Overall Data Frame
multi_outputs <- readRDS("outputs/malaria_model_running.rds")
for (i in 1:(length(scalar_values) + 1)) {
  if (i == 1) {
    temp <- multi_outputs[[i]][, c("id", "t", "prev", "Incidence", "mv")]
  } else {
    temp <- rbind(temp, multi_outputs[[i]][, c("id", "t", "prev", "Incidence", "mv")])
  }
}
temp$id <- as.factor(temp$id)

# Calculating Aggregate Quantities E.g. Degree of Seasonality, Time to 2% Etc
calc_incidence_seasonality <- function(input, num_months) {
  days_in_year <- length(input)
  period <- round(days_in_year * num_months/12)
  incidence_vector <- vector(mode = "numeric", length = days_in_year)
  for (i in 1:days_in_year) {
    index <- i:(i+period-1)
    if (sum(index > days_in_year) > 0) {
      index[index>days_in_year] <- index[index>days_in_year] - days_in_year
    }
    incidence_vector[i] <- sum(input[index])/sum(input)
  }
  return(max(incidence_vector))
}

time_to_2_percent <- c()
time_to_2_percent_yearly_average <- c()
seasonality <- c()
for (i in 1:(length(scalar_values) + 1)) {
  temp <- multi_outputs[[i]]$prev
  above_2_percent <- which(temp > 0.02)[1]/365
  yearly_averages <- colMeans(matrix(temp, nrow = 365))
  yearly_average_above_2_percent <- which(yearly_averages > 0.02)[1]
  
  time_to_2_percent <- c(time_to_2_percent, above_2_percent)
  time_to_2_percent_yearly_average <- c(time_to_2_percent_yearly_average, yearly_average_above_2_percent)
  
  seasonality <- c(seasonality, calc_incidence_seasonality(loc_seasonality_list[[i]], 3))
}

# Plotting Seasonal Vector Profiles
cols <- c("#A0CFD3", "#8D94BA", "#9A7AA0", "#87677B")
ex_ind2 <- c(1, 3, 7, 11, 14, 17)
cols_21 <- colorRampPalette(c(cols))
seasonal_vector_plot <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar %in% c(scalar_values[ex_ind2], 0), ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  scale_colour_manual(values = rev(cols_21(length(ex_ind2) + 1))) +
  theme_bw() +
  labs(x = "Time (Days)", y = "Vector Density") +
  theme(legend.position = "none")
annual_vector_plot <- ggplot(dens_vec_df, aes(x = time, y = vector_density)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(x = "Time (Years)", y = "Annual Vector\nDensity") +
  geom_line() +
  lims(x = c(0, 6)) +
  theme(axis.title = element_text(size = 10),
        panel.background = element_blank())
vector_plot <- seasonal_vector_plot + inset_element(annual_vector_plot, 0.01, 0.52, 0.50, 0.97)

# Plotting the Individual Malaria Profiles
ex_ind <- c(1, 7, 14, 21)
example_malaria <- ggplot(temp[temp$id %in% ex_ind, ], aes(x = t, y = 10000 * Incidence, col = factor(id))) +
  geom_line(size = 1.5) +
  lims(x = c(0, 10000)) +
  scale_x_continuous(limits = c(0,10000), expand = c(0, 0)) +
  scale_color_manual(values = cols) +
  facet_wrap(~factor(id), nrow = 4) +
  labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
  theme_bw() + 
  scale_y_continuous(limits = c(-1,14), expand = c(0, 0), position = "left",
                     breaks = c(0, 5, 10)) +
  theme(legend.position = "none", 
        #panel.grid.major.x = element_line(colour="grey", size=0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.line.x = element_line(colour="black", size=0.5),
        axis.line.x.bottom = element_line(colour="black", size=0.5),
        axis.line.y.left = element_line(colour="black", size=0.5),
        panel.spacing = unit(1, "lines")) +
        #axis.line.y.right = element_line(colour="dark grey", size=0.5)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, colour = "dark grey")

test1 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[1]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[1]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test2 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[2]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[2]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test3 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == scalar_values[ex_ind[3]], ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[3]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
test4 <- ggplot(data = seasonal_profiles[seasonal_profiles$scalar == 0, ], aes(x = time, y = 10 * density, colour = factor(scalar))) +
  geom_line(size = 2) +
  lims(y = c(0, 45)) +
  theme_bw() +
  scale_colour_manual(values = cols[4]) +
  labs(x = "", y = "Vector Density") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        plot.margin=grid::unit(c(0,0,-5,-5), "mm"),
        axis.title.y = element_text(size = 8))
malaria_plots <- example_malaria + inset_element(test1, 0.05, 0.82, 0.22, 0.97) +
  inset_element(test2, 0.05, 0.565, 0.22, 0.715) +
  inset_element(test3, 0.05, 0.315, 0.22, 0.465) +
  inset_element(test4, 0.05, 0.06, 0.22, 0.21)

# Plotting Seasonality vs Time to 2%
seasonality_dynamics_df <- data.frame(mal_seas = seasonality, rel_time = time_to_2_percent/min(time_to_2_percent))
seasonal_establishment_plot <- ggplot(seasonality_dynamics_df, aes(x = mal_seas, y = rel_time, col = mal_seas)) +
  scale_colour_gradientn(colours=cols) +
  geom_line(size = 2) +
  theme_bw() +
  labs(x = "% Malaria Incidence In Any 3 Month period",
       y = "Relative Increase In Time Taken\nTo Reach 2% Prevalence") +
  theme(legend.position = "none")

# Overall Plotting
left <- plot_grid(vector_plot, seasonal_establishment_plot, nrow = 2, align = "v", axis = "lrtl")
overall <- plot_grid(left, malaria_plots, rel_widths = c(1, 2))
#7.5 x 15

# (optional) Specifies whether graphs in the grid should be horizontally ("h") 
# or vertically ("v") aligned. Options are "none" (default), "hv" (align in both
# directions), "h", and "v".

# (optional) Specifies whether graphs should be aligned by the left ("l"), 
# right ("r"), top ("t"), or bottom ("b") margins. Options are "none" 
# (default), or a string of any combination of l, r, t, and b in any order 
# (e.g. "tblr" or "rlbt" for aligning all margins). Must be specified if 
# any of the graphs are complex (e.g. faceted) and alignment is specified 
# and desired. See align_plots() for details.

# par(mfrow = c(1, 1))
# plot(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 2), ylab = "Seasonality factor", xlab = "Time (days)", bty = "l", type = "l", col = "red", lwd = 2)
# lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 1), col = "blue", lwd = 2) # scalar controls how flat
# lines(seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = 0.7), col = "black", lwd = 2)

# plot(density_vec, bty = "l", ylab = "Vector density", xlab = "Time (days)")
# abline(v = seq(365, length(density_vec), by = 365), lty = 2)

# a <- ggplot(dens_vec_df, aes(x = time, y = vector_density)) +
#   geom_line() +
#   theme_bw() +
#   geom_rect(aes(xmin = 5, xmax = 6, ymin = 9, ymax = 11), 
#             fill = "grey", alpha = 0.03) +
#   labs(x = "Time (Years)", y = "Average Annual Vector Density") +
#   geom_line() +
#   lims(x = c(0, 6))

# d <- ggplot(temp, aes(x = t, y = 10000 * Incidence, col = factor(id))) +
#   geom_line() +
#   lims(x = c(0, 11000)) +
#   scale_color_manual(values = cols) +
#   theme_bw() + 
#   labs(x = "Time (Days)", y = "Incidence Per 10,000 Population") +
#   theme(legend.position = "right")
# d

# palette(rev(c("#C0E1CB", "#8FD5A6", "#61BA81", "#329F5B", "#1F9151", "#0C8346", "#0D5D56")))
# palette(gray(seq(0,.9,len = 25)))
# par(mfrow = c(3, 1))
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$prev, type = "l", bty = "l", ylab = "Prevalence", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0,0.055), xlim = c(0, 10000))
#   } else {
#     lines(multi_outputs[[i]]$prev, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$Incidence, type = "l", bty = "l", ylab = "Vector Density", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0, max(multi_outputs[[i]]$Incidence)), xlim = c(0, 7500))
#   } else {
#     lines(multi_outputs[[i]]$Incidence, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(multi_outputs[[i]]$mv, type = "l", bty = "l", ylab = "Vector Density", xlab = "Time (days)", 
#          col = palette()[i], ylim = c(0, max(multi_outputs[[i]]$mv)), xlim = c(0, 7500))
#   } else {
#     lines(multi_outputs[[i]]$mv, col = palette()[i])
#   }
# }
# for (i in 1:(length(scalar_values) + 1)) {
#   if (i == 1) {
#     plot(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)), type = "l", bty = "l", ylab = "Mean Annual Vector Density", xlab = "Time (Years)", 
#          col = palette()[i], ylim = c(0, max(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)))))
#   } else {
#     lines(colMeans(matrix(multi_outputs[[i]]$mv, nrow = 365)), col = palette()[i])
#   }
# }