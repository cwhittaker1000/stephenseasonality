#######################################################################################################
##                                                                                                   ##
##             Loading Required Libraries, Functions & Preprocessing Catch Data                      ##
##                                                                                                   ##
#######################################################################################################
library(sp); library(raster); library(rgeos); library(rgdal); library(maptools); library(dplyr); 
library(tidyr); library(maps);library(scales); library(here); library(sf); library(ggplot2);
library(zoo); library(cowplot); library(ggspatial); library(rnaturalearthdata); library(rnaturalearth)
library(ggrepel)

# Load functions
source(here("functions", "time_series_characterisation_functions.R"))

# Helper function to get ggplot2 default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Plotting some representative time-series
metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds")) 
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = dim(metadata)[1], ncol = (12 * interpolating_points + 1))
lower_realisation <- matrix(nrow = dim(metadata)[1], ncol = (12 * interpolating_points + 1))
upper_realisation <- matrix(nrow = dim(metadata)[1], ncol = (12 * interpolating_points + 1))
counter <- 1
for (i in 1:(dim(metadata)[1])) {

  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", metadata$id[i], ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", metadata$id[i], ".rds"))
  }  
  
  # Extracting the mean fitted time-series
  timepoints <- STAN_output$timepoints
  all_timepoints <- STAN_output$all_timepoints
  ordered_timepoints <- all_timepoints[order(all_timepoints)]
  MCMC_output <- STAN_output$chain_output
  f <- MCMC_output[, grepl("f", colnames(MCMC_output))]
  f_mean <- apply(f, 2, mean)
  f_lower <- apply(f, 2, quantile, 0.125)
  f_upper <- apply(f, 2, quantile, 0.875)
  mean <- as.numeric(exp(f_mean)[order(all_timepoints)])
  lower <- as.numeric(exp(f_lower)[order(all_timepoints)])
  upper <- as.numeric(exp(f_upper)[order(all_timepoints)])
  mean_realisation[counter, ] <- mean
  lower_realisation[counter, ] <- lower
  upper_realisation[counter, ] <- upper
  counter <- counter + 1
  # plot(mean_realisation[i, ], main = paste0(metadata$id[index], "  ", metadata$country[index]))
  # browser()
}

# Plotting all the Fitted Time-Series
mean <- as.data.frame(mean_realisation)
mean$country <- metadata$country
mean$id <- fit_index
mean_pv <- mean %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "mean")

lower <- as.data.frame(lower_realisation)
lower$country <- metadata$country
lower$id <- fit_index
lower_pv <- lower %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "lower")

upper <- as.data.frame(upper_realisation)
upper$country <- metadata$country
upper$id <- fit_index
upper_pv <- upper %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "upper")

catch_data <- metadata %>%
  filter(id %in% fit_index) %>%
  dplyr::select(id, country, Jan:Dec)
colnames(catch_data) <- c("id", "country", seq(2, 24, length.out = 12)) 
total_raw_catch <- catch_data %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "raw_catch") %>%
  group_by(id) %>%
  summarise(total_raw = sum(raw_catch, na.rm = TRUE))
# total_per_country <- mean %>%
#   pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
#   group_by(id, country) %>%
#   summarise(total = sum(catch, na.rm = TRUE))
raw_catch_data <- catch_data %>%
  pivot_longer(-c(country, id), names_to = "timepoint", values_to = "raw_catch") %>%
  left_join(total_raw_catch, by = "id") %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  dplyr::select(country, id, timepoint, raw_catch)

overall <- mean_pv %>%
  left_join(lower_pv, by = c("id", "country", "timepoint")) %>%
  left_join(upper_pv, by = c("id", "country", "timepoint")) %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  left_join(raw_catch_data, by =  c("id", "country", "timepoint"))

all_ts_plot <- ggplot(data = overall) +
  geom_path(aes(x = timepoint, y = mean, col = factor(country)), size = 2) +
  geom_point(aes(x = timepoint, y = raw_catch, group = factor(id)), col = "black") +
  geom_ribbon(aes(x = timepoint, ymin = lower, ymax = upper, fill = country), alpha = 0.2) +
  facet_wrap(~id, scales = "free_y") +
  scale_y_continuous(limits=c(0, NA), position = "left") +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(2, 24, length.out = 12)) +
  scale_fill_manual(values = gg_color_hue(6)) +
  ylab("Monthly Catch") +
  theme_bw() +
  theme(#legend.position = "bottom",
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.justification = c(1, 0), legend.position = c(1, 0)) +
  guides(col = "none", fill = guide_legend(nrow = 1, ncol = 6, title = ""))

ggsave2(file = here("figures", "Supp_Fig_AllTs.pdf"), plot = all_ts_plot, width = 16, height = 7, dpi = 500)

