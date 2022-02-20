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

# Load metadata and admin unit geometries
metadata <- readRDS(here("data", "systematic_review_results", "metadata_and_processed_unsmoothed_counts.rds"))
admin2 <- readRDS(here("data", "admin_units", "simplified_admin2.rds"))
admin1 <- readRDS(here("data", "admin_units", "simplified_admin1.rds"))
admin0 <- readRDS(here("data", "admin_units", "simplified_admin0.rds"))

# Helper function to get ggplot2 default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Collating centroids for all of the locations
for (i in 1:nrow(metadata)) {
  
  # Get relevant admin unit and convert to sf object
  if(!is.na(metadata$admin2[i])) {
    metadata_admin <- metadata$admin2[i]
    admin_2_data <- admin2[admin2$NAME_2 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_2_data)
    centroid <- st_centroid(admin_geometry)
    centroid <- centroid[c("continent", "ISO", "NAME_0", "geometry", "NAME_1", "NAME_2")]
  } else {
    metadata_admin <- metadata$admin1[i]
    admin_1_data <- admin1[admin1$NAME_1 == metadata_admin, ]
    admin_geometry <- st_as_sf(admin_1_data)
    centroid <- st_centroid(admin_geometry)
    centroid <- centroid[c("continent", "ISO", "NAME_0", "geometry", "NAME_1")]
    centroid$NAME_2 <- NA
  }
  
  if (i == 1) {
    centroid_metadata <- centroid
  } else {
    centroid_metadata <- rbind(centroid_metadata, centroid)
  }
  
}

# Plotting the centroids of study locations on a background of the below countries
countries <- c("Afghanistan", "Bangladesh", "Bhutan", "Djibouti", "Egypt", "Ethiopia",
               "India", "Iran", "Iraq", "Israel", "Jordan", "Kenya", "Lebanon", "Myanmar", 
               "Nepal", "Oman", "Pakistan", "Saudi Arabia", "Sudan", "Somalia", "South Sudan", "Yemen",
               "China", "Thailand", "Laos", "Cambodia", "Vietnam", "Syria",
               "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Mongolia", "Armenia", "Azerbaijan",
               "Georgia", "Armenia", "Kazakhstan", "Turkey", "Ukraine", "Romania", "Bulgaria", "Libya",
               "Democratic Republic of the Congo", "Central African Republic", "Rwanda", "Burundi", "Tanzania", "Uganda", "Chad",
               "Malaysia", "Singapore", "Indonesia", "Uzbekistan", "Russia")
pres_countries <- c("Afghanistan", "Djibouti", "India", "Iran", "Myanmar", "Pakistan")

red_admin0 <- admin0[admin0$NAME_0 %in% countries, ]
red_admin0
red_admin0$stephensi <- ifelse(red_admin0$NAME_0 %in% pres_countries, "Yes", "No")

centroids <- st_centroid(red_admin0)
points <- as.data.frame(st_coordinates(centroids))
points$NAME_0 <- red_admin0$NAME_0
points <- points[admin0$NAME_0 %in% pres_countries, ]

set.seed(10)
N <- unname(table(metadata$country))
a <- ggplot() + 
  geom_sf(data = red_admin0, fill = ifelse(red_admin0$stephensi == "Yes", gray(0.98), gray(0.85)),
          col = ifelse(red_admin0$stephensi == "Yes", "black", gray(0.6))) +
  geom_sf(data = red_admin0[red_admin0$stephensi == "Yes", ], fill = gray(0.98),
          col = "black") +
  geom_sf(data = st_jitter(centroid_metadata, factor = .008), 
          aes(fill = NAME_0), col = "black", size = 2, shape = 21) +
  geom_label_repel(data = points, aes(x = X, y=Y, label=NAME_0),
                   color = "black", fontface = "bold", max.overlaps = 2,
                   force = 3, nudge_x = c(0, 0, -11, -1, 5, 9),
                   nudge_y = c(6, 4, -5, 7, 10, 9)) +
  theme_bw() +
  scale_fill_manual(values = gg_color_hue(6),
                    labels = c(paste0("Afghanistan (n = ", N[1],")"), 
                               paste0("Djibouti (n = ", N[2],")"),  
                               paste0("India (n = ", N[3],")"), 
                               paste0("Iran (n = ", N[4],")"), 
                               paste0("Myanmar (n = ", N[5],")"), 
                               paste0("Pakistan (n = ", N[6],")"))) +
  coord_sf(ylim = c(0, 45), xlim = c(30, 105), expand = FALSE)  + 
  theme(panel.background = element_rect(fill = "aliceblue"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14, face = "bold"), 
        #axis.title.x = element_text(margin=margin(t=12)), 
        #axis.title.y = element_text(margin=margin(r=12)),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position = "bottom", 
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13, face = "bold")) +
  #xlab("Longitude") + ylab("Latitude") +
  guides(fill = guide_legend(title = "Country",
                             override.aes = list(size = 5)))
   
ggsave2(file = here("figures", "Fig1A.pdf"), plot = a, width = 11, height = 6, dpi = 500)

# Plotting some representative time-series
for_plot_index <- c(4, 57, 26, 53, 34, 31)
interpolating_points <- 2
prior <- "informative"
mean_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
lower_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
upper_realisation <- matrix(nrow = length(for_plot_index), ncol = (12 * interpolating_points + 1))
counter <- 1
for (i in 1:length(for_plot_index)) {
  
  index <- metadata$id[for_plot_index[i]]
  
  # Loading in and processing the fitted time-series
  if (prior == "informative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/informative"), "/inf_periodic_fit_", index, ".rds"))
  } else if (prior == "uninformative") {
    STAN_output <- readRDS(paste0(here("outputs/neg_binom_gp_fitting/uninformative/"), "/uninf_periodic_fit_", index, ".rds"))
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
  plot(mean_realisation[i, ], main = paste0(metadata$id[for_plot_index[i]], "  ", metadata$country[for_plot_index[i]]))
  #browser()
}

# Plotting all of the outputs to see which to feature in Fig 1A plot
mean <- as.data.frame(mean_realisation)
mean$country <- metadata$country[for_plot_index]
mean_pv <- mean %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "mean")

lower <- as.data.frame(lower_realisation)
lower$country <- metadata$country[for_plot_index]
lower_pv <- lower %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "lower")

upper <- as.data.frame(upper_realisation)
upper$country <- metadata$country[for_plot_index]
upper_pv <- upper %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "upper")

catch_data <- metadata[for_plot_index, ] %>%
  dplyr::select(country, Jan:Dec)
total_raw_catch <- catch_data %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
  group_by(country) %>%
  summarise(total_raw = sum(raw_catch))
total_per_country <- mean %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "catch") %>%
  group_by(country) %>%
  summarise(total = sum(catch))
colnames(catch_data) <- c("country", seq(2, 24, length.out = 12)) 
raw_catch_data <- catch_data %>%
  pivot_longer(-country, names_to = "timepoint", values_to = "raw_catch") %>%
  left_join(total_raw_catch, by = "country") %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  dplyr::select(country, timepoint, raw_catch)

overall <- mean_pv %>%
  left_join(lower_pv, by = c("country", "timepoint")) %>%
  left_join(upper_pv, by = c("country", "timepoint")) %>%
  left_join(total_per_country, by = c("country")) %>%
  mutate(timepoint = as.numeric(gsub("V", "", timepoint))) %>%
  left_join(raw_catch_data, by =  c("country", "timepoint"))

b <- ggplot(data = overall) +
  geom_path(aes(x = timepoint, y = mean, col = country), size = 2) +
  geom_point(aes(x = timepoint, y = raw_catch)) +
  geom_ribbon(aes(x = timepoint, ymin = lower, ymax = upper, fill = country), alpha = 0.2) +
  facet_wrap(~country, scales = "free_y") +
  scale_y_continuous(limits=c(0, NA), position = "right") +
  scale_x_continuous(labels = c("J", "F", "M", "A", "M", "J", 
                                "J", "A", "S", "O", "N", "D"),
                     breaks = seq(2, 24, length.out = 12)) +
  ylab("Monthly Catch") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold")) 
b
ggsave2(file = here("figures", "Fig1B.pdf"), plot = b, width = 8, height = 4.5, dpi = 500)

c <- plot_grid(a, b, ncol = 2, align = "h", axis = "b")
ggsave2(file = here("figures", "Fig1_Overall.pdf"), plot = c, width = 15, height = 5.5, dpi = 500)

a2 <- a +
  theme(legend.position = "right") +
  theme(plot.margin = unit(c(1, 0, 1, 0), "cm"))
b2 <- b + 
  scale_y_continuous(limits=c(0, NA), position = "left") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
c2 <- plot_grid(a2, b2, nrow = 2, align = "v", axis = "lr", rel_heights = c(1.3, 1)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
ggsave2(file = here("figures", "Fig1_Overall_Alt.pdf"), plot = c2, width = 8, height = 8.2, dpi = 500)
