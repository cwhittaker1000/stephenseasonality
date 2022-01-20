
generate_plot_from_cluster_output_runs <- function(shp0_use, processed_runs, fit_to, id){
  
  save_here <- paste0("figs/", id)
  if(!dir.exists(save_here)) dir.create(save_here, recursive = T)
  
  #Plot combined predictions
  png(paste0(save_here, "/", fit_to, "combined_predictions.png"), height = 5, width = 10, units = 'in', res = 300)
  par(oma = c(0, 0, 1, 1), mar = c(0, 0, 2, 1))
  plot(processed_runs[[1]], box = F, axes = F, main = paste0(fit_to, "\nR2: ", round(processed_runs$R2, 3)))
  plot(shp0_use, add = T)
  dev.off()
  
  #Plot of individual runs
  individual_runs <- processed_runs[[4]]
  varimp_all_overall <- processed_runs[[2]]
  overall_metric <- processed_runs[[5]]
  
  png(paste0(save_here, "/", fit_to, "all_model_type_predictions.png"), height = 6, width = 17, units = 'in', res = 300)
  par(oma = c(1, 1, 1, 1), mar = c(0, 0, 2, 0), mfrow = c(plyr::round_any(length(individual_runs)/4, 1, ceiling), min(c(length(individual_runs), 4))))
  sapply(1:length(individual_runs), function(x){
    plot(individual_runs[[x]], box = F, axes = F, main = paste0(fit_to, "\n", names(individual_runs)[x], ": ", round(subset(overall_metric, mod == names(individual_runs)[x])$overall_r2[1], 3)))
    plot(shp0_use, add = T)
  })
  dev.off()
  
  #Plot variable importance
  individual_r2 <- ggplot(data = overall_metric, aes(x = mod, y = Rsquared, color = mod)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(height = 0) +
    theme_minimal() +
    labs(y = "R2", x = "", color = "", title = fit_to) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, round_any(x = max(overall_metric$Rsquared), accuracy = .1, f = ceiling)))
  
  ggsave(paste0("figs/", id, "/", fit_to, "individual_model_R2_by_model_type.png"), individual_r2, height = 4, width = 4, units = "in", dpi = 300)
  
  var_imp <- processed_runs[[2]]
  var_imp_all <- aggregate(x = list(importance = var_imp$importance_weighted),
                           by = list(variable = var_imp$variable),
                           FUN = mean)
  
  var_imp$variable <- factor(var_imp$variable, levels = var_imp_all[order(var_imp_all$importance), ]$variable)

  variable_importance_ggplot2 <- ggplot(data = subset(var_imp, variable %in% rev(levels(var_imp$variable))[1:15]), aes(y = variable, x = importance_weighted, fill = mod)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(y = "", x = "Importance", fill = "", title = fit_to) +
    theme(axis.text.y = element_text(angle = 0),
          legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2))
  
  ggsave(paste0("figs/", id, "/", fit_to, "individual_model_variable_importance.png"), variable_importance_ggplot2, height = 8, width = 5, units = "in", dpi = 300)
  
  
}