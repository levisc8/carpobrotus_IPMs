
# Set default theme for figures
theme.bl <- theme(panel.background = element_rect(fill = NA,
                                                  color = 'black',
                                                  size = 1.25),
                  panel.grid = element_blank(),
                  axis.title.y = element_text(size = 14,
                                              margin = margin(t = 0,
                                                              l = 5,
                                                              r = 10, 
                                                              b = 0)),
                  axis.title.x = element_text(size = 14,
                                              margin = margin(t = 10, 
                                                              b = 5,
                                                              l = 0, 
                                                              r = 0)))

surv_mods <- readRDS("Model_Fits/ramet_surv_list_brms.rds")
grow_mods <- readRDS("Model_Fits/ramet_grow_list_brms.rds")
repr_mods <- readRDS("Model_Fits/ramet_repro_list_brms.rds")
flow_mods <- readRDS("Model_Fits/ramet_flower_list_brms.rds")

ramet_surv_slope_int_cor_brm <- surv_mods$slope_int_cor
ramet_grow_slope_int_uncor_brm <- grow_mods$slope_int_uncor
ramet_repro_slope_int_uncor_brm <- repr_mods$slope_int_uncor
ramet_flower_n_slope_int_uncor_brm <- flow_mods$slope_int_uncor

# loop over populations so that each one gets its own pdf

for(i in seq_along(unique(all_ramets$population))) {
  
  pop <- unique(all_ramets$population)[i]
  
  ramet_temp <- dplyr::filter(all_ramets, population == pop)

  max_flow <- max(ramet_temp$flower_n, na.rm = TRUE)
  temp_min <- min(c(ramet_temp$log_size, ramet_temp$log_size_next), na.rm = TRUE)
  temp_max <- max(c(ramet_temp$log_size, ramet_temp$log_size_next), na.rm = TRUE)
  pred_seq <- seq(temp_min, temp_max, length.out = dim(ramet_temp)[1])
  
  temp_pred_fl <- expand.grid(list(log_size = pred_seq, 
                                   flower_col = unique(ramet_temp$flower_col),
                                   population = pop)) %>%
    .[complete.cases(.), ]
  
  temp_pred_other <- expand.grid(list(log_size   = pred_seq,
                                      population = pop)) %>%
    .[complete.cases(.), ]
  
  flow_pred <- predict(ramet_flower_n_slope_int_uncor_brm,
                       newdata = temp_pred_fl)[ , 1] %>%
    cbind(temp_pred_fl) %>% 
    setNames(c('estimate',
               'log_size',
               'flower_col',
               'population'))
  
  surv_pred <- predict(ramet_surv_slope_int_cor_brm,
                       newdata = temp_pred_other)[ , 1] %>%
    cbind(temp_pred_other) %>%
    setNames(c(
      "estimate",
      "log_size",
      "population"
    ))
  
  grow_pred <- predict(ramet_grow_slope_int_uncor_brm,
                       newdata = temp_pred_other)[ , 1] %>%
    cbind(temp_pred_other) %>%
    setNames(c(
      "estimate",
      "log_size",
      "population"
    ))
  
  repr_pred <- predict(ramet_repro_slope_int_uncor_brm,
                       newdata = temp_pred_other)[ , 1] %>%
    cbind(temp_pred_other) %>%
    setNames(c(
      "estimate",
      "log_size",
      "population"
    ))
  
  n_ram_flower <- sum(ramet_temp$repro)
  n_ram_tot    <- dim(ramet_temp)[1]

  country <- all_sites$Country[all_sites$Site == pop][1]
  
  recr_temp <- filter(recruits, population == pop)
  
  if(!dir_exists(glue('Figures/{country}/{pop}'))) {
    dir_create(glue('Figures/{country}/{pop}'))
  }
  
  ramet_hist <- ggplot(ramet_temp, aes(x = log_size)) + 
    geom_histogram(bins = 20,
                   color = 'red',
                   fill = NA,
                   size = 2) + 
    theme.bl + 
    scale_y_continuous('# of Ramets') + 
    scale_x_continuous(paste('Size of Ramets (n = ',
                             n_ram_tot,
                             ')', 
                             sep = "")) + 
    ggtitle(pop)
  
  n_recr_bins <- ifelse(nrow(recr_temp) > 0,
                        round(dim(recr_temp)[1] / 3),
                        1)
  
  recr_hist <- ggplot(recr_temp, aes(x = log_size_next)) +
    geom_histogram(bins = n_recr_bins,
                   color = 'red',
                   fill = NA,
                   size = 2) +
    theme.bl +
    scale_y_continuous("") +
    scale_x_continuous(paste('Size of Recruits (n = ',
                             nrow(recr_temp),
                             ')', 
                             sep = ""))
  
  ramet_grow_plot <- ggplot(ramet_temp,
                            aes(x = log_size,
                                y = log_size_next)) + 
    geom_point(color = 'black',
               size = 2) +
    geom_line(data = grow_pred, 
              aes(x = log_size,
                  y = estimate),
              color = 'red',
              linetype ='dashed',
              size = 1.5, 
              alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = 'grey', size = 1.5) + 
    theme.bl +
    scale_y_continuous("Size (t + 1)") + 
    scale_x_continuous("Size (t)")
  
  ramet_surv_plot <- ggplot(ramet_temp,
                             aes(x = log_size,
                                 y = alive)) + 
    geom_point(color = 'red',
               size = 1.25) + 
    geom_line(data = surv_pred,
              aes(y = estimate,
                  x = log_size),
              color = 'red',
              linetype = 'dashed',
              size = 1.1,
              alpha = 0.5) + 
    theme.bl + 
    scale_y_continuous('Pr(Survival)') + 
    scale_x_continuous(paste('Size of Ramets (n = ',
                             n_ram_tot,
                             ')',
                             sep = ""))
  
 
  ramet_fec_plot <- ggplot(ramet_temp, 
                           aes(x = log_size,
                               y = flower_n)) + 
    geom_point(aes(color = flower_col)) + 
    geom_line(data = flow_pred,
              aes(x = log_size, 
                  y = estimate,
                  color = flower_col
              )) +
    theme.bl + 
    theme(legend.direction = 'vertical') + 
    scale_y_continuous('# of Flowers',
                       limits = c(0, max_flow + (max_flow/10))) + 
    scale_x_continuous(paste('Size of Ramets (n = ', 
                             n_ram_flower,
                             ')',
                             sep = "")) 
  
  fec_legend <- get_gg_legend(ramet_fec_plot)
  
  ramet_repro_plot <- ggplot(ramet_temp,
                             aes(x = log_size,
                                 y = repro)) + 
    geom_point(color = 'black',
               size = 2) + 
    geom_line(data = repr_pred,
              aes(y = estimate, x = log_size),
              color = 'red',
              linetype = 'dashed',
              size = 1.5,
              alpha = 0.5) + 
    theme.bl + 
    scale_y_continuous('Pr(Reproductive)') + 
    scale_x_continuous("")
  
  
  pdf(glue('Figures/{country}/{pop}/Preliminary_Plots.pdf'),
      width = 8,
      height = 8)
  grid.arrange(ramet_hist, recr_hist,
               ramet_grow_plot, ramet_surv_plot,
               ramet_fec_plot, ramet_repro_plot,
               nrow = 3, ncol = 2)
  dev.off()
  
  png(glue('Figures/{country}/{pop}/Preliminary_Plots.png'),
      width = 8,
      height = 8,
      units = 'in',
      res = 72)
  
  grid.arrange(ramet_hist, recr_hist,
               ramet_grow_plot, ramet_surv_plot,
               ramet_fec_plot, ramet_repro_plot,
               nrow = 3, ncol = 2)
  
  dev.off()
  
  
}
