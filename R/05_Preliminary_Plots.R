
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

# loop over populations so that each one gets its own pdf

for(i in seq_along(unique(all_ramets$population))) {
  
  pop <- unique(all_ramets$population)[i]
  
  ramet_temp <- dplyr::filter(all_ramets, population == pop)
  genet_temp <- dplyr::filter(all_genets, population == pop)
  
  n_ram_flower <- sum(ramet_temp$repro)
  n_gen_flower <- sum(genet_temp$repro)
  n_ram_tot <- dim(ramet_temp)[1]
  n_gen_tot <- dim(genet_temp)[1]
  
  genet_temp$clean_bin <- cut(genet_temp$log_size, 
                              breaks = round(dim(genet_temp)[1]/4)) %>%
    clean_bins()
  
  
  country <- all_sites$Country[all_sites$Site == pop][1]
  
  if(!dir_exists(glue('Figures/{country}/{pop}'))) {
    dir_create(glue('Figures/{country}/{pop}'))
  }
  
  ramet_hist <- ggplot(ramet_temp, aes(x = clean_bin)) + 
    geom_bar(width = 0.05,
             fill = 'red') + 
    theme.bl + 
    scale_y_continuous('# of Ramets') + 
    scale_x_continuous(paste('Size of Ramets (n = ',
                             n_ram_tot,
                             ')', 
                             sep = "")) + 
    ggtitle(pop)
  
  genet_hist <- ggplot(genet_temp, aes(x = clean_bin)) + 
    geom_bar(width = 0.2,
             fill = 'blue') + 
    theme.bl + 
    scale_y_continuous('# of Genets') + 
    scale_x_continuous(paste('Size of Genets (n = ',
                             n_gen_tot,
                             ')', 
                             sep = ""))
  
  ramet_fec_plot <- ggplot(ramet_temp, 
                           aes(x = log_size,
                               y = flower_n)) + 
    geom_point(color = 'red') + 
    geom_line(aes(y = flower_pred),
              color = 'red') +
    theme.bl + 
    scale_y_continuous('# of Flowers') + 
    scale_x_continuous(paste('Size of Ramets (n = ', 
                             n_ram_flower,
                             ')',
                             sep = "")) 
  
  genet_fec_plot <- ggplot(genet_temp,
                           aes(x = log_size,
                               y = flower_n)) + 
    geom_point(color = 'blue') + 
    geom_line(aes(y = flower_pred),
              color = 'blue') +
    theme.bl + 
    scale_y_continuous('# of Flowers') + 
    scale_x_continuous(paste('Size of Genets (n = ',
                             n_gen_flower,
                             ')',
                             sep = ""))
  
  ramet_repro_plot <- ggplot(ramet_temp,
                             aes(x = log_size,
                                 y = repro)) + 
    geom_point(color = 'red',
               size = 1.25) + 
    geom_line(aes(y = repro_pred),
              color = 'red',
              linetype = 'dashed',
              size = 1.1,
              alpha = 0.5) + 
    theme.bl + 
    scale_y_continuous('Pr(Reproductive)') + 
    scale_x_continuous(paste('Size of Ramets (n = ',
                             n_ram_tot,
                             ')',
                             sep = ""))
  
  
  
  genet_repro_plot <- ggplot(genet_temp,
                             aes(x = log_size,
                                 y = repro)) + 
    geom_point(color = 'blue',
               size = 1.25) + 
    geom_line(aes(y = repro_pred),
              color = 'blue',
              linetype = 'dashed',
              size = 1.1,
              alpha = 0.5) + 
    theme.bl + 
    scale_y_continuous('Pr(Reproductive)') + 
    scale_x_continuous(paste('Size of Genets (n = ',
                             n_gen_tot,
                             ')',
                             sep = ""))
  
  pdf(glue('Figures/{country}/{pop}/Preliminary_Plots.pdf'),
      width = 8,
      height = 8)
    grid.arrange(ramet_hist, genet_hist,
                 ramet_fec_plot, genet_fec_plot,
                 ramet_repro_plot, genet_repro_plot,
                 nrow = 3, ncol = 2)
  dev.off()
  
  png(glue('Figures/{country}/{pop}/Preliminary_Plots.png'),
      width = 8,
      height = 8,
      units = 'in',
      res = 72)
    grid.arrange(ramet_hist, genet_hist,
                 ramet_fec_plot, genet_fec_plot,
                 ramet_repro_plot, genet_repro_plot,
                 nrow = 3, ncol = 2)
  dev.off()
  

}
