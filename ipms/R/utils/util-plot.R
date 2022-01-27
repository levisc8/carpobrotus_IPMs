# Plotting utilities


get_gg_legend<-function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



plot_models <- function(mod_list, vr) {
  
  pp_type <- switch(vr,
                    "repro"         = "bars_grouped",
                    "survival"      = "bars_grouped",
                    "log_size_next" = "scatter_avg_grouped",
                    "flower_n"      = "scatter_avg_grouped")
  
  
  for(i in seq_along(mod_list)) {
    
    use_mod <- mod_list[[i]]
    
    file_nm <- names(mod_list)[i]
    
    
    if(vr == "log_size_next") vr <- "growth"
    
    pdf(glue("ipms/Model_Summaries/Trace_Plots/{vr}/{file_nm}.pdf"))
    
    for(j in 1:3) {
      
      plot(use_mod[[j]],
           ask = FALSE)
      
      if(vr %in% c("repro", "survival")){
        
        p <- pp_check(use_mod[[j]],
                      type   = pp_type,
                      group  = "site",
                      freq   = FALSE,
                      ndraws = 100L)
      } else {
        
        p <- pp_check(use_mod[[j]],
                      type   = pp_type,
                      group  = "site",
                      ndraws = 100L)
        
      }
      
      print(p)
      
    }
    
    dev.off()
    
    
    sink(file = glue('ipms/Model_Summaries/{vr}/{file_nm}.txt')) 
    cat('Size only\n\n *********************\n\n')
    print(summary(use_mod[[1]]))
    cat(glue('\n\n*********************\n\nSize + {file_nm}'),
        '\n\n*********************\n\n')
    print(summary(use_mod[[2]]))
    cat(glue('\n\n*********************\n\nSize * {file_nm}'),
        '\n\n*********************\n\n')
    print(summary(use_mod[[3]]))
    cat('\n\n*********************\n\nWAIC Results\n\n')
    print(use_mod[[4]])
    
    cat('\n\nEnd output')
    sink()
    
  }
  
}

.pred_data <- function(native, min_z, max_z, site) {
  
  switch(native,
         "no"  = .pred_clim_data(min_z, max_z, site),
         "yes" = .pred_nat_data(min_z, max_z, site))
}

.pred_nat_data <- function(min_z, max_z, site) {
  
  nat_data <- data.frame(
    log_size = seq(min_z, max_z, length.out = 100),
    native   = rep(1, 100)
    
  )
  inv_data <- data.frame(
    log_size = seq(min_z, max_z, length.out = 100),
    native   = rep(0, 100)
  )
  
  all_data <- list()
  
  for(i in seq_along(site)) {
    
    if(site[i] %in% c("Melkboss",      "Vogelgat", 
                      "St_Francis",    "Struisbaai",
                      "Springfontein", "Rooisand")) {
      
      all_data[[i]] <- cbind(nat_data, site = site[i])
    } else {
      all_data[[i]] <- cbind(inv_data, site = site[i])
    }
    
  }
  
  do.call(rbind, all_data)
}

.pred_clim_data <- function(min_z, max_z, site) {
  
  pred_data <- data.frame(
    log_size     = seq(min_z, max_z, length.out = 100),
    mat_rec      = seq(-1.8, 2.3, length.out = 100),
    map_rec      = seq(-1.32, 1, length.out = 100),
    t_seas_rec   = seq(-1.45, 0.8, length.out = 100),
    p_seas_rec   = seq(-1.9, 2.1, length.out = 100),
    t_co_qu_rec  = seq(-1.04, 1.35, length.out = 100)
  ) 
  
  all_data <- list()
  
  for(i in seq_along(site)) {
    all_data[[i]] <- cbind(pred_data, site = site[i])
  }
  
  do.call(rbind, all_data)
}

plot_preds <- function(models, vr, native = "no") {
  
  z <- c(models[[1]]$size_only$data$log_size,
         models[[1]]$size_only$data$log_size_next)
  site <- unique(models[[1]]$size_only$data$site)
  
  min_z <- min(z, na.rm = TRUE)
  max_z <- max(z, na.rm = TRUE)
  
  pred_data <- .pred_data(native, min_z, max_z, site)
  
  temp <- list()
  
  all_pred <- for(i in seq_along(models)) {
    
    clim_nm <- names(models)[i]
    use_mods <- models[[i]][1:3]
    
    preds <- lapply(use_mods,
                    function(x, pred_data) {
                      temp <- predict(x, newdata = pred_data)
                      cbind(pred_data, temp)
                    },
                    pred_data = pred_data)
    
    for(j in seq_along(preds)) {
      preds[[j]] <- as.data.frame(preds[[j]]) %>%
        mutate(mod_type = names(preds)[j])
    } 
    
    temp[[i]] <- do.call(rbind, preds) %>% 
      mutate(clim = clim_nm)
    
  }
  
  all_data <- do.call(rbind, temp) #%>%
  
  if(native == "no") {
    all_data <- filter(site %in% c("Melkboss", "Rough_Island",
                                   "Foxton", "Praia_de_Areao"))
    
    use_col <- quo(site)
    use_line <- quo(mod_type)
  } else {
    
    all_data <- mutate(all_data,
                       native = case_when(
                         native == 1 ~ "Native",
                         TRUE ~ "Invasive"
                       )) 
    use_col <- quo(site)
    use_line <- quo(native)
  }
  
  
  
  plt <- ggplot(all_data, 
                aes(x = log_size,
                    y = Estimate)) +
    geom_line(aes(linetype = !! use_line, 
                  color = !! use_col),
              size = 1.2,
              alpha = 0.6) +
    theme_bw()
  
  if(native == "no") {
    plt <- plt +
      facet_wrap(~clim, 
                 scales = "free") +
      scale_color_discrete(guide = "none")
  } else {
    plt <- plt +
      facet_wrap(~mod_type, scales = "free")
  }
  
  pdf(glue("ipms/Figures/vr_models/{vr}_model_predictions.pdf"))
  print(plt)
  dev.off()
  
}
