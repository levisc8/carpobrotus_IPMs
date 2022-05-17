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

    rmd <- glue('ipms/Model_Summaries/{vr}/{file_nm}_gam.rmd')
    
    sink(file = rmd)
    
    cat(glue("---\ntitle: '{vr}'"),
        "\noutput:\n  html_document:\n    toc: true\n---\n\n")
    
    
    cat("```{r echo = FALSE}\n\n",
        glue("mod_list <- readRDS('../../Model_Fits/ramet_{vr}_list_krig_gam.rds')"),
        "\n```")
    
    cat('\n\n## WAIC Results\n\n')
    cat("```{r echo = FALSE}\n\n",
        "print(mod_list[[1]]$mod_waic.diffs)",
        "\n\n```")
    
    for(j in 1:13) {
      
      cat("\n\n## ", 
          glue("{names(mod_list[[1]])[j]}"),
          "\n\n") 
      
      if(vr %in% c("repro", "survival")) {
      
        # using 89% intervals, as 10k ESS is needed to compute precise 95%
        # (Kruschke 2014), and we don't want to run our MCMC for that long.
        # Apparently 89% is more stable? (Kruschke 2014). 
        # Kruschke, J. (2014). Doing bayesian data analysis: A tutorial with r,
        # JAGS, and stan. Academic Press.

        
        cat("```{r echo = FALSE}\n\n",
            glue("summary(mod_list[[1]][[{j}]], prob = 0.89)"),
            "\n\n",
            glue("plot(mod_list[[1]][[{j}]], ask = FALSE)\n\n"),
            glue("pp_check(mod_list[[1]][[{j}]],
                      type   = '{pp_type}',
                      group  = 'site',
                      freq   = FALSE,
                      ndraws = 100L)"),
            "\n\n```\n\n")
        
      } else {
        
        cat("\n\n```{r echo = FALSE}\n\n",
            glue("summary(mod_list[[1]][[{j}]], prob = 0.89)"),
            "\n\n",
            glue("plot(mod_list[[1]][[{j}]], ask = FALSE)\n\n"),
            glue("pp_check(mod_list[[1]][[{j}]],
                      type   = '{pp_type}',
                      group  = 'site',
                      ndraws = 100L) + 
                 geom_abline(slope = 1, intercept = 0)"),
            "\n\n```\n\n")
        
      }
      
      # Plot gam smooths for all non-control models
      
      if(j != 1) {
        cat("```{r echo = FALSE}\n\n",
            glue("plot(conditional_smooths(mod_list[[1]][[{j}]],
                                           surface = TRUE), 
                       stype = 'raster',
                       theme = theme_bw(),
                       ask = FALSE)"),
            "\n\n```\n\n")
      }
      
    }
    
    sink()
    
    render(input = rmd)
    
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




# Set default theme for figures
theme.bl <- theme(panel.background = element_rect(fill = NA,
                                                  color = 'black',
                                                  size = 1.25),
                  panel.grid = element_blank(),
                  axis.title.y = element_text(size = 14,
                                              margin = ggplot2::margin(t = 0,
                                                                       l = 5,
                                                                       r = 10, 
                                                                       b = 0)),
                  axis.title.x = element_text(size = 14,
                                              margin = ggplot2::margin(t = 10, 
                                                                       b = 5,
                                                                       l = 0, 
                                                                       r = 0)))

best_mod <- function(mod_list) {
  
  rownames(mod_list[[1]]$mod_waic.diffs)[1]
  
}

mod_clim_terms <- function(mod, clim_terms) {
  
  form <- brmsterms(mod$formula)
  vars <- terms(form$allvars)
  nms  <- attr(vars, "term.labels")
  
  out <- clim_terms[clim_terms %in% nms]
  
  return(out)
}

make_clim_seqs <- function(clim_vars, data) {
  
  temp <- data %>%
    select(!!! syms(clim_vars)) %>%
    mutate(pred_clim_var = NA)
  
  out <- temp[0, ] 
  
  for(i in seq_len(ncol(temp) - 1L)) { # drop "pred_clim_var" from loop
    
    use_seq <- seq(temp[ , i, drop = TRUE] - 0.5,
                   temp[ , i, drop = TRUE] + 0.5, 
                   length.out = 50)
    others  <- temp[ , -c(i)]
    
    all_vars <- list2(!!names(temp)[i] := use_seq, !!! others) %>%
      .[names(out)] %>%
      as.data.frame() %>%
      mutate(pred_clim_var = names(temp)[i])
    
    out <- rbind(all_vars, out)
    
  }
  
  out
  
}


clim_plots <- function(size_data, clim_seqs, mod, clim_terms) {
  
  pred_data <- left_join(size_data, clim_seqs)
  
  preds <- predict(mod, newdata = pred_data)[ , 1] %>%
    cbind(pred = ., pred_data) %>%
    select(pred, log_size, size_qtle, pred_clim_var, !!! syms(clim_terms)) %>%
    pivot_longer(all_of(clim_terms), 
                 names_to = "climate_type", values_to = "clim_val")
  
  plts <- list()
  
  for(i in seq_along(unique(preds$climate_type))) {
    
    clim_nm <- unique(preds$climate_type)[i]
    
    temp_dat <- filter(preds, 
                       climate_type == clim_nm & pred_clim_var == clim_nm) %>%
      group_by(size_qtle) #%>%
      # filter(!duplicated(clim_val))
    
    # if(i == 5) {
    #   leg_plt <- ggplot(temp_dat, 
    #            aes(y     = pred, 
    #                x     = clim_val, 
    #                color = size_qtle)) +
    #       geom_line() +
    #       theme(legend.position = "bottom")
    #   
    #   leg <- get_gg_legend(leg_plt)
    # } 
    
    plts[[i]] <- ggplot(temp_dat, 
                        aes(y     = pred, 
                            x     = clim_val,
                            color = size_qtle)) +
      geom_line() +
      theme.bl + 
      theme_bw() +
      ggtitle(clim_nm, subtitle = unique(size_data$site)) 
    
    
  }
  
  # plts[[7]] <- leg
  
  print(((plts[[1]] + plts[[2]] + plts[[3]]) / 
           (plts[[4]] + plts[[5]] + plts[[6]])) + plot_layout(guides = "collect"))
}


size_clim_plots <- function(size_data, clim_seqs, mod, clim_terms, sites) {
  
  pred_data <- left_join(size_data, clim_seqs) %>%
    filter(site %in% sites)
  
  preds <- predict(mod, newdata = pred_data)[ , 1] %>%
    cbind(pred = ., pred_data) %>%
    select(pred, site, log_size, !!! syms(clim_terms)) %>%
    pivot_longer(all_of(clim_terms), 
                 names_to = "climate_type", values_to = "clim_val")
  
  plts <- list()
  
  for(i in seq_along(unique(preds$climate_type))) {
    
    clim_nm <- unique(preds$climate_type)[i]
    
    temp_dat <- filter(preds, 
                       climate_type == clim_nm) 
    
    # if(i == 3) {
    #   
    #  plts[[i]] <- 
    #    ggplot(temp_dat, 
    #          aes(y = pred, 
    #              x = log_size, 
    #              color = site)) +
    #     geom_line() +
    #     theme.bl + 
    #     theme_bw() +
    #     ggtitle(clim_nm) + 
    #     theme(legend.position = "top")
    #  
    #  leg <- get_gg_legend(plts[[i]])
    #  
    # }
    
    plts[[i]] <- ggplot(temp_dat, 
                        aes(y = pred, 
                            x = log_size, 
                            color = site)) +
      geom_line() +
      theme.bl + 
      theme_bw() +
      ggtitle(clim_nm) 
    
    
  } 
  
  
  print(((plts[[1]] + plts[[2]] + plts[[3]]) / 
           (plts[[4]] + plts[[5]] + plts[[6]])) + plot_layout(guides = "collect"))
}
