# Dependency file

library(finch)
library(purrr)
library(tidyr)
library(rlang)
library(raster)
library(dplyr)
library(CoordinateCleaner)
library(ggplot2)
library(GGally)
library(ggforce)
library(bayesplot)
library(glue)
library(dismo)
library(sf)
library(fs)
library(brms)
library(rstan)
library(gridExtra)

# Helpers for model diagnostics

# Find observations w/ pareto-K's above a certain threshold value from LOO object
find_prob_obs <- function(loo_obj, thresh) {
  
  return(which(loo_obj$diagnostics$pareto_k > thresh))
  
}

# Get exact row IDs for a data set who's pareto-k values are above a
# certain threshold
get_prob_rows <- function(model_obj, loo_obj, thresh) {
  
  stopifnot(length(model_obj) == length(loo_obj))
  
  # index of values for the model_obj@frame slot
  data_ind <- lapply(loo_obj, FUN = find_prob_obs, thresh = thresh)
  
  
  out <- lapply(seq_len(length(model_obj)), function(x) {
    model_obj[[x]]$data[data_ind[[x]], ]
  })
  
  return(out)
  
}

# Horseshoe prior model helpers

hs_diagnostic_plots <- function(clim, vr, func) {
  
  mod_path <- glue("repro_analysis/model_fits/shrink_mods/vr_model_{clim}_{vr}_{func}.rds")
  
  vr_mod <- readRDS(mod_path)
  
  out_path <- glue("repro_analysis/Manuscript/Figures/shrink_mod_{clim}_{vr}_{func}.pdf")
  
  switch(vr,
         "repro"    = .hs_repro_plot(vr_mod, out_path),
         "flower_n" = .hs_flower_plot(vr_mod, out_path))

  invisible(TRUE)
  
}

.hs_repro_plot <- function(model, out_path) {
  
  pdf(out_path)
  plot(model,
       ask = FALSE)
  
  
  p <- pp_check(model,
                type     = 'bars',
                nsamples = 100L,
                freq     = FALSE)
  
  print(p)
  
  p <- pp_check(model,
                type     = 'dens_overlay',
                nsamples = 100L)
  
  print(p)
  
  dev.off()
  invisible(TRUE)
  
}

.hs_flower_plot <- function(model, out_path) {
  
  pdf(out_path)
  plot(model,
       ask = FALSE)
  
  
  p <- pp_check(model,
                type     = 'bars',
                nsamples = 100L,
                freq     = FALSE) + 
    xlim(c(0, 20))
  
  print(p)
  
  p <- pp_check(model,
                type     = 'rootogram',
                nsamples = 100L) 
  
  print(p + xlim(c(0,50)))
  
  print(p + xlim(c(50, 5000)))
  
  dev.off()
  invisible(TRUE)
  
}

# Plotting functions for size X climate image plots

gg_image_plot <- function(data, x_vals,
                          y_vals, z_var,
                          dummy_x, dummy_y) {
  
  dummy_x <- enquo(dummy_x)
  dummy_y <- enquo(dummy_y)
  z_var <- enquo(z_var)
  
  x_vals <- seq(min(x_vals), max(x_vals), length.out = 5) %>% round(digits = 1)
  y_vals <- seq(min(y_vals), max(y_vals), length.out = 5) %>% round(digits = 2)
  
  out <- ggplot(data) +
    geom_tile(
      aes(
        x    = !! dummy_x,
        y    = !! dummy_y,
        fill = !! z_var 
      )
    ) +
    geom_contour(
      aes(x = !! dummy_x,
          y = !! dummy_y,
          z = !! z_var),
      color = "black",
      size  = 0.9
    ) +
    scale_x_continuous("SD(Global Climate Envelope Mean)",
                       breaks = seq(1, 50, length.out = 5),
                       labels = x_vals) +
    scale_y_continuous("Ln(Surface Area)",
                       breaks = seq(1, 50, length.out = 5),
                       labels = y_vals) +
    scale_fill_gradient("",
                        low = "red",
                        high = "yellow") +
    theme(
      panel.background = element_blank()
    )
  
  return(out)
}

get_gg_legend<-function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plot_raw_data <- function(data, preds, ...) {
  
  resp <- rlang::enquo(resp)
  preds <- rlang::enquos(...)
  
  use_data <- data %>%
    dplyr::select(id, population, !! resp, log_size, !!! preds, native) %>%
    tidyr::gather(key = "var", value = "value", -c(id, population, !! resp, native)) %>%
    dplyr::filter(!is.na(!!resp))
  
  ggplot2::ggplot(use_data,
                  ggplot2::aes(x = value,
                               y = !! resp)) +
    ggplot2::geom_point(ggplot2::aes(color = native)) +
    ggplot2::facet_wrap(~var,
               ncol = 1, nrow = 2,
               scales = "free_x")
  
}


pp_bars <- function(object, nsamples, group, freq = FALSE) {
  
  y <- object$data$flower_n
  group <- object$data$population
  yrep <- posterior_predict(object, nsamples = nsamples, group = group)
  
  pp_check_bars(y, yrep, group, freq = freq)
  
}

pp_check_bars <- function (y, yrep, group, ..., facet_args = list(), prob = 0.9, 
                           width = 0.9, size = 1, fatten = 3, freq = FALSE){

  alpha <- (1 - prob)/2
  probs <- sort(c(alpha, 0.5, 1 - alpha))
  yrep_data <- bayesplot:::ppc_bars_yrep_data(y, yrep, probs, freq = freq, 
                                              group = group)
  .pp_check_bars(y_data = data.frame(y, group), yrep_data, grouped = TRUE, 
            facet_args = facet_args, width = width, size = size, 
            fatten = fatten, freq = freq)
  
}


.pp_check_bars <- function (y_data, yrep_data, 
                            facet_args = list(), grouped = FALSE, 
                            width = 0.9, size = 1, fatten = 3, freq = TRUE) {
  
  for(i in 1:5) {
  graph <- ggplot() + 
    geom_bar(data = y_data,
             mapping = if (freq) 
               aes_(x = ~y, fill = "y")
             else 
               aes_(x = ~y, 
                    y = ~..prop.., 
                    fill = "y"), 
             color = bayesplot:::get_color("lh"), 
             width = width) + 
    geom_pointrange(data = yrep_data,
                    mapping = aes_(x = ~x, 
                                   y = ~mid, ymin = ~lo,
                                   ymax = ~hi,
                                   color = "yrep"), 
                    size = size, fatten = fatten) +
    scale_fill_manual("", 
                      values = bayesplot:::get_color("l")) + 
    scale_color_manual("", values = bayesplot:::get_color("dh")) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 2)) + 
    labs(x = NULL,
         y = if (freq) 
           "Count"
         else 
           "Proportion") +
    ggforce::facet_wrap_paginate(
      facets = vars(group),
      nrow = 2,
      ncol = 2,
      page = i,
      scales = "free_y") +
    xlim(c(0, 50))
  

  print(graph)
  
  }
  
  w_y <- filter(y_data, group == "Whirinaki")
  w_y_rep <- filter(yrep_data, group == "Whirinaki")
  
  graph <- ggplot() + 
    geom_bar(data = w_y,
             mapping = if (freq) 
               aes_(x = ~y, fill = "y")
             else 
               aes_(x = ~y, 
                    y = ~..prop.., 
                    fill = "y"), 
             color = bayesplot:::get_color("lh"), 
             width = width) + 
    geom_pointrange(data = w_y_rep,
                    mapping = aes_(x = ~x, 
                                   y = ~mid, ymin = ~lo,
                                   ymax = ~hi,
                                   color = "yrep"), 
                    size = size, fatten = fatten) +
    scale_fill_manual("", 
                      values = bayesplot:::get_color("l")) + 
    scale_color_manual("", values = bayesplot:::get_color("dh")) +
    guides(color = guide_legend(order = 1), 
           fill = guide_legend(order = 2)) + 
    labs(x = NULL,
         y = if (freq) 
           "Count"
         else 
           "Proportion")  +
    ggtitle("whirinaki")
  
  print(graph)
}

pp_dens <- function (object, group, nsamples = 100, 
                     size = 0.25, alpha = 0.7, trim = FALSE, 
                     bw = "nrd0", adjust = 1, kernel = "gaussian", 
                     n_dens = 1024) {

  y <- object$data$repro
  group <- object$data$population
  yrep <- posterior_predict(object, nsamples = nsamples, group = group)
  data <- bayesplot::ppc_data(y, yrep, group = group)
  
  for(i in 1:5){
    p_overlay <- ggplot(data) + 
      aes_(x = ~value) + 
      stat_density(aes_(group = ~rep_id, 
                        color = "yrep"), 
                   data = function(x)
                     dplyr::filter(x, !.data$is_y),
                   geom = "line", position = "identity", 
                   size = size,
                   alpha = alpha, 
                   trim = trim,
                   bw = bw, 
                   adjust = adjust, 
                   kernel = kernel,
                   n = n_dens) + 
      stat_density(aes_(color = "y"), 
                   data = function(x) 
                     dplyr::filter(x, .data$is_y),
                   geom = "line", 
                   position = "identity",
                   lineend = "round", 
                   size = 1, trim = trim, bw = bw,
                   adjust = adjust, kernel = kernel, 
                   n = n_dens) + 
      ggforce::facet_wrap_paginate(
        vars(group),
        nrow = 2,
        ncol = 2,
        page = i,
        scales = "free_y",
      )
    
    print(p_overlay)
  }
  
  data <- filter(data, group == "Whirinaki")
  
  p_overlay <- ggplot(data) + 
    aes_(x = ~value) + 
    stat_density(aes_(group = ~rep_id, 
                      color = "yrep"), 
                 data = function(x)
                   dplyr::filter(x, !.data$is_y),
                 geom = "line", position = "identity", 
                 size = size,
                 alpha = alpha, 
                 trim = trim,
                 bw = bw, 
                 adjust = adjust, 
                 kernel = kernel,
                 n = n_dens) + 
    stat_density(aes_(color = "y"), 
                 data = function(x) 
                   dplyr::filter(x, .data$is_y),
                 geom = "line", 
                 position = "identity",
                 lineend = "round", 
                 size = 1, trim = trim, bw = bw,
                 adjust = adjust, kernel = kernel, 
                 n = n_dens) +
    ggtitle("whirinaki")
  
  print(p_overlay)
}

# The following functions are for reshaping the polygon data and converting it 
# into a data frame. The data.frame format will be made available with the paper,
# so this part is really just for transparency, and is not intended to be 
# fully reproducible.

if(!file_exists("repro_analysis/Data/demography/all_ramets_di.rds")) {
  # Utilities for carpobrotus flower model pipeline
  
  infile_data <- function(path, pop, ...) {
    
    sel_vars <- rlang::enquos(...)
    temp <- sf::st_read(path,
                        stringsAsFactors = FALSE)  %>%
      .type_poly_fields() %>%
      filter(!st_is_empty(.))
    
    # if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
    
    # Now do some initial data manipulation
    
    names(temp) <- tolower(names(temp))
    
    # calculate breaks for size bins. Only used for histogram plotting
    log_breaks <- round(dim(temp)[1]/10, digits = -1)
    
    # sf::st_area has substantially better precision than QGIS's $area macro. 
    # overwrite computed sizes from QGIS. as.numeric converts it from
    # unit-type (m^2) to a simple numeric, as some dplyr verbs don't know
    # how to handle units
    
    temp$size <- as.numeric(sf::st_area(temp))
    
    out <- temp %>%
      dplyr::mutate(population = pop) %>%
      dplyr::select(!!! sel_vars) %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(population = first(population),
                       clone_of   = first(clone_of),
                       size       = sum(size, na.rm = TRUE),
                       flower_n   = sum(flower_n, na.rm = TRUE),
                       flower_col = ifelse(!all(is.na(flower_col)),
                                           flower_col[which(!is.na(flower_col)[1])],
                                           NA_character_)) %>%
      dplyr::mutate(log_size = log(size)) %>%
      dplyr::mutate(size_bin = cut(log_size, breaks = log_breaks)) %>%
      dplyr::mutate(clean_bin = clean_bins(size_bin)) %>%
      dplyr::select(-size_bin)
    
    # add in some additional info and make sure columns are correct
    
    # without flowers, we can't be sure of the flower color 
    out$flower_col[out$flower_n == 0] <- NA_character_
    
    # QGIS automatically sets the NULLS to 0s when doing the auto-counting. Change
    # these back to NAs, as they are NOT 0s
    out$flower_n <- as.integer(out$flower_n)
    out$flower_n[out$flower_n == 0] <- NA_integer_
    out$flower_col[is.na(out$flower_n)] <- NA_character_
    
    # add in reproduction binary variable
    out$repro <- ifelse(is.na(out$flower_n), 0, 1)
    
    # check that id and clone_of are integers
    if(!is.integer(out$id)) out$id <- as.integer(out$id)
    if(!is.integer(out$clone_of)) out$clone_of <- as.integer(out$clone_of)
    
    # drop rows that are mistakenly added to attribute table
    out <- out[!is.na(out$id), ]
    
    
    return(out)
  }
  
  name_pop <- function(pop, level) {
    
    glue('{pop}_{level}')
  }
  
  clean_bins <- function(x) {
    lapply(x,
           function(x) stringr::str_split(x, '\\(|,|\\]')[[1]][3] %>%
             as.numeric()) %>% 
      unlist()
  }
  
  
  # Extras ------
  
  # Correct data types from QGIS polygons for use in C++ merge_ramets/calc_surv
  
  .type_poly_fields <- function(polygon) {
    
    names(polygon) <- tolower(names(polygon))
    
    if(!is.integer(polygon$id))       polygon$id       <- as.integer(polygon$id)
    if(!is.integer(polygon$flower_n)) polygon$flower_n <- as.integer(polygon$flower_n)
    if(!is.integer(polygon$clone_of)) polygon$clone_of <- as.integer(polygon$clone_of)
    
    return(polygon)
  }
  
  
  .ground_truth_data <- function(data, gcps, pop, time, method = 'mean_dif') {
    
    out <- switch(method,
                  'mean_dif' = .ground_truth_data_mean_dif(data, gcps, pop, time))
    
    return(out)
    
  }
  
  .ground_truth_data_mean_dif <- function(data, gcps, pop, t_n) {
    
    ground_truth_sizes_pop <- dplyr::filter(data, id >= 22000)
    
    ground_truth_sizes_tru <- dplyr::filter(gcps, pops == pop & time == t_n)
    
    # If there aren't any ground truth points in the data, then I'm not sure
    # what we can do about it. return early for now.
    
    if(nrow(ground_truth_sizes_pop) == 0) { #||
       # nrow(ground_truth_sizes_tru) == 0) {
      
      warning("No ground truth points found for population: ", pop)
      
      return(data)
      
    }
    
    for_tru                <- dplyr::right_join(ground_truth_sizes_tru,
                                                ground_truth_sizes_pop,
                                                by = 'id') %>%
      mutate(dif = (size - gcp_size) / gcp_size) %>%
      filter(!is.na(pops))
    
    # Next, probably use a distance based re-scaling, but perhaps makes
    # more sense to do it based on relative elevation of each gcp and plant?
    # I think this is where the DSMs come in handy. For now, using the mean
    # error to compute a "true" size per plant
    
    data$size     <- data$size * (1 + mean(for_tru$dif))
    data$log_size <- log(data$size)
    
    msg <- paste(
      "\nAverage error for ", pop, 
      ' at time ', t_n, 
      ' with ', length(for_tru$dif),
      ' points is: ',
      mean(for_tru$dif),
      '\n'
    )
    
    message(msg, appendLF = FALSE)
    
    return(data)
    
  }
  
  update_flower_col <- function(data, pop, id, new_col, time) {
    
    time_col <- switch(time,
                       "t_1" = 'flower_col',
                       "t_2" = 'flower_col_next')
    
    data[data$population == pop & data$id == id, time_col] <- new_col
    
    return(data)
    
  }
  
  
  source("repro_analysis/R/data_restructure/reshape_data.R")
  
}