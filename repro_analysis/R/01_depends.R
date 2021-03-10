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

hs_diagnostic_plots <- function(clim, func) {
  
  mod_path <- glue("repro_analysis/model_fits/shrink_mods/model_objects_{clim}_{func}.rds")
  
  vr_mod_list <- readRDS(mod_path)
  
  repro_mod <- vr_mod_list$repro_mod
  flower_n_mod <- vr_mod_list$flower_n_mod
  
  pdf(glue("repro_analysis/Manuscript/Figures/shrink_mod_{clim}_{func}.pdf"))
  
  plot(repro_mod,
       ask = FALSE)
  
  
  p <- pp_check(repro_mod,
                type     = 'bars',
                nsamples = 100L,
                freq     = FALSE)
  
  print(p)
  
  p <- pp_check(repro_mod,
                type     = 'dens_overlay',
                nsamples = 100L)
  
  print(p)
  
  
  
  
  plot(flower_n_mod,
       ask = FALSE)
  
  
  p <- pp_check(flower_n_mod,
                type     = 'bars',
                nsamples = 100L,
                freq     = FALSE) + 
    xlim(c(0, 20))
  
  print(p)
  
  p <- pp_check(flower_n_mod,
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