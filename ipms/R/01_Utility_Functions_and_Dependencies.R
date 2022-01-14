# Dependencies
library(sf)        # Handle spatial data
# library(s2)
library(fs)        # file manipulation
library(stringr)   # string manipulation
library(ggplot2)   # plotting
library(dplyr)     # data manipulation
library(gridExtra) # plotting
library(tidyr)     # Climate data munging
library(glue)      # string manipulation
library(purrr)
# library(lme4)
library(brms)      # attempt to get fecundity models converging
library(mgcv)
library(Rcpp)

# loading ggplot2 and gridExtra here causes RStudio to crash for some reason.
# these are now loaded in 04_Data_Exploration.R, which for some reason
# causes no such crashes(?????????)


# Utilities for carpobrotus IPM pipeline

infile_data <- function(path, pop, 
                        type = c('remote', 'local_t_1', 'local_t_2'),
                        ...,
                        merged_ramets = NULL) {
  
  dots <- rlang::enquos(...)
  
  out <- switch(type,
                'remote'    = .infile_remote(path, pop),
                'local_t_1' = .infile_data_local_t_1(path, 
                                                     pop, 
                                                     !!! dots, 
                                                     merged_ramets = merged_ramets),
                'local_t_2' = .infile_data_local_t_2(path, 
                                                     pop, 
                                                     !!! dots)
  )
  return(out)
}



# Without paste(getwd()) to operate on locally stored polygons
.infile_data_local_t_1 <- function(path, pop, ..., merged_ramets) {
  
  sel_vars <- rlang::enquos(...)
  temp <- sf::st_read(path,
                      stringsAsFactors = FALSE)  %>%
    .type_poly_fields() %>%
    st_make_valid()
  
  # if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
  # Substitute merged ramet IDs so values can be recomputed 
  # C++ source is in R/Cpp/cpp_utils.cpp
  merged_ramets <- dplyr::filter(merged_ramets, population == pop)
  
  # temp$id <- merge_ramets(temp$id, 
  #                         merged_ramets$absorbed_ramet,
  #                         merged_ramets$big_ramet)
  
  for(i in merged_ramets$absorbed_ramet) {
    
    temp$id[temp$id == i] <- merged_ramets[merged_ramets$absorbed_ramet == i, "big_ramet"]
    
  }
  
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
    dplyr::mutate(log_size = log(size))
  
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


.infile_data_local_t_2 <- function(path, pop, ...) {
  
  sel_vars <- rlang::enquos(...)
  temp <- sf::st_read(path,
                  stringsAsFactors = FALSE) %>%
    .type_poly_fields() %>%
    st_make_valid()
  
  # if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
  # Now do some initial data manipulation
  
  names(temp) <- tolower(names(temp))
  
  if('survival' %in% names(temp)) {
    names(temp) <- gsub('survival', 'alive', names(temp))
    names(temp) <- gsub('surv', 'alive', names(temp))
  }
  
  # calculate breaks for size bins. Only used for histogram plotting
  log_breaks <- round(dim(temp)[1]/10, digits = 0)
  
  # sf::st_area has substantially better precision than QGIS's $area macro. 
  # to overwrite computed sizes from QGIS. as.numeric converts it from
  # unit-type (m^2) to a simple numeric, as some dplyr verbs don't know
  # how to handle units
  
  temp$size <- as.numeric(sf::st_area(temp))
  
  out <- temp %>%
    dplyr::mutate(population = pop) %>%
    dplyr::select(!!! sel_vars) %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(population = first(population),
              clone_of   = first(clone_of),
              size       = sum(size),
              flower_n   = sum(flower_n, na.rm = TRUE),
              flower_col = ifelse(!all(is.na(flower_col)),
                                  flower_col[which(!is.na(flower_col)[1])],
                                  NA_character_),
              alive = mean(alive, na.rm = TRUE)) %>%
    dplyr::mutate(log_size = log(size)) 
  
  # add in some additional info and make sure columns are correct
  
  # without flowers, we can't be sure of the flower color 
  out$flower_col[out$flower_n == 0] <- NA_character_
  
  # QGIS automatically sets the NULLS to 0s when doing the auto-counting. Change
  # these back to NAs, as they are NOT 0s
  
  out$flower_n                        <- as.integer(out$flower_n)
  out$flower_n[out$flower_n == 0]     <- NA_integer_
  out$flower_col[is.na(out$flower_n)] <- NA_character_
  out$alive                           <- as.integer(out$alive)
  
  # add in reproduction binary variable
  
  out$repro <- ifelse(is.na(out$flower_n), 0, 1)
  
  # drop rows that are mistakenly added to attribute table
  out <- out[!is.na(out$id), ]
  
  
  return(out)
}

.infile_remote <- function(path, pop) {
  
  temp <- sf::st_read(paste(getwd(), path, sep = '/'),
                  stringsAsFactors = FALSE)  
  
  names(temp) <- tolower(names(temp)) 
  temp        <- .type_poly_fileds(temp)
  
  if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
  # calculate breaks for size bins. Only used for histogram plotting
  log_breaks <- round(dim(temp)[1]/10, digits = -1)
  
  # We only have precision to 7 0s for size in Qgis, and some ramets are smaller
  # Trying sf::st_area to see if it has higher precision
  temp$size <- as.numeric(sf::st_area(temp))
  
  out <- temp %>%
    dplyr::mutate(population = pop) %>%
    dplyr::select(population, id, clone_of,
                  size, flower_n, flower_col) %>%
    dplyr::mutate(log_size = log(size)) %>%
    dplyr::mutate(size_bin = cut(.$log_size, breaks = log_breaks)) %>%
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

ramet_to_genet <- function(data, ...) {
  group_var <- rlang::enquos(...)
  
  out <- data %>%
    dplyr::group_by(!!! group_var) %>%
    dplyr::summarise(n_ramets = length(id),
              country = first(country),
              size = sum(size),
              flower_n = sum(flower_n, na.rm = TRUE),
              flower_col = first(flower_col)) %>%
    dplyr::mutate(log_size = log(size)) 
  
  out$flower_n[out$flower_n == 0] <- NA_integer_
  
  out$repro <- ifelse(is.na(out$flower_n), 0, 1)
  
  return(out)
  
}

clean_bins <- function(x) {
  lapply(x,
         function(x) stringr::str_split(x, '\\(|,|\\]')[[1]][3] %>%
           as.numeric()) %>% 
    unlist()
}


get_gg_legend<-function(plot){
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


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

sourceCpp(file = 'ipms/R/Cpp/cpp_utils.cpp')

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
  
  if(dim(ground_truth_sizes_pop)[1] == 0) {
    
    return(data)
    
  }
  
  for_tru                <- dplyr::right_join(ground_truth_sizes_tru,
                                              ground_truth_sizes_pop,
                                              by = 'id') %>%
    mutate(dif = (size - gcp_size) / gcp_size) 
  
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

check_outliers <- function(all_ramets, population, type) {
  
  switch(type,
         "growth"  = check_growth_outliers(all_ramets, population),
         "big_new" = check_big_new_outerliers(all_ramets, population))
  
}

fit_vr_model <- function(data, vr, clim) {
  
  fam <- switch(
    vr,
    "repro"     = bernoulli(),
    "flower_n"  = negbinomial(),
    "alive"     = bernoulli(),
    "size_next" = gaussian()
  )
  
  form_1 <- as.formula(glue("{vr} ~ log_size + (1|site)"))
  form_2 <- as.formula(glue("{vr} ~ log_size + {clim} + (1|site)"))
  form_3 <- as.formula(glue("{vr} ~ log_size * {clim} + (1|site)"))
  
  be <- "cmdstanr"
  
  if(vr %in% c("flower_n")) {
    inits <- "0"
    
  } else {
    
    inits <- "random"
  }
  
  if(vr == "alive") vr <- "survival"
  
  # NB: Sigma has log link by default! Remember to exponentiate at IPM
  # build time.
  if(vr == "log_size_next") {
    
    sigma_form_1 <- sigma ~ log_size
    sigma_form_2 <- as.formula(glue("sigma ~ log_size + {clim}"))
    sigma_form_3 <- as.formula(glue("sigma ~ log_size * {clim}"))
    
    
    form_1 <- bf(form_1, sigma_form_1)
    form_2 <- bf(form_2, sigma_form_2)
    form_3 <- bf(form_3, sigma_form_3)
    
    vr <- "growth"
  }
  
  
  message("Model 1: size only----------\n\n")
  mod_1 <- brm(form_1,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_int_only.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 2: Size Plus climate-------\n\n")
  
  mod_2 <- brm(form_2,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_size_plus_clim.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 3: Size Times Climate------\n\n")
  
  mod_3 <- brm(form_3,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_size_times_clim.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))

  mod_waic <- waic(mod_1, mod_2, mod_3)
  
  out <- list(size_only       = mod_1,
              size_plus_clim  = mod_2,
              size_times_clim = mod_3,
              waic            = mod_waic)
  
  return(out)

  
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
    
    pdf(glue("ipms/Model_Summaries/Trace_Plots/{vr}/{file_nm}_no_is.pdf"))
    
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
    
    
    sink(file = glue('ipms/Model_Summaries/{vr}/{file_nm}_no_is.txt')) 
    cat('Size only\n\n *********************\n\n')
    print(summary(use_mod[[1]]))
    cat('\n\n*********************\n\nSize + Climate',
        '\n\n*********************\n\n')
    print(summary(use_mod[[2]]))
    cat('\n\n*********************\n\nSize * Climate',
        '\n\n*********************\n\n')
    print(summary(use_mod[[3]]))
    cat('\n\n*********************\n\nWAIC Results\n\n')
    print(use_mod[[4]])
    
    cat('\n\nEnd output')
    sink()
    
  }
  
}

plot_preds <- function(models, vr) {
  
  z <- c(models[[1]]$size_only$data$log_size,
         models[[1]]$size_only$data$log_size_next)
  site <- unique(models[[1]]$size_only$data$site)
  
  min_z <- min(z, na.rm = TRUE)
  max_z <- max(z, na.rm = TRUE)
  
  pred_data <- data.frame(
    log_size     = seq(min_z, max_z, length.out = 100),
    mat_rec      = seq(-1.8, 2.3, length.out = 100),
    map_rec      = seq(-1.32, 1, length.out = 100),
    t_seas_rec   = seq(-1.45, 0.8, length.out = 100),
    p_seas_rec   = seq(-1.9, 2.1, length.out = 100),
    t_co_qu_rec  = seq(-1.04, 1.35, length.out = 100)
  ) 
  
  all_data <- list()
  
  temp <- list()
  
  for(i in seq_along(site)) {
    all_data[[i]] <- cbind(pred_data, site = site[i])
  }
  
  pred_data <- do.call(rbind, all_data)
  
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
  
  all_data <- do.call(rbind, temp) %>%
    filter(site %in% c("Melkboss", "Rough_Island", "Foxton", "Praia_de_Areao"))
  
  plt <- ggplot(all_data, 
                aes(x = log_size,
                    y = Estimate)) + 
    geom_line(aes(linetype = mod_type, 
                  color = site),
              size = 1.2,
              alpha = 0.6) + 
    facet_wrap(~clim, 
               scales = "free") +
    theme_bw() +
    scale_color_discrete(guide = "none")
    
  pdf(glue("ipms/Figures/vr_models/{vr}_model_predictions_no_is.pdf"))
    print(plt)
  dev.off()
  
}

# Initialize output directories

if(!dir_exists('ipms/Model_Fits'))                  dir_create('ipms/Model_Fits')
if(!dir_exists('ipms/Model_Summaries'))             dir_create('ipms/Model_Summaries')
if(!dir_exists('ipms/Model_Summaries/Trace_Plots')) dir_create('ipms/Model_Summaries/Plots')
if(!dir_exists('ipms/Stan'))                        dir_create('ipms/Stan')
if(!dir_exists('ipms/Figures'))                     dir_create('ipms/Figures')

