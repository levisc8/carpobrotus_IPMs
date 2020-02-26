# Dependencies
library(sf) # Handle spatial data
library(fs) # file manipulation
library(stringr) # string manipulation
library(ggplot2) # plotting 
library(dplyr) # data manipulation
library(gridExtra) # plotting
library(glue) # string manipulation
library(purrr)
library(lme4)
library(brms) # attempt to get fecundity models converging
library(mgcv)
library(Rcpp)


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
    .type_poly_fields()
  
  if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
  # Substitute merged ramet IDs so values can be recomputed 
  # C++ source is in R/Cpp/cpp_utils.cpp
  merged_ramets <- dplyr::filter(merged_ramets, population == pop)
  
  temp$id <- merge_ramets(temp$id, 
                          merged_ramets$absorbed_ramet,
                          merged_ramets$big_ramet)
  
  # Now do some initial data manipulation
  
  names(temp) <- tolower(names(temp))
  
  # calculate breaks for size bins. Only used for histogram plotting
  log_breaks <- round(dim(temp)[1]/10, digits = -1)
  
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


.infile_data_local_t_2 <- function(path, pop, ...) {
  
  sel_vars <- rlang::enquos(...)
  temp <- sf::st_read(path,
                  stringsAsFactors = FALSE) %>%
    .type_poly_fields() 
  
  if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
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
    dplyr::mutate(log_size = log(size)) %>%
    dplyr::mutate(size_bin = cut(log_size, breaks = log_breaks)) %>%
    dplyr::mutate(clean_bin = clean_bins(size_bin)) %>%
    dplyr::select(-size_bin)
  
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

sourceCpp(file = 'R/Cpp/cpp_utils.cpp')

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
    ' with ', ength(for_tru$dif),
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

# Initialize output directories

if(!dir_exists('Model_Fits')) dir_create('Model_Fits')
if(!dir_exists('Model_Summaries')) dir_create('Model_Summaries')
if(!dir_exists('Model_Summaries/Trace_Plots')) dir_create('Model_Summaries/Plots')
if(!dir_exists('Stan')) dir_create('Stan')
if(!dir_exists('Figures')) dir_create('Figures')

