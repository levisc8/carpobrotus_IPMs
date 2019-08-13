# Dependencies
library(sf) # Handle spatial data
library(fs) # file manipulation
library(stringr) # string manipulation
library(ggplot2) # plotting 
library(dplyr) # data manipulation
library(gridExtra) # plotting
library(glue) # string manipulation
library(tidyr)
library(rlang)
library(purrr)
library(mgcv)
library(Rcpp)
library(ipmr)
library(nlme)


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
  temp <- st_read(path,
                  stringsAsFactors = FALSE)  
  
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
    mutate(population = pop) %>%
    dplyr::select(!!! sel_vars) %>%
    group_by(id) %>%
    summarise(population = first(population),
              clone_of   = first(clone_of),
              size       = sum(size),
              flower_n   = sum(flower_n, na.rm = TRUE),
              flower_col = ifelse(!all(is.na(flower_col)),
                                  flower_col[which(!is.na(flower_col)[1])],
                                  NA_character_)) %>%
    mutate(log_size = log(size)) %>%
    mutate(size_bin = cut(log_size, breaks = log_breaks)) %>%
    mutate(clean_bin = clean_bins(size_bin)) %>%
    select(-size_bin)
  
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
  temp <- st_read(path,
                  stringsAsFactors = FALSE)  
  
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
    mutate(population = pop) %>%
    dplyr::select(!!! sel_vars) %>%
    group_by(id) %>%
    summarise(population = first(population),
              clone_of   = first(clone_of),
              size       = sum(size),
              flower_n   = sum(flower_n, na.rm = TRUE),
              flower_col = ifelse(!all(is.na(flower_col)),
                                  flower_col[which(!is.na(flower_col)[1])],
                                  NA_character_),
              alive = mean(alive, na.rm = TRUE)) %>%
    mutate(log_size = log(size)) %>%
    mutate(size_bin = cut(log_size, breaks = log_breaks)) %>%
    mutate(clean_bin = clean_bins(size_bin)) %>%
    select(-size_bin)
  
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

.infile_remote <- function(path, pop) {
  
  temp <- st_read(paste(getwd(), path, sep = '/'),
                  stringsAsFactors = FALSE)  
  
  names(temp) <- tolower(names(temp))
  
  # calculate breaks for size bins. Only used for histogram plotting
  log_breaks <- round(dim(temp)[1]/10, digits = -1)
  
  # We only have precision to 7 0s for size in Qgis, and some ramets are smaller
  # Trying sf::st_area to see if it has higher precision
  temp$size <- as.numeric(sf::st_area(temp))
  
  out <- temp %>%
    mutate(population = pop) %>%
    dplyr::select(population, id, clone_of,
                  size, flower_n, flower_col) %>%
    mutate(log_size = log(size)) %>%
    mutate(size_bin = cut(.$log_size, breaks = log_breaks)) %>%
    mutate(clean_bin = clean_bins(size_bin)) %>%
    select(-size_bin)
  
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
  group_var <- enquos(...)
  
  out <- data %>%
    group_by(!!! group_var) %>%
    summarise(n_ramets = length(id),
              country = first(country),
              size = sum(size),
              flower_n = sum(flower_n, na.rm = TRUE),
              flower_col = first(flower_col)) %>%
    mutate(log_size = log(size)) 
  
  out$flower_n[out$flower_n == 0] <- NA_integer_
  
  out$repro <- ifelse(is.na(out$flower_n), 0, 1)
  
  return(out)
  
}

clean_bins <- function(x) {
  lapply(x,
         function(x) str_split(x, '\\(|,|\\]')[[1]][3] %>%
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

sourceCpp(file = 'Ana_Israel_IPM/Cpp/cpp_utils.cpp')

# Top level wrappers for sensitivity and elasticity computations. Pseudo-generic,
# but dispatches on arguments rather than classes

sensitivity <- function(K, h, level = c('kernel', 'vr', 'param'),
                        ...) {
  
  switch(level,
         "kernel" = .kernel_sens(K, h),
         "vr"     = .vr_sens(K, h, ...))
}

elasticity <- function(K, h, level = c('kernel', 'vr', 'param'),
                       ...) {
  
  switch(level,
         "kernel" = .kernel_elas(K, h),
         "vr"     = .vr_elas(K, h, ...))
}


# Sensitivity methods for kernel, vital rates, and parameters
.kernel_sens <- function(K, h) {
  w <- Re(eigen(K)$vectors[ , 1])
  v <- Re(eigen(t(K))$vectors[ , 1])
  
  out <- outer(v, w) / sum(v * w * h)
  
  return(out)
}


# Elasticity methods for kernel, vital rates, and parameters

.kernel_elas <- function(K, h) {
  
  sens <- sensitivity(K, h, 'kernel')
  
  lambda <- Re(eigen(K)$values[1])
  
  out <- sens * (K / h) / lambda
  
  return(out)
}


mat_to_df <- function(mat, meshp) {
  
  out <- expand.grid(list(x = meshp, y = meshp)) %>%
    mutate(value = NA_real_)
  
  it <- 1
  
  for(i in seq_len(dim(mat)[1])) {
    for(j in seq_len(dim(mat)[2])) {
      
      out[it, 3] <- mat[i, j]
      it <- it + 1
      
    }
  }
  
  return(out)
}
