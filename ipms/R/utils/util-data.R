# data reading/writing utilities

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
  
  if(sf::st_crs(temp)$epsg != 4326) temp <- sf::st_transform(temp, 4326)
  
  # Substitute merged ramet IDs so values can be recomputed 
  # C++ source is in R/Cpp/cpp_utils.cpp
  merged_ramets <- dplyr::filter(merged_ramets, population == pop)
  
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


.type_poly_fields <- function(polygon) {
  
  names(polygon) <- tolower(names(polygon))
  
  if(!is.integer(polygon$id))       polygon$id       <- as.integer(polygon$id)
  if(!is.integer(polygon$flower_n)) polygon$flower_n <- as.integer(polygon$flower_n)
  if(!is.integer(polygon$clone_of)) polygon$clone_of <- as.integer(polygon$clone_of)
  
  return(polygon)
}


ground_truth_data <- function(data, gcps, pop, time, method = 'mean_dif') {
  
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
