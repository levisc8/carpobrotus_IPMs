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


# Utilities for carpobrotus IPM pipeline

infile_data <- function(path, pop) {
  
  temp <- st_read(paste(getwd(), path, sep = '/'),
                  stringsAsFactors = FALSE)  
  
  names(temp) <- tolower(names(temp))
  
  # calculate breaks for size bins. Only used for histogram plotting
  log_breaks <- round(dim(temp)[1]/10, digits = -1)
  
  # We only have precision to 7 0s for size in Qgis, and some ramets are smaller
  # Trying sf::st_area to see if it has higher precision
  temp$size[temp$size == 0] <- sf::st_area(temp[temp$size == 0, ])
  
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
  
  # Q Automatically sets the NULLS to 0s when doing the auto-counting. Change
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
  group_var <- quos(...)
  
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

