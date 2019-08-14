## Prototype survival and growth for Israel 

# First, figure out where the most recent polygons are. I do the digitization on
# my laptop, but these files should really live on the server so we can all access
# them. The next few lines make sure that the rest of the script has the most
# recent set of polygons from Israel

local_path <- 'Ana_Israel_IPM/Data/Havatselet/'

polys <- paste(local_path, c("Polygon_2019_Final",
                             'Polygon_Final'), 
               '.shp', 
               sep = "")

# This file tracks which ramets get merged together as a result of overlapping
# growth between T and T+1

merged_ramets <- read.csv("Ana_Israel_IPM/Data/merged_ramets.csv",
                          stringsAsFactors = FALSE)

rm_ramets <- read.csv("Ana_Israel_IPM/Data/t_2_omissions.csv",
                      stringsAsFactors = FALSE)

# _local deals with t_2 because the data are structured a little differently
# and require different manipulations. These functions will eventually get 
# refactored to work on a switch() call internally so that you only need to call
# infile_data(), but I'm too lazy to deal with that right now ;)

is_t_2 <- infile_data(polys[1], 'Havatselet',
                      type = 'local_t_2',
                      # Variables to keep
                      population, id, clone_of,
                      size, flower_n, flower_col, alive)

is_t_1 <- infile_data(polys[2], 'Havatselet',
                      type = 'local_t_1',
                      # Variables to keep
                      population, id, clone_of,
                      size, flower_n, flower_col,
                      merged_ramets = merged_ramets) %>%
  filter(!id %in% rm_ramets$ramet) 

# Join data together for modeling

data_t_2 <- setNames(is_t_2,
                     c(
                       "id", 'population', 'clone_of', 'size_next', 'flower_n_next', 
                       'flower_col', 'survival',
                       'log_size_next', 'clean_bin_next', 'geometry',
                       'repro_next'
                     )
) %>%
  as_tibble() %>%
  select(-geometry)

data_t_1 <- setNames(is_t_1,
                     c(
                       'id', 'population', 'clone_of',
                       'size', 'flower_n', 'flower_col',
                       'log_size', 'clean_bin',
                       'geometry', 'repro'
                     )
) %>%
  as_tibble() %>%
  select(-geometry)

all_data <- full_join(data_t_1, data_t_2, by = 'id')

# update_survival inserts 0s for plants that died and 1s for everything else.
# For plants that are new this year, a 1 is incorrect - these should be NA.

all_data$survival <- update_survival(all_data$population.y, 
                                     all_data$survival)

all_data$survival[is.na(all_data$log_size)] <- NA_integer_
