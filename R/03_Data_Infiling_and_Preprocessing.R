# Basic shape file analysis

# This creates a two batches of data - one that ignores spatial relationships
# (all_ramets, i.e. drops geometry column) and 2 that preserve them (all_polygons_t_1
# and all_polygons_t_2). infile_data automatically converts all polygons
# that aren't in EPSG 4326 to that coordinate system, so these could eventually
# be used to make models explicitly including conspecific density (more on that 
# in 04_Fecundity_Models)

all_sites     <- read.csv('Data/All_Field_Sites.csv', 
                          stringsAsFactors = FALSE)

all_gr_tr     <- read.csv('Data/all_ground_truth_points.csv',
                          stringsAsFactors = FALSE)

merged_ramets <- read.csv('Data/merged_ramets.csv',
                          stringsAsFactors = FALSE)

rm_ramets     <- read.csv('Data/t_2_omissions.csv',
                          stringsAsFactors = FALSE)

pops <- all_sites[all_sites$Demography == 1 &
                    all_sites$Polygon_checked == 1 &
                    all_sites$Polygon_2019_checked == 1 &
                    all_sites$Notes == 'ready', 
                  c('Country', 'Site')]

all_polygons_t_1 <- all_polygons_t_2 <- list()

for(i in seq_along(pops[ ,1])) {
  
  pop <- pops$Site[i]
  
  # This will generate both ramet and genet level data frames for each population
  temp_t_1 <- infile_data(
    glue('Polygons/{pops$Country[i]}/{pop}/Polygon_Final.shp'),
    pop,
    type = 'local_t_1',
    # Variables to keep
    population, id, clone_of,
    size, flower_n, flower_col,
    merged_ramets = merged_ramets
  ) %>%
    .ground_truth_data(gcps = all_gr_tr,
                       pop = pop,
                       time = 1) %>%
    filter(id < 22000)
  
  temp_t_2 <- infile_data(
    glue('Polygons/{pops$Country[i]}/{pop}/Polygon_2019_Final.shp'),
    pop,
    type = 'local_t_2',
    population, id, clone_of,
    size, flower_n, flower_col, alive
   ) %>%
     .ground_truth_data(gcps = all_gr_tr,
                        pop = pop,
                        time = 2) %>%
    filter(id < 22000)

  temp_t_1$country <- all_sites$Country[all_sites$Site == pop][1]
  all_polygons_t_1 <- splice(all_polygons_t_1, list(temp_t_1)) 
  names(all_polygons_t_1)[i] <- pop
  
  temp_t_2$country <- all_sites$Country[all_sites$Site == pop][1]
  all_polygons_t_2 <- splice(all_polygons_t_2, list(temp_t_2)) 
  names(all_polygons_t_2)[i] <- pop
  
}

no_polygons_t_1 <- map(all_polygons_t_1, 
                   .f = function(x) {
                     st_geometry(x) <- NULL
                     
                     return(x)
                   }
)

all_ramets_t_1 <- do.call(rbind, no_polygons_t_1)

no_polygons_t_2 <- map(all_polygons_t_2, 
                       .f = function(x) {
                         st_geometry(x) <- NULL
                         
                         return(x)
                       }
)

all_ramets_t_2 <- do.call(rbind, no_polygons_t_2)

all_ramets <- full_join(all_ramets_t_1,
                        all_ramets_t_2, 
                        by = c('population', 'id'))


# Fix some inconsistency in the classification of hybrids from QGIS

all_ramets$flower_col.x[grepl('Y_P', 
                              all_ramets$flower_col.x, 
                              ignore.case = TRUE)] <- 'P_Y'
all_ramets$flower_col.y[grepl('Y_P', 
                              all_ramets$flower_col.y, 
                              ignore.case = TRUE)] <- 'P_Y'

# A few plants that flowered didn't have enough visual information at 
# time T to ID flower color. Some had flower colors at T+1, so those are easy.
# For the others, I'm going to use P for St. Francis, P_Y for Springfontein, 
# Vogelgat, and Whirinaki

all_ramets <- all_ramets %>%
  select(id, population, clone_of.y, size.x, flower_n.x, flower_col.x,
         log_size.x, repro.x, size.y, flower_n.y, flower_col.y, alive,
         log_size.y, repro.y, country.y) %>%
  setNames(
    c(
      'id', 'population', 'clone_of',
      'size', 'flower_n', 'flower_col',
      'log_size', 'repro', 'size_next', 
      'flower_n_next', 'flower_col_next',
      'alive', 'log_size_next', 'repro_next',
      'country'
    )
  )
  

# If they had a size at t_1 but do not at t_2 and they weren't recorded in 
# t_2_omissions, then we are assuming they are dead. New recruits won't have
# a size at t_1, but have a size at t_2, and so they shouldn't be included in 
# in this

all_ramets$alive[!is.na(all_ramets$log_size) & 
                   is.na(all_ramets$log_size_next)] <- 0

# A few plants that flowered didn't have enough visual information at 
# time T to ID flower color. Some had flower colors at T+1, so those are easy.
# For the others, I'm going to use P for St. Francis, P_Y for Springfontein, 
# Vogelgat, and Whirinaki

all_ramets <- update_flower_col(all_ramets, 
                                "Struisbaai",
                                66,
                                "P",
                                't_1') %>%
  update_flower_col('Vogelgat',
                    1,
                    'Y_P',
                    't_1') %>%
  update_flower_col('Vogelgat',
                    1,
                    'Y_P',
                    't_2') %>%
  update_flower_col('Springfontein',
                    44,
                    'P_Y',
                    't_1') %>%
  update_flower_col('St_Francis',
                    110,
                    'P',
                    't_1') %>%
  update_flower_col('St_Francis',
                    110,
                    'P',
                    't_2') %>%
  update_flower_col('Rough_Island',
                    54,
                    'P_Y',
                    't_1') %>%
  update_flower_col('Rough_Island',
                    262,
                    'P_Y',
                    't_1') %>%
  update_flower_col('Rough_Island',
                    290,
                    'P',
                    't_1') %>%
  update_flower_col('Whirinaki',
                    142,
                    'P',
                    't_1') %>%
  update_flower_col('Whirinaki',
                    147,
                    'P',
                    't_1') %>%
  update_flower_col('Struisbaai',
                    72,
                    'P',
                    't_2')
  
# Store data and sned it to the RStudio server project folder.
saveRDS(all_ramets, file = 'Data/all_ramets_di.rds')

dest <- 'I:/sie/102_Data_SL/carpobrotus_ipms/Data/'

file_copy('Data/all_ramets_di.rds', dest, overwrite = TRUE)
