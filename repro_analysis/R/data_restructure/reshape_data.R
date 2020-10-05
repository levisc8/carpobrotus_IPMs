
# Convert polygons to data.frames. Called from within R/01_depends.R

# This creates a data frame that does not contain spatial relationships. 
# subsequent analysis is density independent!

all_sites     <- read.csv('Data/All_Field_Sites.csv', 
                          stringsAsFactors = FALSE)

all_gr_tr     <- read.csv('Data/all_ground_truth_points.csv',
                          stringsAsFactors = FALSE)

skips <- c("Pt_Reyes", "Fort_Ord", "Omaha")

pops <- all_sites[!all_sites$Site %in% skips, 
                  c('Country', 'Site')]

all_polygons_t_1 <- list()

for(i in seq_along(pops[ ,1])) {
  
  pop <- pops$Site[i]
  
  # This will generate both ramet and genet level data frames for each population
  temp_t_1 <- infile_data(
    glue('Polygons/{pops$Country[i]}/{pop}/Polygon_Final.shp'),
    pop,
    # Variables to keep
    population, id, clone_of,
    size, flower_n, flower_col
  ) %>%
    .ground_truth_data(gcps = all_gr_tr,
                       pop = pop,
                       time = 1) %>%
    filter(id < 22000)
  
  temp_t_1$country <- all_sites$Country[all_sites$Site == pop][1]
  all_polygons_t_1 <- c(all_polygons_t_1, list(temp_t_1)) 
  names(all_polygons_t_1)[i] <- pop
  
}

no_polygons_t_1 <- map(all_polygons_t_1, 
                       .f = function(x) {
                         st_geometry(x) <- NULL
                         
                         return(x)
                       }
)

all_ramets <- do.call(rbind, no_polygons_t_1)

# Fix some inconsistency in the classification of hybrids from QGIS

all_ramets$flower_col[grepl('Y_P', 
                            all_ramets$flower_col, 
                            ignore.case = TRUE)] <- 'P_Y'


# A few plants that flowered didn't have enough visual information at 
# time T to ID flower color. Some had flower colors at T+1, so those are easy.
# For the others, I'm going to use P for St. Francis, P_Y for Springfontein, 
# Vogelgat, and Whirinaki

all_ramets <- all_ramets %>%
  select(id, population, 
         size, flower_n, flower_col,
         log_size, repro) %>%
  setNames(
    c(
      'id', 'population',
      'size', 'flower_n', 'flower_col',
      'log_size', 'repro'
    )
  )


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
                    'P_Y',
                    't_1') %>%
  update_flower_col('Springfontein',
                    44,
                    'P_Y',
                    't_1') %>%
  update_flower_col('St_Francis',
                    110,
                    'P',
                    't_1') %>%
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
                    't_1')

all_ramets <- mutate(all_ramets, id_pop = paste(id, population, sep = "_"))


saveRDS(all_ramets, file = 'repro_analysis/Data/demography/all_ramets_di.rds')
