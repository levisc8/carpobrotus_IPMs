# Basic shape file analysis. For some reason, this script crashes RStudio, but
# runs fine from the R command line. See debugging/data-preprocessing-w-command-line.R.

# This creates a two batches of data - one that ignores spatial relationships
# (all_ramets, i.e. drops geometry column) and 2 that preserve them (all_polygons_t_1
# and all_polygons_t_2). infile_data automatically converts all polygons
# that aren't in EPSG 4326 to that coordinate system, so these could eventually
# be used to make models explicitly including conspecific density (more on that 
# in 04_Fecundity_Models)

all_sites     <- read.csv('ipms/Data/All_Field_Sites.csv', 
                          stringsAsFactors = FALSE)

all_gr_tr     <- read.csv('ipms/Data/all_ground_truth_points.csv',
                          stringsAsFactors = FALSE)

merged_ramets <- read.csv('ipms/Data/merged_ramets.csv',
                          stringsAsFactors = FALSE)

rm_ramets     <- read.csv('ipms/Data/t_2_omissions.csv',
                          stringsAsFactors = FALSE)

buried_ramets <- read.csv('ipms/Data/buried_ramets.csv',
                          stringsAsFactors = FALSE) %>%
  mutate(pop_id = paste(population, id, sep = "_"))

pops <- all_sites[all_sites$Demography == 1 &
                    all_sites$Polygon_checked == 1 &
                    all_sites$Polygon_2019_checked == 1 &
                    # all_sites$Site != "Melkboss" &
                    all_sites$Notes == 'ready', # &
                    #all_sites$Country != 'Israel', 
                  c('Country', 'Site')] %>%
  arrange(Site)

out <- list()

for(i in pops$Site) {
  
  pop <- i
  country <- pops$Country[pops$Site == i]
  
  debug_log <- fs::file_create(glue("ipms/debugging/{pop}.txt"))
  
  sink(file = ,
       type = c("output", "message"))
  
  # This will generate both ramet and genet level data frames for each population
  
  temp <- infile_data(
    glue('ipms/Polygons/{country}/{pop}/Polygon_Final.shp'),
    pop,
    type = 'local_t_1',
    # Variables to keep
    population, id, clone_of,
    size, flower_n, flower_col,
    merged_ramets = merged_ramets
  ) %>%
    .ground_truth_data(
      gcps = all_gr_tr,
      pop = pop,
      time = 1) %>%
    filter(id < 22000)
  
  st_geometry(temp) <- NULL
  
  temp_t_2 <- infile_data(
    glue('ipms/Polygons/{country}/{pop}/Polygon_2019_Final.shp'),
    pop,
    type = 'local_t_2',
    population, id, clone_of,
    size, flower_n, flower_col, alive
  ) %>%
    .ground_truth_data(gcps = all_gr_tr,
                       pop = pop,
                       time = 2) %>%
    filter(id < 22000)
  # st_geometry(temp_t_2) <- NULL
  
  temp$country <- all_sites$Country[all_sites$Site == pop][1]
  temp_t_2$country <- all_sites$Country[all_sites$Site == pop][1]
  

  out <- c(
    out, 
    list(full_join(temp,
                   temp_t_2,
                   by = c('population', 'id'))
    )
  )

  rm(temp, temp_t_2)
  
  sink()
}

out <- setNames(out, pops$Site)


# no_polygons <- map(  
#   pops$Site,
#   function(pop, all_sites, all_gr_tr, merged_ramets) {
#   sink(file = glue("ipms/debugging/{pop}_t_1.txt"),
#        type = c("output", "message"))
#     
#     on.exit(sink())
#     
#    country <- all_sites$Country[all_sites$Site == pop]
# 
#   # This will generate both ramet and genet level data frames for each population
#    temp <- infile_data(
#      glue('ipms/Polygons/{country}/{pop}/Polygon_Final.shp'),
#      pop,
#      type = 'local_t_1',
#      # Variables to keep
#      population, id, clone_of,
#      size, flower_n, flower_col,
#      merged_ramets = merged_ramets
#    ) %>%
#      .ground_truth_data(
#        gcps = all_gr_tr,
#        pop = pop,
#        time = 1) %>%
#      filter(id < 22000)
# 
#    st_geometry(temp) <- NULL
#    
#    temp_t_2 <- infile_data(
#      glue('ipms/Polygons/{country}/{pop}/Polygon_2019_Final.shp'),
#      pop,
#      type = 'local_t_2',
#      population, id, clone_of,
#      size, flower_n, flower_col, alive
#    ) %>%
#      .ground_truth_data(gcps = all_gr_tr,
#                         pop = pop,
#                         time = 2) %>%
#      filter(id < 22000)
#    # st_geometry(temp_t_2) <- NULL
#    
#    temp$country <- all_sites$Country[all_sites$Site == pop][1]
#    temp_t_2$country <- all_sites$Country[all_sites$Site == pop][1]
#    
#    out <- full_join(temp,
#                     temp_t_2,
#                     by = c('population', 'id'))
#    
#    rm(temp, temp_t_2)
#    
#   return(out)
#   },
#   all_sites = all_sites,
#   all_gr_tr = all_gr_tr, 
#   merged_ramets = merged_ramets
# ) %>%
#   setNames(pops$Site)


all_ramets <- do.call(rbind, out)

# Correct survival for clones/new recruits and a data entry error

all_ramets$alive[is.na(all_ramets$size)] <- NA_integer_
all_ramets$alive[all_ramets$alive == 2]  <- 1L

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
                    'P_Y',
                    't_1') %>%
  update_flower_col('Vogelgat',
                    1,
                    'P_Y',
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
  
# For those that flowered twice but have differing flower colors,
# use T+1 designations. I trust those a bit more.

all_ramets <- mutate(all_ramets, id_pop = paste(id, population, sep = "_"))

fl_col_diffs <- filter(all_ramets, 
                       !is.na(flower_col) & !is.na(flower_col_next)) %>%
  filter(flower_col != flower_col_next)

fl_col_diffs$flower_col <- fl_col_diffs$flower_col_next


# Now, replace the old ones w/ the new ones

all_ramets <- filter(all_ramets,
                     ! id_pop %in% fl_col_diffs$id_pop) %>% # remove one's we updated 
  bind_rows(fl_col_diffs) %>%                               # add them back in
  select(-id_pop)


# Some had flowers at both times, but flower color at T was unidentifiable because
# the petals had already dropped. In cases where we know flower color at T+1, 
# we can use that to infer flower color at T


all_ramets$flower_col[is.na(all_ramets$flower_col) &
                        !is.na(all_ramets$flower_col_next) &
                        !is.na(all_ramets$flower_n)] <- 
  all_ramets$flower_col_next[is.na(all_ramets$flower_col) & 
                               !is.na(all_ramets$flower_col_next) &
                               !is.na(all_ramets$flower_n)]

# Buried plants ------------ 
# A number of plants were buried/covered by trees at T and exposed at T+1 (or vice versa), 
# giving them unrealistic growth increments. We can't really use those for growth, as they'll
# skew our estimates of the mean and variance in growth as a function of size. 
# However, we can use them for survival and fecundity because those only rely on
# estimates of size at time T! Therefore, going to create a separate file for 
# growth data that will get used in those models, versus a general file with ALL
# THE RAMETS for use in other models.
# NB: adding 0s to flower_n for non-reproductive plants so that we can also include
# a flower_n * size interaction in the growth model as flowering branches dessicate
# and appear as shrinkage after fertilization. (Upon visual inspection, this 
# effect seems minimal, but will formally test anyway).
# 
# all_ramets <- all_ramets %>%
#   mutate(incr = abs(log_size_next - log_size),
#          col = case_when(
#            population == "Praia_de_Areao" & log_size < 3.5 & incr > 2.5 ~ "out",
#            TRUE ~ "normal")) %>%
#   filter(col != "normal") %>%
#   select(-(incr:col))

grow_data <- all_ramets %>%
  mutate(pop_id = paste(population, id, sep = "_")) %>%
  filter(!pop_id %in% buried_ramets$pop_id) %>%
  mutate(flower_n_grow = ifelse(is.na(flower_n), 0, flower_n)) %>%
  select(-pop_id)

# Store data

saveRDS(all_ramets, file = 'ipms/Data/all_ramets_di.rds')
saveRDS(grow_data, file  = 'ipms/Data/growth_data.rds')


# Now, get Seedling -> plant information

sdl_pts <- glue("ipms/Polygons/South_Africa/seedlings/rooisand_{c(1,3)}.shp")
sdls    <- lapply(sdl_pts, st_read) %>% 
  do.call(what = "rbind", args = .) %>%
  mutate(size_next = as.numeric(st_area(.)),
         log_size_next = log(size_next),
         plot = c(rep(1, 7), rep(3, 8))) %>%
  select(-size)

st_geometry(sdls) <- NULL

sdl_surv <- read.csv("ipms/Data/sdl_plots.csv") %>%
  right_join(sdls, by = "plot") %>%
  summarise(surv = n()/sum(unique((n_sdl))))

