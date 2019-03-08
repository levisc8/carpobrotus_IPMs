# Basic shape file analysis

# Unfortunately, populations use different CRS's localized to their lat/longs
# within WGS 84. These are more accurate, but also mean that we can't just
# combine them. Thus, for the hierarchical models for each vital rate, I'm
# going to create a separate data frame that drops the spatial info (but still
# designates who belongs to what population), and retain a list of data frames
# that still have the polygons embedded in them. These can be used for any 
# kind of spatially explicit analysis later (e.g. density dependence)

all_sites <- read.csv('Data/All_Field_Sites.csv', stringsAsFactors = FALSE)

pops <- all_sites[all_sites$Demography == 1 &
                    all_sites$Notes == 'ready', 
                  c('Country', 'Site')]

all_polygons <- list()
for(i in seq_along(pops[ ,1])) {
  
  pop <- pops$Site[i]
  
  # This will generate both ramet and genet level data frames for each population
  temp <- infile_data(
    glue('Polygons/{pops$Country[i]}/{pop}/Polygon_Final.shp'),
    pop
  )
  temp$country <- all_sites$Country[all_sites$Site == pop][1]
  all_polygons <- splice(all_polygons, temp) 
  names(all_polygons)[i] <- pop
  
}

no_polygons <- map(all_polygons, 
                   .f = function(x) {
                     st_geometry(x) <- NULL
                     
                     return(x)
                   }
)

all_ramets <- do.call(rbind, no_polygons)

all_genets <- all_ramets %>%
  ramet_to_genet(population, clone_of)

# next, start fitting models and producing preliminary figures
