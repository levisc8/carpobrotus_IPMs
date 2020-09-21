# Climate envelopes. Start with long-run averages

# clean_occ <- read.csv("repro_analysis/Data/gbif/gbif_occ_clean.csv",
#                       stringsAsFactors = FALSE)
# 
# all_sites <- read.csv("Data/All_Field_Sites.csv",
#                       stringsAsFactors = FALSE)
#
# mat     <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_01.tif")
# t_seas  <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_04.tif")
# t_co_mo <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_06.tif")
# t_co_qu <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_11.tif")
# map     <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_12.tif")
# p_seas  <- raster("repro_analysis/Data/chelsa_clim/CHELSA_bio10_15.tif")
# 
# coords <- cbind(lon = clean_occ$decimalLongitude,
#                 lat = clean_occ$decimalLatitude) %>%
#   as.matrix()

# clean_occ$mat     <- raster::extract(mat, coords, method = "bilinear")
# clean_occ$t_seas  <- raster::extract(t_seas, coords, method = "bilinear")
# clean_occ$t_co_mo <- raster::extract(t_co_mo, coords, method = "bilinear")
# clean_occ$t_co_qu <- raster::extract(t_co_qu, coords, method = "bilinear")
# clean_occ$map     <- raster::extract(map, coords, method = "bilinear")
# clean_occ$p_seas  <- raster::extract(p_seas, coords, method = "bilinear")
# 
# 
# site_coords <- cbind(all_sites$Lon,
#                      all_sites$Lat) %>%
#   as.matrix()
# 
# all_sites$mat     <- raster::extract(mat, site_coords, method = "bilinear")
# all_sites$t_seas  <- raster::extract(t_seas, site_coords, method = "bilinear")
# all_sites$t_co_mo <- raster::extract(t_co_mo, site_coords, method = "bilinear")
# all_sites$t_co_qu <- raster::extract(t_co_qu, site_coords, method = "bilinear")
# all_sites$map     <- raster::extract(map, site_coords, method = "bilinear")
# all_sites$p_seas  <- raster::extract(p_seas, site_coords, method = "bilinear")

# write.csv(clean_occ, file = "repro_analysis/Data/gbif/gbif_occ_clean_clim.csv",
#           row.names = FALSE)
# write.csv(all_sites, file = "repro_analysis/Data/all_sites_clim.csv",
#           row.names = FALSE)

clean_occ <- read.csv("repro_analysis/Data/gbif/gbif_occ_clean_clim.csv",
                      stringsAsFactors = FALSE) %>%
  select(species,
         decimalLongitude,
         decimalLatitude,
         mat:p_seas) %>%
  setNames(
    c(
      "species",
      "lon",
      "lat",
      "Mean_Annual_Temp",
      "Temp_Seasonality",
      "Temp_cold_Month",
      "Temp_Cold_Quarter",
      "Mean_Annual_Prec",
      "Prec_Seasonality"
    )
  ) %>%
  mutate(Sampled = "No")

all_sites <- read.csv("repro_analysis/Data/all_sites_clim.csv",
                      stringsAsFactors = FALSE) %>%
  mutate(species = "Carpobrotus spp.") %>%
  select(species, 
         Lon,
         Lat,
         mat:p_seas)%>%
  setNames(
    c(
      "species",
      "lon",
      "lat",
      "Mean_Annual_Temp",
      "Temp_Seasonality",
      "Temp_cold_Month",
      "Temp_Cold_Quarter",
      "Mean_Annual_Prec",
      "Prec_Seasonality"
    )
  ) %>%
  mutate(Sampled = "Yes")

for_plot <- rbind(clean_occ, all_sites) %>%
  select(Mean_Annual_Temp:Sampled)


pair_plot <- ggpairs(for_plot,
        aes(color = Sampled,
            shape = Sampled,
            size  = Sampled,
            alpha = Sampled),
        columns = 1:6) +
  scale_alpha_manual(
    values = c(0.2, 1)
  ) +
  scale_size_manual(values = c(0.5, 1.75)) +
  theme_bw()

print(pair_plot)


# Next, compute these values using weather data for the observed sites.
# dismo::biovars can do this!

