# Climate envelopes. Start with long-run averages

# clean_occ <- read.csv("repro_analysis/Data/gbif/gbif_occ_clean.csv",
#                       stringsAsFactors = FALSE)
# 
# all_sites <- read.csv("Data/All_Field_Sites.csv",
#                       stringsAsFactors = FALSE) %>%
#   filter(Demography == 1)
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
  mutate(Site = NA_character_) %>%
  select(species,
         Site,
         decimalLongitude,
         decimalLatitude,
         mat:p_seas) %>%
  setNames(
    c(
      "species",
      "site",
      "lon",
      "lat",
      "Mean_Annual_Temp_hist",
      "Temp_Seasonality_hist",
      "Temp_cold_Month_hist",
      "Temp_Cold_Quarter_hist",
      "Mean_Annual_Prec_hist",
      "Prec_Seasonality_hist"
    )
  ) %>%
  mutate(Sampled = "No")

all_sites <- read.csv("repro_analysis/Data/all_sites_clim.csv",
                      stringsAsFactors = FALSE) %>%
  mutate(species = "Carpobrotus spp.") %>%
  select(species, 
         Site,
         Lon,
         Lat,
         mat:p_seas) %>%
  setNames(
    c(
      "species",
      "site",
      "lon",
      "lat",
      "Mean_Annual_Temp_hist",
      "Temp_Seasonality_hist",
      "Temp_cold_Month_hist",
      "Temp_Cold_Quarter_hist",
      "Mean_Annual_Prec_hist",
      "Prec_Seasonality_hist"
    )
  ) %>%
  mutate(Sampled = "Yes")

for_plot <- rbind(clean_occ, all_sites) %>%
  select(Mean_Annual_Temp_hist:Sampled)


old_plot <- ggpairs(for_plot,
        aes(color = Sampled,
            shape = Sampled,
            size  = Sampled,
            alpha = Sampled),
        columns = 1:6) +
  scale_alpha_manual(
    values = c(0.2, 1)
  ) +
  scale_size_manual(values = c(0.5, 2)) +
  scale_color_manual(values = c("black", "red")) + 
  theme_bw()


png(filename = "repro_analysis/Manuscript/Figures/clim_coverage_1979-2013.png",
    height = 10,
    width = 10,
    units = "in",
    res = 300)

  print(old_plot + ggtitle("Climate 1979-2013"))

dev.off()

# Next, compute these values for the range of occurrence data + observed
# sites using worldClim data, and dismo::biovars(). The code below is run on
# iDiv's RStudio server, but keeping here for documentation.
# Worldclim URL: https://www.worldclim.org/data/bioclim.html
# 
# years <- 2016:2018
# 
# months <- c(paste(0, 1:9, sep  = ""), 10:12)
# 
# full_grid <- expand.grid(years, months, stringsAsFactors = FALSE)
# 
# base_path <- c("repro_analysis/Data/weather/wc2.1_2.5m")
# 
# var <- "tmin"
# t_min_paths <- glue("{base_path}_{var}_2010-2018/wc2.1_2.5m_{var}_{full_grid[,1]}-{full_grid[,2]}.tif") %>%
#   sort()
# 
# var <- "tmax"
# t_max_paths <- glue("{base_path}_{var}_2010-2018/wc2.1_2.5m_{var}_{full_grid[,1]}-{full_grid[,2]}.tif") %>%
#   sort()
# 
# var <- "prec"
# prec_paths <- glue("{base_path}_{var}_2010-2018/wc2.1_2.5m_{var}_{full_grid[,1]}-{full_grid[,2]}.tif") %>%
#   sort()
# 
# # Compute annual values for each layer. Next, we'll average over those
# # to get a single 2016-2018 value for the variables we computed above.
# 
# bio_vars <- list()
# 
# for(i in seq_along(years)) {
#   
#   
#   yr <- years[i]
#   
#   message("\nBegin processing for: ", yr, "\n")
# 
#   # All file paths have "2018" in them, so we use a more laborious but
#   # exact process of constructing file paths for them to make sure our stacks
#   # only contain the 12 months of 2018
# 
#   if(yr == 2018) {
#     
#     base_path <- "repro_analysis/Data/weather"
#     months <- c(paste("0", 1:9, sep = ""), 10:12)
#     months <- paste("2018-", months, ".tif", sep = "")
#     
#     base_tmax <- glue("{base_path}/wc2.1_2.5m_tmax_2010-2018/wc2.1_2.5m_tmax_")
#     base_tmin <- glue("{base_path}/wc2.1_2.5m_tmin_2010-2018/wc2.1_2.5m_tmin_")
#     base_prec <- glue("{base_path}/wc2.1_2.5m_prec_2010-2018/wc2.1_2.5m_prec_")
#     
#     use_tmax <- as.list(paste(base_tmax, months, sep = ""))
#     use_tmin <- as.list(paste(base_tmin, months, sep = ""))
#     use_prec <- as.list(paste(base_prec, months, sep = ""))
#     
#   } else{
#     
#     use_tmin <- as.list(t_min_paths[grepl(yr, t_min_paths)])
#     use_tmax <- as.list(t_max_paths[grepl(yr, t_max_paths)])
#     use_prec <- as.list(prec_paths[grepl(yr, prec_paths)])
#     
#   }
#   
#   min_stack <- raster::stack(use_tmin)
#   max_stack <- raster::stack(use_tmax)
#   pre_stack <- raster::stack(use_prec)
#   
#   nm <- paste("bioclim_", yr, sep = "")
#   
#   bio_vars[[i]] <- biovars(pre_stack, min_stack, max_stack)
#   
#   names(bio_vars)[i] <- nm
#   
#   message("\nEnd processing for: ", yr, "\n")
#   
# }
# 
# # Save computed values 
# 
# saveRDS(bio_vars, 
#         file = "repro_analysis/Data/weather/computed_bioclim_from_weather_data.rds")

# bio_vars <- readRDS("repro_analysis/Data/weather/computed_bioclim_from_weather_data.rds")
# 
# 
# # Next, we want to pluck out the variables we need. These will be layers
# # 1, 4, 6, 11, 12, 15 from the brick that comes out biovars()
# 
# layer_to_stack <- function(layer_list, layer) {
#   
#   lapply(layer_list,
#          function(x, layer) {
#            subset(x, subset = layer)
#          },
#          layer = layer) %>%
#     stack(x = .)
#   
# }
# 
# 
# mat     <- layer_to_stack(bio_vars, 1) 
# t_seas  <- layer_to_stack(bio_vars, 4) 
# t_co_mo <- layer_to_stack(bio_vars, 6) 
# t_co_qu <- layer_to_stack(bio_vars, 11) 
# map     <- layer_to_stack(bio_vars, 12)
# p_seas  <- layer_to_stack(bio_vars, 15) 
# 
# # Once we have those layers arranged such that each variable has it's own list,
# # we can use calc(var_list, fun = mean)
# 
# mat <- calc(mat, fun = mean)
# t_seas <- calc(t_seas, fun = mean)
# t_co_mo <- calc(t_co_mo, fun = mean)
# t_co_qu <- calc(t_co_qu, fun = mean)
# map <- calc(map, fun = mean)
# p_seas <- calc(p_seas, fun = mean)
# 
# mat     <- raster("repro_analysis/Data/weather/mat_2016-2018.tif")
# t_seas  <- raster("repro_analysis/Data/weather/t_seas_2016-2018.tif")
# t_co_mo <- raster("repro_analysis/Data/weather/t_co_mo_2016-2018.tif")
# t_co_qu <- raster("repro_analysis/Data/weather/t_co_qu_2016-2018.tif")
# map     <- raster("repro_analysis/Data/weather/map_2016-2018.tif")
# p_seas  <- raster("repro_analysis/Data/weather/p_seas_2016-2018.tif")
# 
# # Finally, extract values for each computed variable layer X occurrence/field site
# 
# all_data <- rbind(clean_occ, all_sites)
# 
# coords <- cbind(all_data$lon,
#                 all_data$lat)
# 
# all_data$mat_rec     <- raster::extract(mat, coords, method = "bilinear")
# all_data$t_seas_rec  <- raster::extract(t_seas, coords, method = "bilinear")
# all_data$t_co_mo_rec <- raster::extract(t_co_mo, coords, method = "bilinear")
# all_data$t_co_qu_rec <- raster::extract(t_co_qu, coords, method = "bilinear")
# all_data$map_rec     <- raster::extract(map, coords, method = "bilinear")
# all_data$p_seas_rec  <- raster::extract(p_seas, coords, method = "bilinear")
# 
# write.csv(all_data,
#           file = "repro_analysis/Data/all_gbif_field_sites.csv",
#           row.names = FALSE)

all_data <- read.csv("repro_analysis/Data/all_gbif_field_sites.csv",
                     stringsAsFactors = FALSE)

# Now, we can plot everything!

for_plot <- select(all_data, site, Sampled:p_seas_rec)

recent_plot <- ggpairs(for_plot,
                       aes(color = Sampled,
                           shape = Sampled,
                           size  = Sampled,
                           alpha = Sampled),
                       columns = 3:8) +
  scale_alpha_manual(
    values = c(0.2, 1)
  ) +
  scale_size_manual(values = c(0.5, 2)) +
  scale_color_manual(values = c("black", "red")) + 
  theme_bw()

png(filename = "repro_analysis/Manuscript/Figures/clim_coverage_2016-2018.png",
    height = 10,
    width = 10,
    units = "in",
    res = 300)

  print(recent_plot + ggtitle("Climate 2016-2018"))

dev.off()

