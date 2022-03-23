# Kriging environmental data for carpobrotus

source("ipms/R/utils/util-packages.R")
source("ipms/R/utils/util-clim.R")


# Downloads ----

all_occ <- read.csv("ipms/Data/all_gbif_field_sites.csv",
                    stringsAsFactors = FALSE) 

coords  <- data.frame(Lon = all_occ$lon, 
                      Lat = all_occ$lat, 
                      ID = paste0("p", seq_len(nrow(all_occ))))

# download_ERA citation:
# 1. Krigr - Kusch et al.
# 2. Muñoz Sabater, J., (2019): ERA5-Land hourly data from 1981 to present.
# Copernicus Climate Change Service (C3S) Climate Data Store (CDS). (Accessed on
# < 02-02-2022 >), 10.24381/cds.e2161bac
# Documented at:
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview

# temp <- download_ERA(Variable    = "2m_temperature", 
#                      DateStart   = "2017-04-01", 
#                      DateStop    = "2020-03-31",
#                      TResolution = "month", 
#                      Extent      = coords, 
#                      Buffer      = 0.1,
#                      ID          = "ID", 
#                      Dir         = "ipms/Data/weather/krigr", 
#                      API_User    = getOption("cds_user"), 
#                      API_Key     = getOption("cds_key"))
# 
# temp <- brick("ipms/Data/weather/krigr/2m_temperature_2017-04-01_2020-03-31_month.nc")
# # # NB: Prec is in mean precip/day for each month. 
# # 
# # prec <- download_ERA(Variable    = "total_precipitation", 
# #                      DateStart   = "2017-04-01", 
# #                      DateStop    = "2020-03-31",
# #                      TResolution = "month", 
# #                      Extent      = coords, 
# #                      Buffer      = 0.1,
# #                      ID          = "ID", 
# #                      Dir         = "ipms/Data/weather/krigr", 
# #                      API_User    = getOption("cds_user"), 
# #                      API_Key     = getOption("cds_key"))
# 
# prec <- brick("ipms/Data/weather/krigr/total_precipitation_2017-04-01_2020-03-31_month.nc")
# 
# # soil_water_1 <- download_ERA(Variable    = "volumetric_soil_water_layer_1", 
# #                              DateStart   = "2017-04-01", 
# #                              DateStop    = "2020-03-31",
# #                              TResolution = "month", 
# #                              Extent      = coords, 
# #                              Buffer      = 0.1,
# #                              ID          = "ID", 
# #                              Dir         = "ipms/Data/weather/krigr", 
# #                              API_User    = getOption("cds_user"), 
# #                              API_Key     = getOption("cds_key"))
# # 
# soil_water_1 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_1_2017-04-01_2020-03-31_month.nc")
# 
# # 
# # soil_water_2 <- download_ERA(Variable    = "volumetric_soil_water_layer_2", 
# #                              DateStart   = "2017-04-01", 
# #                              DateStop    = "2020-03-31",
# #                              TResolution = "month", 
# #                              Extent      = coords, 
# #                              Buffer      = 0.1,
# #                              ID          = "ID", 
# #                              Dir         = "ipms/Data/weather/krigr", 
# #                              API_User    = getOption("cds_user"), 
# #                              API_Key     = getOption("cds_key"))
# # 
# soil_water_2 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_2_2017-04-01_2020-03-31_month.nc")
# 
# # 
# # soil_water_3 <- download_ERA(Variable    = "volumetric_soil_water_layer_3", 
# #                              DateStart   = "2017-04-01", 
# #                              DateStop    = "2020-03-31",
# #                              TResolution = "month", 
# #                              Extent      = coords, 
# #                              Buffer      = 0.1,
# #                              ID          = "ID", 
# #                              Dir         = "ipms/Data/weather/krigr", 
# #                              API_User    = getOption("cds_user"), 
# #                              API_Key     = getOption("cds_key"))
# soil_water_3 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_3_2017-04-01_2020-03-31_month.nc")
# 
# 
# # covs <- download_DEM(
# #   Train_ras      = temp,
# #   Target_res     = 0.02,
# #   Keep_Temporary = TRUE,
# #   Shape          = coords,
# #   Buffer         = 0.1,
# #   ID             = "ID",
# #   Dir            = "ipms/Data/weather/krigr"
# # )
# 
# covs <- list(
#   brick("ipms/Data/weather/krigr/GMTED2010_Train.nc"),
#   brick("ipms/Data/weather/krigr/GMTED2010_Target.nc")  
# )
# 
# 
# month_dict <- expand.grid(month = month.name,
#                           year = 2017:2020,
#                           stringsAsFactors = FALSE)
# month_dict <- month_dict[-c(1:3, 40:48), ]
# 
# mnts <- apply(month_dict, 1, function(x) paste0(x, collapse = "_"))
# 
# start <- Sys.time()
# temp_krig_data <- krigR(
#   Data              = temp,
#   Covariates_coarse = covs[[1]],
#   Covariates_fine   = covs[[2]],
#   KrigingEquation   = Data ~ layer,
#   Keep_Temporary    = FALSE,
#   Cores             = 4L,
#   Dir               = "ipms/Data/weather/krigr",
#   FileName          = "temp_kriged.nc",
#   nmax              = 15
# )
# 
# temp_krig_pts <- extract_krig(temp_krig_data, 
#                               coords, 
#                               "ID",
#                               mnts)
# 
# tim <- round(Sys.time() - start, 3)
# print(tim)
# 
# saveRDS(temp_krig_pts, file = "ipms/Data/weather/krigr/kriged_temp_pts.rds")
# saveRDS(temp_krig_data, file = "ipms/Data/weather/krigr/kriged_temp_data.rds")
# 
# 
# start <- Sys.time()
# prec_krig_data <- krigR(
#   Data              = prec,
#   Covariates_coarse = covs[[1]],
#   Covariates_fine   = covs[[2]],
#   KrigingEquation   = Data ~ layer,
#   Keep_temporary    = FALSE,
#   Cores             = 4L,
#   Dir               = "ipms/Data/weather/krigr",
#   FileName          = "prec_kriged.nc",
#   nmax              = 15
# )
# 
# 
# prec_krig_pts <- extract_krig(prec_krig_data, 
#                               coords, 
#                               "ID",
#                               mnts)
# 
# tim <- round(Sys.time() - start, 3)
# print(tim)
# 
# saveRDS(prec_krig_pts, file = "ipms/Data/weather/krigr/kriged_prec_pts.rds")
# saveRDS(prec_krig_data, file = "ipms/Data/weather/krigr/kriged_prec_data.rds")
# 
# start <- Sys.time()
# 
# soil_water_1_krig_data <- krigR(
#   Data              = soil_water_1,
#   Covariates_coarse = covs[[1]],
#   Covariates_fine   = covs[[2]],
#   KrigingEquation   = Data ~ layer,
#   Keep_Temporary    = FALSE,
#   Cores             = 4L,
#   Dir               = "ipms/Data/weather/krigr",
#   FileName          = "soil_water_1_kriged.nc",
#   nmax              = 15
# )
# 
# 
# sw1_krig_pts <- extract_krig(soil_water_1_krig_data, 
#                               coords, 
#                               "ID",
#                               mnts)
# 
# tim <- round(Sys.time() - start, 3)
# print(tim)
# 
# saveRDS(sw1_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw1_pts.rds")
# saveRDS(soil_water_1_krig_data, file = "ipms/Data/weather/krigr/kriged_sw1_data.rds")
# 
# 
# start <- Sys.time()
# soil_water_2_krig_data <- krigR(
#   Data              = soil_water_2,
#   Covariates_coarse = covs[[1]],
#   Covariates_fine   = covs[[2]],
#   KrigingEquation   = Data ~ layer,
#   Keep_Temporary    = FALSE,
#   Cores             = 4L,
#   Dir               = "ipms/Data/weather/krigr",
#   FileName          = "soil_water_2_kriged.nc",
#   nmax              = 15
# )
# 
# 
# sw2_krig_pts <- extract_krig(soil_water_2_krig_data, 
#                              coords, 
#                              "ID",
#                              mnts)
# 
# tim <- round(Sys.time() - start, 3)
# print(tim)
# 
# 
# saveRDS(sw2_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw2_pts.rds")
# saveRDS(soil_water_2_krig_data, file = "ipms/Data/weather/krigr/kriged_sw2_data.rds")
# 
# start <- Sys.time()
# 
# soil_water_3_krig_data <- krigR(
#   Data              = soil_water_3,
#   Covariates_coarse = covs[[1]],
#   Covariates_fine   = covs[[2]],
#   KrigingEquation   = Data ~ layer,
#   Keep_Temporary    = FALSE,
#   Cores             = 4L,
#   Dir               = "ipms/Data/weather/krigr",
#   FileName          = "soil_water_3_kriged.nc",
#   nmax              = 15
# )
# 
# sw3_krig_pts <- extract_krig(soil_water_3_krig_data, 
#                              coords, 
#                              "ID",
#                              mnts)
# 
# tim <- round(Sys.time() - start, 3)
# print(tim)"
# 
# saveRDS(soil_water_3_krig_data, file = "ipms/Data/weather/krigr/kriged_sw3_data.rds")
# saveRDS(sw3_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw3_pts.rds")


# tidy_clim_data() centers/scales each variable w/r/t to global range, and 
# computes variables we need.
# I can't remember why, but the ID column of the kriged values is integer (
# missing its "p", so that's an intermediate step). Next, we join with coordinates
# so that we can use latitude to determine wet/dry seasons. After that,
# we use tidy_clim_data to:
# 1. convert to long form, and split month/year into separate columns
# 2. aggregate values w/ summarise by year and by year * season
# 3. rejoin w/ IDs
# 4. Voila. see util-clim.R - tidy_climd_data() and .aggregate_data() for 
#    details.

source("ipms/R/utils/util-packages2.R")

temp_krig_pts <- readRDS("ipms/Data/weather/krigr/kriged_temp_pts.rds") %>%
  lapply(function(x) mutate(x, ID = paste0("p", seq_len(nrow(x))))) %>%
  lapply(function(x, coords) left_join(x, coords, by = "ID"), coords = coords) %>%
  tidy_clim_data("temp")

prec_krig_pts <- readRDS("ipms/Data/weather/krigr/kriged_prec_pts.rds") %>%
  lapply(function(x) mutate(x, ID = paste0("p", seq_len(nrow(x))))) %>%
  lapply(function(x, coords) left_join(x, coords, by = "ID"), coords = coords) %>%
  tidy_clim_data("prec")

sw1_krig_pts <- readRDS("ipms/Data/weather/krigr/kriged_sw1_pts.rds") %>%
  lapply(function(x) mutate(x, ID = paste0("p", seq_len(nrow(x))))) %>%
  lapply(function(x, coords) left_join(x, coords, by = "ID"), coords = coords) %>%
  tidy_clim_data("sw1")

sw2_krig_pts <- readRDS("ipms/Data/weather/krigr/kriged_sw2_pts.rds") %>%
  lapply(function(x) mutate(x, ID = paste0("p", seq_len(nrow(x))))) %>%
  lapply(function(x, coords) left_join(x, coords, by = "ID"), coords = coords) %>%
  tidy_clim_data("sw2")

sw3_krig_pts <- readRDS("ipms/Data/weather/krigr/kriged_sw3_pts.rds") %>%
  lapply(function(x) mutate(x, ID = paste0("p", seq_len(nrow(x))))) %>%
  lapply(function(x, coords) left_join(x, coords, by = "ID"), coords = coords) %>%
  tidy_clim_data("sw3")

# Now, we need to join all this data together and filter it down to our sites.

all_seas_data <- left_join(temp_krig_pts$seasonal, 
                           prec_krig_pts$seasonal,
                           by = "ID") %>%
  left_join(sw1_krig_pts$seasonal, by = "ID") %>%
  left_join(sw2_krig_pts$seasonal, by = "ID") %>%
  left_join(sw3_krig_pts$seasonal, by = "ID")

all_annual_data <- left_join(temp_krig_pts$annual, 
                             prec_krig_pts$annual,
                             by = "ID") %>%
  left_join(sw1_krig_pts$annual, by = "ID") %>%
  left_join(sw2_krig_pts$annual, by = "ID") %>%
  left_join(sw3_krig_pts$annual, by = "ID")

`%between%` <- function(x, y) {
  
  x >= y[1] & x <= y[2]
  
}

seas_17 <- select(all_seas_data, c(ID, Lat, Lon), ends_with("2017")) %>%
  filter(complete.cases(.)) %>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  ))
seas_18 <- select(all_seas_data, c(ID, Lat, Lon), ends_with("2018")) %>%
  filter(complete.cases(.)) %>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  ))
seas_19 <- select(all_seas_data, c(ID, Lat, Lon), ends_with("2019")) %>%
  filter(complete.cases(.))%>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  ))

ann_17 <- select(all_annual_data, c(ID, Lat, Lon), ends_with("2017")) %>%
  filter(complete.cases(.)) %>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  )) %>%
  select(-c(Lat, Lon))
ann_18 <- select(all_annual_data, c(ID, Lat, Lon), ends_with("2018")) %>%
  filter(complete.cases(.)) %>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  )) %>%
  select(-c(Lat, Lon))
ann_19 <- select(all_annual_data, c(ID, Lat, Lon), ends_with("2019")) %>%
  filter(complete.cases(.)) %>% 
  mutate(native = case_when(
    Lat %between% c(-40, -20) & Lon %between% c(0, 45) ~ "native",
    TRUE ~ "invasive"
  )) %>%
  select(-c(Lat, Lon))

# PCAs -----
# TODO: Figure out how to change shape/color/size/something for the sites
# that we actually sampled for demography.

library(ggplot2)
library(ggfortify)

seas_17_pca <- princomp(seas_17[ , 4:13])
seas_18_pca <- princomp(seas_18[ , 4:13])
seas_19_pca <- princomp(seas_19[ , 4:13])

ann_17_pca <- princomp(ann_17[ , 2:11])
ann_18_pca <- princomp(ann_18[ , 2:11])
ann_19_pca <- princomp(ann_19[ , 2:11])


seas_17_pca_gg <-  autoplot(seas_17_pca,
                            data = seas_17,
                            loadings = TRUE,
                            loadings.label = TRUE, 
                            loadings.label.size = 6, 
                            colour = "native",
                            alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                       values = viridis::inferno(2, end = 0.7)) +
  theme(legend.position = "none") +
  ggtitle("Seasonally aggregated data",
          subtitle = "2017")
  
seas_18_pca_gg <- autoplot(seas_18_pca,
                        data = seas_18,
                        loadings = TRUE,
                        loadings.label = TRUE, 
                        loadings.label.size = 6, 
                        colour = "native",
                        alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                     values = viridis::inferno(2, end = 0.7)) +
  theme(legend.position = "none") +
  ggtitle("",
          subtitle = "2018")

seas_19_pca_gg <- autoplot(seas_19_pca,
                        data = seas_19,
                        loadings = TRUE,
                        loadings.label = TRUE, 
                        loadings.label.size = 6, 
                        colour = "native",
                        alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                     values = viridis::inferno(2, end = 0.7))  +
  ggtitle("",
          subtitle = "2019")

ann_17_pca_gg <- autoplot(ann_17_pca,
                       data = ann_17,
                       loadings = TRUE,
                       loadings.label = TRUE, 
                       loadings.label.size = 6, 
                       colour = "native",
                       alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                     values = viridis::inferno(2, end = 0.7))  +
  theme(legend.position = "none") +
  ggtitle("Annually aggregated data",
          subtitle = "2017")

ann_18_pca_gg <- autoplot(ann_18_pca,
                       data = ann_18,
                       loadings = TRUE,
                       loadings.label = TRUE, 
                       loadings.label.size = 6, 
                       colour = "native",
                       alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                     values = viridis::inferno(2, end = 0.7))  +
  theme(legend.position = "none") +
  ggtitle("",
          subtitle = "2018")

ann_19_pca_gg <- autoplot(ann_19_pca,
                       data = ann_19,
                       loadings = TRUE,
                       loadings.label = TRUE, 
                       loadings.label.size = 6, 
                       colour = "native",
                       alpha = 0.3) + 
  theme_bw() +
  scale_color_manual(breaks = c("native", "invasive"),
                     values = viridis::inferno(2, end = 0.7))  +
  ggtitle("",
          subtitle = "2019")

pdf("ipms/Figures/clim/seasonal_annual_krig_pcas.pdf",
    width = 18, height = 9)
  grid.arrange(seas_17_pca_gg, seas_18_pca_gg, seas_19_pca_gg,
               ann_17_pca_gg,  ann_18_pca_gg,  ann_19_pca_gg,
               layout_matrix = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE))

dev.off()
# Not sure how to approach 2020 data yet - these only apply to PT sites, and
# only for a short period of time....



seas_20 <- select(all_seas_data, c(ID, Lat, Lon), ends_with("2020")) 
ann_20 <- select(all_annual_data, c(ID), ends_with("2020"))

seas_17$seas_pc_1 <- seas_17_pca$scores[ , 1]
seas_18$seas_pc_1 <- seas_18_pca$scores[ , 1]
seas_19$seas_pc_1 <- seas_19_pca$scores[ , 1]
seas_17$seas_pc_2 <- seas_17_pca$scores[ , 2]
seas_18$seas_pc_2 <- seas_18_pca$scores[ , 2]
seas_19$seas_pc_2 <- seas_19_pca$scores[ , 2]

ann_17$ann_pc_1 <- ann_17_pca$scores[ , 1] 
ann_18$ann_pc_1 <- ann_18_pca$scores[ , 1]
ann_19$ann_pc_1 <- ann_19_pca$scores[ , 1]
ann_17$ann_pc_2 <- ann_17_pca$scores[ , 2]
ann_18$ann_pc_2 <- ann_18_pca$scores[ , 2]
ann_19$ann_pc_2 <- ann_19_pca$scores[ , 2]  


# Extract site data. We no longer need the envelope stuff since everything is
# now computed.

site_data <- select(all_occ, species:lat, Sampled) %>%
  mutate(ID = paste0("p", seq_len(nrow(.)))) %>%
  filter(Sampled == "Yes") 

all_17 <- left_join(seas_17, ann_17, by = "ID") %>%
  filter(ID %in% site_data$ID) %>%
  left_join(site_data, by = "ID")
all_18 <- left_join(seas_18, ann_18, by = "ID") %>%
  filter(ID %in% site_data$ID) %>%
  left_join(site_data, by = "ID")
all_19 <- left_join(seas_19, ann_19, by = "ID") %>%
  filter(ID %in% site_data$ID) %>%
  left_join(site_data, by = "ID")

# Create temporary
temp_t_1 <- filter(all_19, site == "acsa") %>%
  select(ID, site, temp_dry_2019:sw3_wet_2019, mean_temp_2019:seas_sw3_2019) %>%
  setNames(c(gsub("2019", "t_1", names(.)))) %>%
  mutate(
    seas_pc_1_t_1 = NA,
    seas_pc_2_t_1 = NA,
    ann_pc_1_t_1  = NA,
    ann_pc_2_t_1  = NA
  )

temp_t <- temp_t_1 %>%
  setNames(gsub("t_1$", "t", names(.))) 

for(i in seq_len(nrow(site_data))) {

  site_id <- site_data$ID[i]
  
  if(site_data$site[i] == "Havatselet") {
    
    use_t_1 <- all_17 %>% 
      filter(ID == site_id) %>% 
      select(ID:sw3_wet_2017, 
             mean_temp_2017:seas_sw3_2017) %>%
      setNames(gsub("2017", "t_1", names(.)))
    
    use_t <- all_18 %>% 
      filter(ID == site_id) %>% 
      select(ID:sw3_wet_2018, 
             mean_temp_2018:seas_sw3_2018) %>%
      setNames(gsub("2018", "t", names(.)))
    
    use_t_1$seas_pc_1_t_1 <- all_17$seas_pc_1[all_17$ID == site_id]
    use_t_1$ann_pc_1_t_1  <- all_17$ann_pc_1[all_17$ID == site_id]
    use_t$seas_pc_1_t     <- all_18$ann_pc_1[all_18$ID == site_id]
    use_t$ann_pc_1_t      <- all_18$ann_pc_1[all_18$ID == site_id]

    use_t_1$seas_pc_2_t_1 <- all_17$seas_pc_2[all_17$ID == site_id]
    use_t_1$ann_pc_2_t_1  <- all_17$ann_pc_2[all_17$ID == site_id]
    use_t$seas_pc_2_t     <- all_18$ann_pc_2[all_18$ID == site_id]
    use_t$ann_pc_2_t      <- all_18$ann_pc_2[all_18$ID == site_id]
    
        
  } else {
    use_t_1 <- all_18 %>% 
      filter(ID == site_id) %>% 
      select(ID:sw3_wet_2018, 
             mean_temp_2018:seas_sw3_2018) %>%
      setNames(gsub("2018", "t_1", names(.)))
    
    use_t <- all_19 %>%
      filter(ID == site_id) %>% 
      select(ID:sw3_wet_2019, 
             mean_temp_2019:seas_sw3_2019) %>%
      setNames(gsub("2019", "t", names(.)))
    
    use_t_1$seas_pc_1_t_1 <- all_18$seas_pc_1[all_18$ID == site_id]
    use_t_1$ann_pc_1_t_1  <- all_18$ann_pc_1[all_18$ID == site_id]
    use_t$seas_pc_1_t     <- all_19$ann_pc_1[all_19$ID == site_id]
    use_t$ann_pc_1_t      <- all_19$ann_pc_1[all_19$ID == site_id]
    
    use_t_1$seas_pc_2_t_1 <- all_18$seas_pc_2[all_18$ID == site_id]
    use_t_1$ann_pc_2_t_1  <- all_18$ann_pc_2[all_18$ID == site_id]
    use_t$seas_pc_2_t     <- all_19$ann_pc_2[all_19$ID == site_id]
    use_t$ann_pc_2_t      <- all_19$ann_pc_2[all_19$ID == site_id]
  } 
  
  temp_t_1 <- rbind(temp_t_1, use_t_1)
  temp_t   <- rbind(temp_t, use_t)

}

all_clim <- left_join(site_data, temp_t_1,  by = "ID") %>%
  left_join(temp_t, by = "ID") %>% 
  select(-c(Lon.x, Lon.y, Lat.x, Lat.y))

write.csv(all_clim, 
          file = "ipms/Data/all_kriged_field_sites.csv",
          row.names = FALSE)
