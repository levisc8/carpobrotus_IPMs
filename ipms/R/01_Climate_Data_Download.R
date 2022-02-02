# Kriging environmental data for carpobrotus

source("ipms/R/utils/util-packages.R")
source("ipms/R/utils/util-clim.R")


# Downloads ----

all_occ <- read.csv("ipms/Data/all_gbif_field_sites.csv",
                    stringsAsFactors = FALSE) #%>%
  # select(species:lat, Sampled)

coords  <- data.frame(Lon = all_occ$lon, 
                      Lat = all_occ$lat, 
                      ID = paste0("p", seq_len(nrow(all_occ))))

# download_ERA citation:
# 1. Krigr - Kusch et al.
# 2. Muñoz Sabater, J., (2019): ERA5-Land hourly data from 1981 to present.
# Copernicus Climate Change Service (C3S) Climate Data Store (CDS). (Accessed on
# < 02-02-2022 >), 10.24381/cds.e2161bac

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
temp <- brick("ipms/Data/weather/krigr/2m_temperature_2017-04-01_2020-03-31_month.nc")
# # NB: Prec is in mean precip/day for each month. 
# 
# prec <- download_ERA(Variable    = "total_precipitation", 
#                      DateStart   = "2017-04-01", 
#                      DateStop    = "2020-03-31",
#                      TResolution = "month", 
#                      Extent      = coords, 
#                      Buffer      = 0.1,
#                      ID          = "ID", 
#                      Dir         = "ipms/Data/weather/krigr", 
#                      API_User    = getOption("cds_user"), 
#                      API_Key     = getOption("cds_key"))

prec <- brick("ipms/Data/weather/krigr/total_precipitation_2017-04-01_2020-03-31_month.nc")

# soil_water_1 <- download_ERA(Variable    = "volumetric_soil_water_layer_1", 
#                              DateStart   = "2017-04-01", 
#                              DateStop    = "2020-03-31",
#                              TResolution = "month", 
#                              Extent      = coords, 
#                              Buffer      = 0.1,
#                              ID          = "ID", 
#                              Dir         = "ipms/Data/weather/krigr", 
#                              API_User    = getOption("cds_user"), 
#                              API_Key     = getOption("cds_key"))
# 
soil_water_1 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_1_2017-04-01_2020-03-31_month.nc")

# 
# soil_water_2 <- download_ERA(Variable    = "volumetric_soil_water_layer_2", 
#                              DateStart   = "2017-04-01", 
#                              DateStop    = "2020-03-31",
#                              TResolution = "month", 
#                              Extent      = coords, 
#                              Buffer      = 0.1,
#                              ID          = "ID", 
#                              Dir         = "ipms/Data/weather/krigr", 
#                              API_User    = getOption("cds_user"), 
#                              API_Key     = getOption("cds_key"))
# 
soil_water_2 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_2_2017-04-01_2020-03-31_month.nc")

# 
# soil_water_3 <- download_ERA(Variable    = "volumetric_soil_water_layer_3", 
#                              DateStart   = "2017-04-01", 
#                              DateStop    = "2020-03-31",
#                              TResolution = "month", 
#                              Extent      = coords, 
#                              Buffer      = 0.1,
#                              ID          = "ID", 
#                              Dir         = "ipms/Data/weather/krigr", 
#                              API_User    = getOption("cds_user"), 
#                              API_Key     = getOption("cds_key"))
soil_water_3 <- brick("ipms/Data/weather/krigr/volumetric_soil_water_layer_3_2017-04-01_2020-03-31_month.nc")


# covs <- download_DEM(
#   Train_ras      = temp,
#   Target_res     = 0.02,
#   Keep_Temporary = TRUE,
#   Shape          = coords,
#   Buffer         = 0.1,
#   ID             = "ID",
#   Dir            = "ipms/Data/weather/krigr"
# )

covs <- list(
  brick("ipms/Data/weather/krigr/GMTED2010_Train.nc"),
  brick("ipms/Data/weather/krigr/GMTED2010_Target.nc")  
)


month_dict <- expand.grid(month = month.name,
                          year = 2017:2020,
                          stringsAsFactors = FALSE)
month_dict <- month_dict[-c(1:3, 40:48), ]

mnts <- apply(month_dict, 1, function(x) paste0(x, collapse = "_"))

start <- Sys.time()
temp_krig_data <- krigR(
  Data              = temp,
  Covariates_coarse = covs[[1]],
  Covariates_fine   = covs[[2]],
  KrigingEquation   = Data ~ layer,
  Keep_Temporary    = FALSE,
  Cores             = 4L,
  Dir               = "ipms/Data/weather/krigr",
  FileName          = "temp_kriged.nc",
  nmax              = 15
)

temp_krig_pts <- extract_krig(temp_krig_data, 
                              coords, 
                              "ID",
                              mnts)

tim <- round(Sys.time() - start, 3)
print(tim)

saveRDS(temp_krig_pts, file = "ipms/Data/weather/krigr/kriged_temp_pts.rds")
saveRDS(temp_krig_data, file = "ipms/Data/weather/krigr/kriged_temp_data.rds")


start <- Sys.time()
prec_krig_data <- krigR(
  Data              = prec,
  Covariates_coarse = covs[[1]],
  Covariates_fine   = covs[[2]],
  KrigingEquation   = Data ~ layer,
  Keep_temporary    = FALSE,
  Cores             = 4L,
  Dir               = "ipms/Data/weather/krigr",
  FileName          = "prec_kriged.nc",
  nmax              = 15
)


prec_krig_pts <- extract_krig(prec_krig_data, 
                              coords, 
                              "ID",
                              mnts)

tim <- round(Sys.time() - start, 3)
print(tim)

saveRDS(prec_krig_pts, file = "ipms/Data/weather/krigr/kriged_prec_pts.rds")
saveRDS(prec_krig_data, file = "ipms/Data/weather/krigr/kriged_prec_data.rds")

start <- Sys.time()

soil_water_1_krig_data <- krigR(
  Data              = soil_water_1,
  Covariates_coarse = covs[[1]],
  Covariates_fine   = covs[[2]],
  KrigingEquation   = Data ~ layer,
  Keep_Temporary    = FALSE,
  Cores             = 4L,
  Dir               = "ipms/Data/weather/krigr",
  FileName          = "soil_water_1_kriged.nc",
  nmax              = 15
)


sw1_krig_pts <- extract_krig(soil_water_1_krig_data, 
                              coords, 
                              "ID",
                              mnts)

tim <- round(Sys.time() - start, 3)
print(tim)

saveRDS(sw1_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw1_pts.rds")
saveRDS(soil_water_1_krig_data, file = "ipms/Data/weather/krigr/kriged_sw1_data.rds")


start <- Sys.time()
soil_water_2_krig_data <- krigR(
  Data              = soil_water_2,
  Covariates_coarse = covs[[1]],
  Covariates_fine   = covs[[2]],
  KrigingEquation   = Data ~ layer,
  Keep_Temporary    = FALSE,
  Cores             = 4L,
  Dir               = "ipms/Data/weather/krigr",
  FileName          = "soil_water_2_kriged.nc",
  nmax              = 15
)


sw2_krig_pts <- extract_krig(soil_water_2_krig_data, 
                             coords, 
                             "ID",
                             mnts)

tim <- round(Sys.time() - start, 3)
print(tim)


saveRDS(sw2_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw2_pts.rds")
saveRDS(soil_water_2_krig_data, file = "ipms/Data/weather/krigr/kriged_sw2_data.rds")

start <- Sys.time()

soil_water_3_krig_data <- krigR(
  Data              = soil_water_3,
  Covariates_coarse = covs[[1]],
  Covariates_fine   = covs[[2]],
  KrigingEquation   = Data ~ layer,
  Keep_Temporary    = FALSE,
  Cores             = 4L,
  Dir               = "ipms/Data/weather/krigr",
  FileName          = "soil_water_3_kriged.nc",
  nmax              = 15
)

sw3_krig_pts <- extract_krig(soil_water_3_krig_data, 
                             coords, 
                             "ID",
                             mnts)

tim <- round(Sys.time() - start, 3)
print(tim)

saveRDS(soil_water_3_krig_data, file = "ipms/Data/weather/krigr/kriged_sw3_data.rds")
saveRDS(sw3_krig_pts, file = "ipms/Data/weather/krigr/kriged_sw3_pts.rds")
# source("ipms/R/utils/utils-packages2.R")