# Download CHELSA Monthly data from chelsa website, store TIFFs locally, and then
# create file that is actually useful for Carpobrotus analysis


# Downloads ----
mnts <- c(paste0("0", 1:9), 10:12)
yrs  <- 2017:2019

mnts_years <- expand.grid(mnts, yrs) %>% 
  mutate(mnt_yr = paste(.[,1], .[,2], sep = "_")) %>%
  select(mnt_yr) %>% 
  unlist(use.names = FALSE) %>%
  paste("_V.2.1.tif", sep = "")

clim_vars <- c(
  # "tasmax", "tasmin",
  "pr", "pet"
)

chelsa_base <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly"


for(i in seq_along(clim_vars)) {
  
  clim <- clim_vars[i]
  
  if(clim == "pet"){
    
    clim_1 <- gsub("pet", "pet_penman", clim)
    
    base_pat <- glue("{clim}/CHELSA_{clim_1}")
    
  } else {
    base_pat <- glue("{clim}/CHELSA_{clim}")
    
  }
  
  for(j in seq_along(mnts_years)) {
    
    message("\n", "Downloading ", clim, " for ", mnts_years[j], "\n")
    
    
    use_pat <- glue("{chelsa_base}/{base_pat}_{mnts_years[j]}")
    
    out_pat <- gsub("_V\\.2\\.1\\.tif$", "\\.tif", mnts_years[j])
    
    out_pat <- glue("ipms/Data/weather/chelsa_{clim}_{out_pat}")
    
    download.file(use_pat,
                  destfile = out_pat,
                  method = "auto")
  }
  
}

