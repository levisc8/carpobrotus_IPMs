# Climate utils 

# This function is no longer used in the IPM pipeline but I love it too much to
# delete. 

giovanni_dl <- function(link_file, out_dir, user_name = "levisc8") {
  
  links <- as.character(readLines(link_file)) %>%
    Filter(f = function(x) x[[1]] != "", x = .)
  
  if(length(links) == 1 && links == "NULL" || rlang::is_empty(links)) {
    
    return(TRUE)
    
  } else if(any(grepl("^https", links))) {
    
    system2(
      "wget",
      args = glue("
          --load-cookies .urs_cookies --save-cookies .urs_cookies --keep-session-cookies --auth-no-challenge=on --user={user_name} --password={getOption('gio_pw')} --content-disposition -i {link_file} -P {out_dir}
    ")
    )
    
    downloaded_files <- fs::dir_ls(out_dir) %>%
      str_remove(pattern = paste0(out_dir,"/"))
    
    poss_downloads <- stringr::str_split(links, "/") %>%
      Filter(f = function(x) !x[[1]] == "", x = .) %>%
      lapply(function(x) x[length(x)]) %>% 
      unlist() 
    
    missing_files <- setdiff(poss_downloads, downloaded_files)
    
    if(rlang::is_empty(missing_files)) {
      
      message("Downloads complete.")
      
      sink(file = link_file, append = FALSE)
        cat("NULL", sep = "\n\n")
      sink()      
      
      return(TRUE)
    }
    
    missing_links <- stringr::str_extract_all(links, 
                                              paste(missing_files, 
                                                    collapse = "|")) %>%
      lapply(function(x) !rlang::is_empty(x)) %>%
      unlist()
    
    missing_links <- links[missing_links]
    
    sink(file = link_file, append = FALSE)
      cat(missing_links, sep = "\n\n")
    sink()
    
    giovanni_dl(link_file, out_dir, user_name)
    
  } else {
    stop("HTTPS links not detected in 'link_file'.",
         " Are you sure these are from GIOVANNI?")
  }

}


# Unit conversions -----------

# Precip is kg/m^2/second in FLDAS. 1 Kg of water will cover 1 square meter
# in 1mm of water e.g. 1 kg/m2/s = 1 mm/s. so 60 sec/m * 60 min/h * 24h/d *
# 30d/month yields mm/month

kg_m_s_to_mm_mnt <- function(values, n_days = 30) {
  
  values * 60 * 60 * 24 * n_days
  
}

k_to_c <- function(values) {
  values - 273.15
}


# 
extract_krig <- function(krig_out, coords, id_col, mnts) {
  
  out <- list()
  
  out$data <- setNames(
    as.data.frame(
      raster::extract(
        krig_out$Kriging_Output, 
        coords[ , 1:2])
    ),
    mnts
  )
  
  out$data$ID <- coords[ , id_col]
  
  out$se <- setNames(
    as.data.frame(
      raster::extract(
        krig_out$Kriging_SE, 
        coords[ , 1:2])
    ),
    paste(mnts, "SE", sep = "_")
  )
  out$se$ID <- coords[ , id_col]
  
  return(out)
}


tidy_clim_data <- function(dat, clim_var) {
  
  dat <- dat$data
  ses <- dat$se
  
  mt_yr_df <- expand.grid(month.name, 2017:2020, stringsAsFactors = FALSE) %>%
    slice(-c(1:3, 40:48)) %>% 
    mutate(dat_name = paste(Var1, Var2, sep = '_'))
  
  id_df <- select(dat, ID:Lat)
  
  # Make long, then split to wet/dry seasons
  long_dat <- pivot_longer(dat, -(ID:Lat), names_to = "month_year", values_to = "val") 
  
  long_dat$year <- vapply(long_dat$month_year, function(x){
    strsplit(x, "_")[[1]][2]
  }, character(1L))
  
  long_dat$month <- vapply(long_dat$month_year, function(x){
    strsplit(x, "_")[[1]][1]
  }, character(1L))
  
  seas_dat <- mutate(long_dat, season = case_when(
    Lat > 0 & month %in% c( "May", "June", "July", 
                           "August", "September", "October") ~ "dry",
    Lat > 0 & month %in% c("November", "December", "January",
                           "February", "March", "April") ~ "wet",
    Lat < 0 & month %in% c("April", "May", "June", "July", 
                           "August", "September") ~ "wet",
    TRUE ~ "dry"
  ))
  
  agg_dat <- .aggregate_data(seas_dat, clim_var) %>%
    lapply(function(x, ids) left_join(x, ids), ids = id_df)
  
  return(agg_dat)
}


.aggregate_data <- function(long_dat, clim_var) {
  
  if(clim_var == "prec") {
    f <- sum
    nm_1 <- sym("total_prec")
    nm_2 <- sym("seas_prec")
  } else {
    f = mean
    nm_1 <- sym(glue("mean_{clim_var}"))
    nm_2 <- sym(glue("seas_{clim_var}"))
  }
  
  seas_dat <- long_dat %>%
    group_by(ID, season, year) %>%
    summarise(
      !! nm_1 := f(val, na.rm = TRUE)
    ) %>% 
    ungroup() %>%
    pivot_wider(id_cols = ID,
                names_from = c(season, year),
                names_glue = "{season}_{year}",
                values_from = 4) %>%
    map_dfc(.f = function(x) { 
      if(is.numeric(x)) {
        as.vector(scale(x))
      } else {
        x
      }
    }) %>%
    setNames(c("ID", paste(clim_var, names(.)[2:ncol(.)], sep = "_"))) %>% 
    remove_outliers()
  
  ann_dat <- long_dat %>%
    group_by(ID,year) %>%
    summarise(
      !! nm_1 := f(val, na.rm = TRUE),
      !! nm_2 := sd(val) / mean(val, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    pivot_wider(id_cols = ID,
                names_from = c(year),
                values_from = 3:4) %>%
    map_dfc(.f = function(x) { 
      if(is.numeric(x)) {
        as.vector(scale(x))
      } else {
        x
      }
    }) %>% 
    remove_outliers()
  
  list(seasonal = seas_dat,
       annual   = ann_dat)
}

# Remove huge outliers in some of the Kriged values, but retains all of our
# field sites (we have to have those values).

remove_outliers <- function(dat) {
  
  site_ids <- paste0("p", 9198:9221)
  site_dat <- filter(dat, ID %in% site_ids)
  rng_dat  <- filter(dat, !ID %in% site_ids)
  
  map_dfc(rng_dat,
          function(x) {
            if(is.numeric(x)) {
              .rm_outs(x)
            } else {
              x
            }
          }) %>%
    rbind(site_dat)
  
}

.rm_outs <- function(x, na.rm = TRUE, ...) {
  
  qnt <- quantile(x, probs=c(.025, .975), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
