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
