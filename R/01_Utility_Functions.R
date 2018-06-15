# Utilities for carpobrotus IPM pipeline

infile_data <- function(path, name, log.breaks) {
  
  expr <- expression(st_read(paste(getwd(), path, sep = '/')) %>%
                       st_buffer(0) %>% # splits mistakenly intersected polygons
                       mutate(LogSize = log(Size)) %>%
                       mutate(SizeBin = cut(.$LogSize, breaks = log.breaks)) %>%
                       mutate(CleanBin = clean_bins(SizeBin)) %>%
                       select(-SizeBin)
  )
  
  assign(name, eval(expr), envir = .GlobalEnv)
  
}


ramet_to_genet <- function(x, ..., size.breaks) {
  group_var <- quos(...)
  
  x %>%
    group_by(!!! group_var) %>%
    summarise(N_Ramets = length(id),
              Size = sum(Size),
              Flower_N = sum(Flower_N, na.rm = TRUE)) %>%
    mutate(LogSize = log(Size)) %>% 
    mutate(SizeBin = cut(.$LogSize, breaks = size.breaks)) %>%
    mutate(CleanBin = clean_bins(SizeBin)) %>%
    select(-SizeBin)
}

clean_bins <- function(x) {
  lapply(x,
         function(x) str_split(x, '\\(|,|\\]')[[1]][3] %>%
           as.numeric()) %>% 
    unlist()
}


get_gg_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
