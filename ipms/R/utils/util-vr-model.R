
# Correct data types from QGIS polygons for use in C++ merge_ramets/calc_surv

check_outliers <- function(all_ramets, population, type) {
  
  switch(type,
         "growth"  = check_growth_outliers(all_ramets, population),
         "big_new" = check_big_new_outerliers(all_ramets, population))
  
}

# Model fitting and checking ---------------

fit_vr_model <- function(data, vr, gam) {

  .fit_native_model(data, vr, gam)
  
}

.fit_native_model <- function(data, vr, gam) {
  
  stan_gam <- ifelse(gam, "gam", "lin")
  
  fam <- switch(
    vr,
    "repro"     = bernoulli(),
    "flower_n"  = negbinomial(),
    "alive"     = bernoulli(),
    "log_size_next" = gaussian()
  )
  
  seas_forms <- switch(stan_gam,
                       "gam" = .make_bfs(vr, "seas"),
                       "lin" = .make_lin_bfs(vr, "seas"))
  ann_forms  <- switch(stan_gam,
                       "gam" = .make_bfs(vr, "ann"),
                       "lin" = .make_lin_bfs(vr, "ann"))
  
  # Drop duplicate native-only formula
  all_forms <- c(seas_forms, ann_forms[2:length(ann_forms)])
  
  be <- "cmdstanr"
  
  if(vr %in% c("flower_n")) {
    
    inits <- "0"
    
  } else {
    
    inits <- "random"
  }
  
  if(vr == "alive") vr <- "survival"
  
  if(vr == "log_size_next") vr <- "growth"
  
  
  
  out <- list()
  
  for(i in seq_along(all_forms)) {
    
    rds_fp <- glue("ipms/Model_Fits/{stan_gam}/{vr}/{names(all_forms)[i]}.rds")
    
    if(file_exists(rds_fp)) {
      
      out[[i]] <- readRDS(rds_fp)
      names(out)[i] <- names(all_forms)[i]
      message(glue("Skipping Model {i} for {vr}: already modeled :))"))
      next
      
    }
    
    # Save gam stan code
    if(i == 2) {
      mod_fp <- glue('ipms/Stan/{vr}_{stan_gam}_seas_krig_data.stan')
    } else if (i == 8) {
      mod_fp <- glue('ipms/Stan/{vr}_{stan_gam}_ann_krig_data.stan')
    } else {
      mod_fp <- NULL
    }
    message(glue("Model {i} for {vr}: {names(all_forms)[i]}----------\n\n"))
    
    out[[i]] <- brm(all_forms[[i]],
                    data       = data,
                    family     = fam,
                    chains     = 4,    
                    backend    = be,
                    init       = inits,
                    cores      = getOption("mc.cores", 4L),
                    save_model = mod_fp,
                    control    = list(adapt_delta = 0.999,
                                      max_treedepth = 15))
    names(out)[i] <- names(all_forms)[i]
    
    saveRDS(out[[i]],
            file = rds_fp)
    
  }
 
  mod_waic <- waic(out[[1]],
                   out[[2]],
                   out[[3]],
                   out[[4]],
                   out[[5]],
                   out[[6]],
                   out[[7]],
                   out[[8]],
                   out[[9]],
                   out[[10]],
                   out[[11]],
                   out[[12]],
                   out[[13]],
                   model_names = names(out))
  
  out <- c(out, mod_waic = mod_waic)
  
  return(out)
}


.make_bfs <- function(vr, clim) {
  
  seas_clim_vars <- c("temp_dry", "temp_wet", 
                      "prec_dry", "prec_wet")
  ann_clim_vars  <- c("mean_temp", "seas_temp",
                      "total_prec", "seas_prec")
  
  seas_sw1_vars  <- c(seas_clim_vars,  
                      "sw1_dry","sw1_wet")
  
  ann_sw1_vars   <-  c(ann_clim_vars,
                       "mean_sw1", "seas_sw1")
  
  seas_sw2_vars  <- c(seas_clim_vars,  
                      "sw2_dry","sw2_wet")
  
  ann_sw2_vars   <-  c(ann_clim_vars,
                       "mean_sw2", "seas_sw2")
  
  seas_sw3_vars  <- c(seas_clim_vars,  
                      "sw3_dry","sw3_wet")
  
  ann_sw3_vars   <-  c(ann_clim_vars,
                       "mean_sw3", "seas_sw3")
  
  # *_t_1 = *_(t minus 1)
  # *_t   = *_(time of sample start) repro + flower_n are lagged by a year
  # because we use vr(t) ~ size (t) to estimate the relationship, and so weather
  # that would matter to it is from the previous year. growth/survival are
  # conditional on the weather from the transition itself, as those are
  # evaluated at t+1.
  
  pc_1_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw1_vars, "_t_1"),
      "ann"  = paste0(ann_sw1_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw1_vars, "_t"),
      "ann"  = paste0(ann_sw1_vars, "_t")
    )
  )
  
  pc_1_add_term <- paste(
    paste0(
      "s(", pc_1_terms, ", k = 4, bs = 'cs')"), 
    collapse = " + ")
  pc_1_int_term <- paste(
      paste0(
        "t2(", pc_1_terms, ", log_size, k = 4, bs = 'cs')"
      ),
    collapse = " + "
  )
  
  pc_2_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw2_vars, "_t_1"),
      "ann"  = paste0(ann_sw2_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw2_vars, "_t"),
      "ann"  = paste0(ann_sw2_vars, "_t")
    )
  )
  
  pc_2_add_term <- paste(
    paste0(
      "s(", pc_2_terms, ", k = 4, bs = 'cs')"
      ),
    collapse = " + ")
  pc_2_int_term <- paste(
    paste0(
      "t2(", pc_2_terms, ", log_size, k = 4, bs = 'cs')"
    ),   
    collapse = " + "
  )
  
  pc_3_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw3_vars, "_t_1"),
      "ann"  = paste0(ann_sw3_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw3_vars, "_t"),
      "ann"  = paste0(ann_sw3_vars, "_t")
    )
  )
  
  pc_3_add_term <- paste(
    paste0(
      "s(", pc_3_terms, ", k = 4, bs = 'cs')"
    ),
    collapse = " + ")
  pc_3_int_term <- paste(
    paste0(
      "t2(", pc_3_terms, ", log_size, k = 4, bs = 'cs')"
    ),    collapse = " + "
  )
  
  native_only <- as.formula(glue("{vr} ~ log_size + native + (1 | site)"))
  add_sw1     <- as.formula(glue("{vr} ~ log_size + {pc_1_add_term} + native + (1 | site)"))
  times_sw1   <- as.formula(glue("{vr} ~ log_size + {pc_1_int_term} + native + (1 | site)"))
  add_sw2     <- as.formula(glue("{vr} ~ log_size + {pc_2_add_term} + native + (1 | site)"))
  times_sw2   <- as.formula(glue("{vr} ~ log_size + {pc_2_int_term} + native + (1 | site)"))
  add_sw3     <- as.formula(glue("{vr} ~ log_size + {pc_3_add_term} + native + (1 | site)"))
  times_sw3   <- as.formula(glue("{vr} ~ log_size + {pc_3_int_term} + native + (1 | site)"))
  
  if(vr == "log_size_next") {
    
    sig_1_add_term <- paste(pc_1_terms, collapse = " + ")
    sig_2_add_term <- paste(pc_2_terms, collapse = " + ")
    sig_3_add_term <- paste(pc_3_terms, collapse = " + ")
    

    sigma_form <- as.formula(
      glue(
        "sigma ~ log_size + native + (1 | site)"
      )
    )
    
    sigma_1_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_1_add_term} + native + (1 | site)"
      )
    )
    
    sigma_2_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_2_add_term} + native + (1 | site)"
      )
    )
    
    sigma_3_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_3_add_term} + native + (1 | site)"
      )
    )
    
    
    native_only <- bf(native_only, sigma_form)
    add_sw1     <- bf(add_sw1, sigma_1_form)
    times_sw1   <- bf(times_sw1, sigma_1_form)
    add_sw2     <- bf(add_sw2, sigma_2_form)
    times_sw2   <- bf(times_sw2, sigma_2_form)
    add_sw3     <- bf(add_sw3, sigma_3_form)
    times_sw3   <- bf(times_sw3, sigma_3_form)
    
  }
  
  base_nms <- c("add_sw1", "times_sw1", 
                "add_sw2", "times_sw2",
                "add_Sw3", "times_sw3")
  
  out_nms <- c("native_only", glue("{base_nms}_{clim}"))
  
  out <- list(native_only,
              add_sw1, times_sw1,
              add_sw2, times_sw2, 
              add_sw3, times_sw3) %>%
    setNames(out_nms)
  
  return(out)
  
  
  
}

.make_lin_bfs <- function(vr, clim) {
  
  seas_clim_vars <- c("temp_dry", "temp_wet", 
                      "prec_dry", "prec_wet")
  ann_clim_vars  <- c("mean_temp", "seas_temp",
                      "total_prec", "seas_prec")
  
  seas_sw1_vars  <- c(seas_clim_vars,  
                      "sw1_dry","sw1_wet")
  
  ann_sw1_vars   <-  c(ann_clim_vars,
                       "mean_sw1", "seas_sw1")
  
  seas_sw2_vars  <- c(seas_clim_vars,  
                      "sw2_dry","sw2_wet")
  
  ann_sw2_vars   <-  c(ann_clim_vars,
                       "mean_sw2", "seas_sw2")
  
  seas_sw3_vars  <- c(seas_clim_vars,  
                      "sw3_dry","sw3_wet")
  
  ann_sw3_vars   <-  c(ann_clim_vars,
                       "mean_sw3", "seas_sw3")
  
  # *_t_1 = *_(t minus 1)
  # *_t   = *_(time of sample start) repro + flower_n are lagged by a year
  # because we use vr(t) ~ size (t) to estimate the relationship, and so weather
  # that would matter to it is from the previous year. growth/survival are
  # conditional on the weather from the transition itself, as those are
  # evaluated at t+1.
  
  pc_1_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw1_vars, "_t_1"),
      "ann"  = paste0(ann_sw1_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw1_vars, "_t"),
      "ann"  = paste0(ann_sw1_vars, "_t")
    )
  )
  
  pc_1_add_term <- paste(
    pc_1_terms, 
    collapse = " + "
  )
  pc_1_int_term <- paste(
    paste0( pc_1_terms, " * log_size"
    ),
    collapse = " + "
  )
  
  pc_2_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw2_vars, "_t_1"),
      "ann"  = paste0(ann_sw2_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw2_vars, "_t"),
      "ann"  = paste0(ann_sw2_vars, "_t")
    )
  )
  
  pc_2_add_term <- paste(
    pc_2_terms,
    collapse = " + ")
  pc_2_int_term <- paste(
    paste0(pc_2_terms, " * log_size"
    ),   
    collapse = " + "
  )
  
  pc_3_terms <- switch(
    vr,
    "repro" = ,
    "flower_n" = switch(
      clim,
      "seas" = paste0(seas_sw3_vars, "_t_1"),
      "ann"  = paste0(ann_sw3_vars, "_t_1")
    ),
    "alive" = ,
    "log_size_next" = switch(
      clim,
      "seas" = paste0(seas_sw3_vars, "_t"),
      "ann"  = paste0(ann_sw3_vars, "_t")
    )
  )
  
  pc_3_add_term <- paste(
    pc_3_terms, 
    collapse = " + ")
  pc_3_int_term <- paste(
    paste0(pc_3_terms, " * log_size"),
    collapse = " + "
  )
  
  native_only <- as.formula(glue("{vr} ~ log_size + native + (1 | site)"))
  add_sw1     <- as.formula(glue("{vr} ~ log_size + {pc_1_add_term} + native + (1 | site)"))
  times_sw1   <- as.formula(glue("{vr} ~ log_size + {pc_1_int_term} + native + (1 | site)"))
  add_sw2     <- as.formula(glue("{vr} ~ log_size + {pc_2_add_term} + native + (1 | site)"))
  times_sw2   <- as.formula(glue("{vr} ~ log_size + {pc_2_int_term} + native + (1 | site)"))
  add_sw3     <- as.formula(glue("{vr} ~ log_size + {pc_3_add_term} + native + (1 | site)"))
  times_sw3   <- as.formula(glue("{vr} ~ log_size + {pc_3_int_term} + native + (1 | site)"))
  
  if(vr == "log_size_next") {
    
    sig_1_add_term <- paste(pc_1_terms, collapse = " + ")
    sig_2_add_term <- paste(pc_2_terms, collapse = " + ")
    sig_3_add_term <- paste(pc_3_terms, collapse = " + ")
    
    
    sigma_form <- as.formula(
      glue(
        "sigma ~ log_size + native + (1 | site)"
      )
    )
    
    sigma_1_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_1_add_term} + native + (1 | site)"
      )
    )
    
    sigma_2_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_2_add_term} + native + (1 | site)"
      )
    )
    
    sigma_3_form <- as.formula(
      glue(
        "sigma ~ log_size + {sig_3_add_term} + native + (1 | site)"
      )
    )
    
    
    native_only <- bf(native_only, sigma_form)
    add_sw1     <- bf(add_sw1, sigma_1_form)
    times_sw1   <- bf(times_sw1, sigma_1_form)
    add_sw2     <- bf(add_sw2, sigma_2_form)
    times_sw2   <- bf(times_sw2, sigma_2_form)
    add_sw3     <- bf(add_sw3, sigma_3_form)
    times_sw3   <- bf(times_sw3, sigma_3_form)
    
  }
  
  base_nms <- c("add_sw1", "times_sw1", 
                "add_sw2", "times_sw2",
                "add_Sw3", "times_sw3")
  
  out_nms <- c("native_only", glue("{base_nms}_{clim}"))
  
  out <- list(native_only,
              add_sw1, times_sw1,
              add_sw2, times_sw2, 
              add_sw3, times_sw3) %>%
    setNames(out_nms)
  
  return(out)
  
}

