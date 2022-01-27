

# Find observations w/ pareto-K's above a certain threshold value from LOO object
find_prob_obs <- function(loo_obj, thresh) {
  
  return(which(loo_obj$diagnostics$pareto_k > thresh))
  
}

# Get exact row IDs for a data set who's pareto-k values are above a
# certain threshold
get_prob_rows <- function(model_obj, loo_obj, thresh) {
  
  stopifnot(length(model_obj) == length(loo_obj))
  
  # index of values for the model_obj@frame slot
  data_ind <- lapply(loo_obj, FUN = find_prob_obs, thresh = thresh)
  
  
  out <- lapply(seq_len(length(model_obj)), function(x) {
    model_obj[[x]]$data[data_ind[[x]], ]
  })
  
  return(out)
  
}

sourceCpp(file = 'ipms/R/Cpp/cpp_utils.cpp')

# Extras ------

# Correct data types from QGIS polygons for use in C++ merge_ramets/calc_surv


check_outliers <- function(all_ramets, population, type) {
  
  switch(type,
         "growth"  = check_growth_outliers(all_ramets, population),
         "big_new" = check_big_new_outerliers(all_ramets, population))
  
}

# Model fitting and checking ---------------

fit_vr_model <- function(data, vr, clim, native = "no") {
  
  switch(native,
         "yes" = .fit_native_model(data, vr),
         "no"  = .fit_clim_model(data, vr, clim))
  
}

.fit_native_model <- function(data, vr) {
  
  fam <- switch(
    vr,
    "repro"     = bernoulli(),
    "flower_n"  = negbinomial(),
    "alive"     = bernoulli(),
    "size_next" = gaussian()
  )
  
  form_1 <- as.formula(glue("{vr} ~ log_size + (1 | site)"))
  form_2 <- as.formula(glue("{vr} ~ log_size + native + (1 | site)"))
  form_3 <- as.formula(glue("{vr} ~ log_size * native + (1 | site)"))
  
  be <- "cmdstanr"
  
  if(vr %in% c("flower_n")) {
    inits <- "0"
    
  } else {
    
    inits <- "random"
  }
  
  if(vr == "alive") vr <- "survival"
  
  # NB: Sigma has log link by default! Remember to exponentiate at IPM
  # build time.
  if(vr == "log_size_next") {
    
    sigma_form_1 <- sigma ~ log_size
    sigma_form_2 <- sigma ~ log_size + native
    sigma_form_3 <- sigma ~ log_size * native
    
    
    form_1 <- bf(form_1, sigma_form_1)
    form_2 <- bf(form_2, sigma_form_2)
    form_3 <- bf(form_3, sigma_form_3)
    
    vr <- "growth"
  }
  
  
  message("Model 1: size only: Native----------\n\n")
  mod_1 <- brm(form_1,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_size_only_native.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 2: Size Plus Native-------\n\n")
  
  mod_2 <- brm(form_2,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_size_plus_native.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 3: Size Times Native------\n\n")
  
  mod_3 <- brm(form_3,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_size_times_native.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  mod_waic <- waic(mod_1, mod_2, mod_3)
  
  out <- list(size_only         = mod_1,
              size_plus_native  = mod_2,
              size_times_native = mod_3,
              waic              = mod_waic)
  
  return(out)
}


.fit_clim_model <- function(data, vr, clim) {
  
  fam <- switch(
    vr,
    "repro"     = bernoulli(),
    "flower_n"  = negbinomial(),
    "alive"     = bernoulli(),
    "size_next" = gaussian()
  )
  
  form_1 <- as.formula(glue("{vr} ~ log_size + (1|site)"))
  form_2 <- as.formula(glue("{vr} ~ log_size + {clim} + (1|site)"))
  form_3 <- as.formula(glue("{vr} ~ log_size * {clim} + (1|site)"))
  
  be <- "cmdstanr"
  
  if(vr %in% c("flower_n")) {
    inits <- "0"
    
  } else {
    
    inits <- "random"
  }
  
  if(vr == "alive") vr <- "survival"
  
  # NB: Sigma has log link by default! Remember to exponentiate at IPM
  # build time.
  if(vr == "log_size_next") {
    
    sigma_form_1 <- sigma ~ log_size
    sigma_form_2 <- as.formula(glue("sigma ~ log_size + {clim}"))
    sigma_form_3 <- as.formula(glue("sigma ~ log_size * {clim}"))
    
    
    form_1 <- bf(form_1, sigma_form_1)
    form_2 <- bf(form_2, sigma_form_2)
    form_3 <- bf(form_3, sigma_form_3)
    
    vr <- "growth"
  }
  
  
  message("Model 1: size only----------\n\n")
  mod_1 <- brm(form_1,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_int_only.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 2: Size Plus climate-------\n\n")
  
  mod_2 <- brm(form_2,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_size_plus_clim.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  message("\n\nModel 3: Size Times Climate------\n\n")
  
  mod_3 <- brm(form_3,
               data       = data,
               family     = fam,
               chains     = 4,    
               backend    = be,
               inits      = inits,
               cores      = getOption("mc.cores", 4L),
               save_model = glue('ipms/Stan/{vr}_{clim}_size_times_clim.stan'),
               control    = list(adapt_delta = 0.99,
                                 max_treedepth = 15))
  
  mod_waic <- waic(mod_1, mod_2, mod_3)
  
  out <- list(size_only       = mod_1,
              size_plus_clim  = mod_2,
              size_times_clim = mod_3,
              waic            = mod_waic)
  
  return(out)
  
  
}

