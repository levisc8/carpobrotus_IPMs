# rename fixed effects from brms to ipmr names

name_fixed_pars <- function(model, vr, sw = 3) {
  
  mod_nm <- deparse(substitute(model))
  
  base_nms <- paste(vr, 
                    c(
                      "int", "z",
                      "temp_dry", "temp_wet", 
                      "prec_dry", "prec_wet",
                      paste(paste0("sw", sw), 
                            c("dry", "wet"), sep = "_"),
                      "native"
                    ), sep = "_"
  )
  
  interaction_nms <- rownames(fixef(model)) %>%
    .[grepl("\\:", .)] %>%
    gsub(pattern = "\\:", replacement = "_x_", x = .) %>%
    gsub(pattern = "_t$|_t_1$", replacement = "", x = .)
  interaction_nms <- paste(vr, interaction_nms, sep = "_")
  
  out <- c(base_nms, interaction_nms)
  
  if(mod_nm == "grow_mod") {
    
    base_nms <- c(base_nms[1],
                  paste0(vr, "_sigma_int"), 
                  base_nms[2:length(base_nms)])
    
    sig_nms <- rownames(fixef(model)) %>%
      .[grepl("sigma_", .) & !grepl("Intercept", .)] %>%
      gsub(pattern = "_t$$", replacement = "", x = .)
    sig_nms <- paste0("grow_", sig_nms)
    
    out <- c(base_nms, interaction_nms, sig_nms)
  }
  
  out <- setNames(fixef(model)[ , 1], out)
  
  return(as.list(out))
  
}

# rename random effects from brms to ipmr names

name_ran_pars <- function(model, vr) {
  
  mod_nm <- deparse(substitute(model))
  
  pops <- rownames(ranef(model)[[1]][ , , 1])
  
  base_nm <- paste(vr, "intercept", pops, sep = "_")
  
  out <- setNames(ranef(model)[[1]][ , 1, 1], base_nm)
  
  if(mod_nm == "grow_mod") {
    sig_nms <- rownames(ranef(model)[[1]][ , , 2])
    sig_nms <- paste0("grow_sigma_intercept_", sig_nms)
    sig_out <- setNames(ranef(model)[[1]][ , 1, 2], sig_nms)
    
    out <- c(out, sig_out)
  }
  
  return(as.list(out))
  
}

# Rename posterior draws data.frames to ipmr

rename_draws <- function(draws, vr) {
  
  nm <- deparse(substitute(draws))
  
  # Fix random effects names
  names(draws)[grepl("\\,Intercept\\]", names(draws))] <- 
    gsub("\\,Intercept\\]", 
         "", 
         names(draws)[grepl("\\,Intercept\\]", names(draws))])
  
  names(draws)[grepl("r_site\\[", names(draws))] <-
    gsub("r_site\\[",
         paste0(vr, "_intercept_"), 
         names(draws)[grepl("r_site\\[", names(draws))])
  
  # Site random effect for growth variance model
  names(draws)[grepl("r_site__sigma\\[", names(draws))] <-
    gsub("r_site__sigma\\[",
         paste0(vr, "_sigma_intercept_"), 
         names(draws)[grepl("r_site", names(draws))])
  
  # Fix the fixed effects now.
  names(draws)[!grepl(vr, names(draws))] <- 
    paste(vr,
          names(draws)[!grepl(vr, names(draws))], sep = "_")
  
  # Remove leading 'b's (now no longer leading)
  
  names(draws) <- gsub("_b_", "_", names(draws))
  
  # Interaction terms
  names(draws) <- gsub("\\:", "_x_", names(draws))
  
  # climate terms
  names(draws) <- gsub("_t$|_t_1$", "", names(draws))
  
  # Fix intercepts
  
  names(draws) <- gsub("Intercept", "int", names(draws))
  # Fix the name of the size coefficient
  z_nm <- paste0(vr, "_log_size")
  
  names(draws)[names(draws) == z_nm] <- paste0(vr, "_z")
  
  
  return(draws)
}

# Simple switch function that takes a symbol site name and returns native status

is_native <- function(x) {
  x <- deparse(substitute(x))
  if(x %in% c("Struisbaai", "Rooisand",
              "Vogelgat", "Springfontein",
              "Melkboss", "St_Francis")) {
    return(1)
  } else {
    return(0)
  }
}

# Returns predicted values w/o CIs from brms

predict_carp_ipm <- function(model, data) {
  x <- predict(model, newdata = data, cores = 4L)
  unname(x[ , 1, drop = TRUE])
}

# Either summarizes the draws to mean posterior value, or selects a specific
# draw (e.g. for uncertainty simulation)
prep_draws <- function(par, draws, draw_id, f = mean) {
  
  if(!is.null(draw_id)){
    
    draws[grepl(par, names(draws))] %>%
      .[draw_id, ] %>%
      unlist()
    
  } else {
    
    draws[grepl(par, names(draws))] %>%
      apply(2, f)
    
  }
}

# Functions that vectorize prediction from GAMs

pred_surv_gam <- function(mod, z, clim_data, site,
                          draw_id = NULL, sum_fun = mean) {
  
  site <- deparse(substitute(site))
  
  log_size <- expand.grid(log_size = z,
                          site = unique(clim_data$site)) %>%
    mutate(
      native = case_when(
        site %in% c("Melkboss",      "Vogelgat", 
                    "St_Francis",    "Struisbaai",
                    "Springfontein", "Rooisand")   ~ 1,
        TRUE ~ 0
      ),
      alive = 0.25 # Create dummy variable for survival so standata doesn't complain
    )  %>%
    left_join(clim_data, by = "site")
  
  ind <- which(log_size$site == site)
  
  data <- make_standata(mod$formula$formula,
                        log_size,
                        family = bernoulli())
  code <- make_stancode(mod$formula$formula,
                        log_size,
                        family = bernoulli(),
                        save_model = "ipms/Stan/surv_times_sw2_seas_gam.stan")
  
  # Initialize the spline coefficients and the basis matrices
  draws <- as.data.frame(mod)
  Z_1 <- data$Zs_1_1[ind, ]
  S_1 <- prep_draws("s_t2mean_temp_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_2 <- data$Zs_2_1[ind, ]
  S_2 <- prep_draws("s_t2seas_temp_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_3 <- data$Zs_3_1[ind, ]
  S_3 <- prep_draws("s_t2total_prec_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_4 <- data$Zs_4_1[ind, ]
  S_4 <- prep_draws("s_t2seas_prec_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_5 <- data$Zs_5_1[ind, ]
  S_5 <- prep_draws("s_t2mean_sw2_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_6 <- data$Zs_6_1[ind, ]
  S_6 <- prep_draws("s_t2seas_sw2_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  # Fixed Effects
  globe_int <- prep_draws("b_Intercept", draws, draw_id, sum_fun)
  Bs        <- c(prep_draws("b_log_size", draws, draw_id, sum_fun),
                 prep_draws("b_native", draws, draw_id, sum_fun))
  X         <- data$X[ind , 2:3]
  site_int  <- ifelse(is.null(draw_id), 
                      mean(draws[ ,grepl(site, names(draws))]),
                      draws[draw_id, grepl(site, names(draws))])
  
  
  lin_p <- globe_int +
    (Z_1 %*% S_1) + 
    (Z_2 %*% S_2) + 
    (Z_3 %*% S_3) + 
    (Z_4 %*% S_4) + 
    (Z_5 %*% S_5) +
    (Z_6 %*% S_6) +
    (X %*% Bs) +
    site_int
  
  as.vector(inv_logit(lin_p))
  
}

pred_repr_gam <- function(mod, z, clim_data, site,
                          draw_id = NULL, sum_fun = mean) {
  
  site <- deparse(substitute(site))
  
  log_size <- expand.grid(log_size = z,
                          site = unique(clim_data$site)) %>%
    mutate(
      native = case_when(
        site %in% c("Melkboss",      "Vogelgat", 
                    "St_Francis",    "Struisbaai",
                    "Springfontein", "Rooisand")   ~ 1,
        TRUE ~ 0
      ),
      repro = 0.25 # Create dummy variable for survival so standata doesn't complain
    )  %>%
    left_join(clim_data, by = "site")
  
  ind <- which(log_size$site == site)
  
  data <- make_standata(mod$formula$formula,
                        log_size,
                        family = bernoulli())
  
  code <- make_stancode(mod$formula$formula,
                        log_size,
                        family = bernoulli(),
                        save_model = "ipms/Stan/repr_mix_lin_gam.stan")
  
  # Initialize the spline coefficients and the basis matrices
  draws <- as.data.frame(mod)
  Z_1 <- data$Zs_1_1[ind, ]
  S_1 <- prep_draws("s_smean_temp_t_1_1\\[", draws, draw_id, sum_fun)
  
  Z_2 <- data$Zs_2_1[ind, ]
  S_2 <- prep_draws("s_stotal_prec_t_1_1\\[", draws, draw_id, sum_fun)
  
  globe_int <- prep_draws("^b_Intercept$", draws, draw_id, sum_fun)
  Bs        <- c(prep_draws("^b_log_size$", draws, draw_id, sum_fun),
                 prep_draws("^b_seas_temp_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_seas_prec_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_mean_sw2_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_seas_sw2_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:seas_temp_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:seas_prec_t_1$", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:mean_sw2_t_1$" , draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:seas_sw2_t_1$", draws, draw_id, sum_fun))
  
  X         <- data$X[ind , -1]
  site_int  <- ifelse(is.null(draw_id), 
                      mean(draws[ ,grepl(site, names(draws))]),
                      draws[draw_id, grepl(site, names(draws))])
  
  lin_p <- globe_int +
    (Z_1 %*% S_1) + 
    (Z_2 %*% S_2) + 
    (X %*% Bs) +
    site_int
  
  as.vector(inv_logit(lin_p))
}

pred_flow_gam <- function(mod, z, clim_data, site,
                          draw_id = NULL, sum_fun = mean) {
  
  site <- deparse(substitute(site))
  
  log_size <- expand.grid(log_size = z,
                          site = unique(clim_data$site)) %>%
    mutate(
      native = case_when(
        site %in% c("Melkboss",      "Vogelgat", 
                    "St_Francis",    "Struisbaai",
                    "Springfontein", "Rooisand")   ~ 1,
        TRUE ~ 0
      ),
      flower_n = 5 # Create dummy variable for survival so standata doesn't complain
    )  %>%
    left_join(clim_data, by = "site")
  
  ind <- which(log_size$site == site)
  
  data <- make_standata(mod$formula$formula,
                        log_size,
                        family = negbinomial())
  code <- make_stancode(mod$formula$formula,
                        log_size,
                        family = negbinomial(),
                        save_model = "ipms/Stan/flow_n_times_sw2_seas_gam.stan")
  
  # Initialize the spline coefficients and the basis matrices
  draws <- as.data.frame(mod)
  Z_1 <- data$Zs_1_1[ind, ]
  S_1 <- prep_draws("s_t2temp_dry_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  Z_2 <- data$Zs_2_1[ind, ]
  S_2 <- prep_draws("s_t2temp_wet_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  Z_3 <- data$Zs_3_1[ind, ]
  S_3 <- prep_draws("s_t2prec_dry_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  Z_4 <- data$Zs_4_1[ind, ]
  S_4 <- prep_draws("s_t2prec_wet_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  Z_5 <- data$Zs_5_1[ind, ]
  S_5 <- prep_draws("s_t2sw2_dry_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  Z_6 <- data$Zs_6_1[ind, ]
  S_6 <- prep_draws("s_t2sw2_wet_t_1log_size_1\\[", draws, draw_id, sum_fun)
  
  # Fixed Effects
  globe_int <- prep_draws("b_Intercept", draws, draw_id, sum_fun)
  Bs        <- c(prep_draws("b_log_size", draws, draw_id, sum_fun),
                 prep_draws("b_native", draws, draw_id, sum_fun))
  X         <- data$X[ind , 2:3]
  site_int  <- ifelse(is.null(draw_id), 
                      mean(draws[ ,grepl(site, names(draws))]),
                      draws[draw_id, grepl(site, names(draws))])
  
  
  lin_p <- globe_int +
    (Z_1 %*% S_1) + 
    (Z_2 %*% S_2) + 
    (Z_3 %*% S_3) + 
    (Z_4 %*% S_4) + 
    (Z_5 %*% S_5) +
    (Z_6 %*% S_6) +
    (X %*% Bs) +
    site_int
  
  as.vector(exp(lin_p))
  
}

pred_grow_gam <-  function(mod, z, clim_data, site,
                           draw_id = NULL, sum_fun = mean) {
  
  site <- deparse(substitute(site))
  
  log_size <- expand.grid(log_size = z,
                          site = unique(clim_data$site)) %>%
    mutate(
      native = case_when(
        site %in% c("Melkboss",      "Vogelgat", 
                    "St_Francis",    "Struisbaai",
                    "Springfontein", "Rooisand")   ~ 1,
        TRUE ~ 0
      ),
      log_size_next = 5 # Create dummy variable for survival so standata doesn't complain
    )  %>%
    left_join(clim_data, by = "site")
  
  ind <- which(log_size$site == site)
  
  form <- bf(mod$formula$formula,
             mod$formula$pforms)
  
  data <- make_standata(form,
                        log_size,
                        family = gaussian())
  
  
  # Initialize the spline coefficients and the basis matrices
  draws <- as.data.frame(mod)
  
  Z_1 <- data$Zs_1_1[ind, ]
  S_1 <- prep_draws("s_t2seas_temp_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_2 <- data$Zs_2_1[ind, ]
  S_2 <- prep_draws("s_t2total_prec_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_3 <- data$Zs_3_1[ind, ]
  S_3 <- prep_draws("s_t2seas_prec_tlog_size_1\\[", draws, draw_id, sum_fun)
  
  Z_4 <- data$Zs_4_1[ind, ]
  S_4 <- prep_draws("^s_smean_sw2_t_1\\[", draws, draw_id, sum_fun)
  
  Z_5 <- data$Zs_5_1[ind, ]
  S_5 <- prep_draws("^s_smean_sw1_t_1\\[", draws, draw_id, sum_fun)
  
  globe_int <- prep_draws("^b_Intercept$", draws, draw_id, sum_fun)
  Bs        <- c(prep_draws("^b_log_size$", draws, draw_id, sum_fun),
                 prep_draws("^b_mean_temp_t$", draws, draw_id, sum_fun),
                 prep_draws("^b_seas_sw2_t$", draws, draw_id, sum_fun),
                 prep_draws("^b_seas_sw1_t", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:mean_temp_t$", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:seas_sw2_t$", draws, draw_id, sum_fun),
                 prep_draws("^b_log_size:seas_sw1_t$" , draws, draw_id, sum_fun))
  X         <- data$X[ind , -1]
  
  Bs_sigma  <- c(prep_draws("^b_sigma_Intercept$", draws, draw_id, sum_fun),
                 prep_draws("^b_sigma_log_size$", draws, draw_id, sum_fun))
  X_sigma   <- data$X_sigma[ind, ]
  
  # Random effects for site and sigma
  site_int_nm     <- paste0("^r_site\\[", site)
  site_sig_int_nm <- paste0("^r_site__sigma\\[", site)
  
  site_int  <- ifelse(is.null(draw_id), 
                      mean(draws[ ,grepl(site_int_nm, names(draws))]),
                      draws[draw_id, grepl(site_int_nm, names(draws))])
  
  site_sig_int  <- ifelse(is.null(draw_id), 
                          mean(draws[ ,grepl(site_sig_int_nm, names(draws))]),
                          draws[draw_id, grepl(site_sig_int_nm, names(draws))])
  
  lin_p <- globe_int +
    (Z_1 %*% S_1) + 
    (Z_2 %*% S_2) + 
    (Z_3 %*% S_3) + 
    (Z_4 %*% S_4) + 
    (Z_5 %*% S_5) + 
    (X %*% Bs) +
    site_int
  
  sigma_lin_p <- (X_sigma %*% Bs_sigma) + site_sig_int
  
  list(mean_est = lin_p,
       sigma    = exp(sigma_lin_p))
  
}

inv_logit <- function(x) {
  1 / (1 + exp(-(x)))
}
