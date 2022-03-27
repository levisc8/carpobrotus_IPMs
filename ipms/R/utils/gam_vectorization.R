source("ipms/R/02_Utility_Functions_and_Dependencies.R")

surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig_gam.rds")
surv_mod <- surv_mods$clim_krig$times_sw2_ann

repr_mod <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")
grow_mod <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
flow_mod <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig.rds") %>%
  .[[1]] %>%
  .[["times_sw2_ann"]]


all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")

clim_data <- all_ramets %>%
  select(site:seas_sw3_t) %>%
  group_by(site) %>%
  distinct() %>%
  ungroup() 

# Want to get Z_<i>_1 %*% s_log_size_X_<i> so that we can vectorize
# calculations of predicted values and avoid costly additional computations
# The general approach is to extract basis coefficient vector (S_<i>) and 
# generate dummy basis matrices (Z_<i>).

# As best I can tell, predict.brms(gam,...) (or any other kind of model, for
# that matter) generates predictions (i.e. 0,1) for every `newdata` point, then
# averages those to get the predicted value/probability. In general, in IPMs, we
# use a point estimate of each coefficient to get a predicted probability, and
# then use the complete set of posterior draws to generate n_iterations number
# of other predicted probabilities to define a distribution of possible values.
# Thus, we can define predict function that uses the spline/ranefs and fixefs
# from each iteration of the mcmc, or that uses the mean of these draws (point
# estimate) to implement an IPM.

z <- -7:5

inv_logit <- function(x) {
  1 / (1 + exp(-(x)))
}

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

p_g <- pred_grow_gam(grow_mod, z,  clim_data, Havatselet)
