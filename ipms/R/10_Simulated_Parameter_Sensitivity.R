# Compute sensitivity of deterministic models to parameter value assumptions
# for p_e and s_sb


# Load all the selected models and get them into the right list format

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")
all_sdls   <- readRDS("ipms/Data/seedlings.rds")
all_surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
all_flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_lin_slope_ran.rds")
all_grow_mods <- readRDS("ipms/Model_Fits/ramet_growth_list_krig.rds")
all_repr_mods <- readRDS("ipms/Model_Fits/ramet_repro_list_krig.rds")


surv_mod <- all_surv_mods[[1]][["times_sw3_seas"]]
flow_mod <- all_flow_mods[[1]][["times_sw3_seas"]]
grow_mod <- all_grow_mods[[1]][["times_sw3_seas"]]
repr_mod <- all_repr_mods[[1]][["times_sw1_seas"]]

recr_mod <- readRDS("ipms/Model_Fits/recr_size_brm.rds")

seed_pars <- readRDS("ipms/Model_Fits/discrete_pars/seed_pars.rds")
p_m_pars  <- readRDS("ipms/Model_Fits/discrete_pars/pod_maturation_pars.rds")

recr_list <- list(recr_mu = -8.62,
                  recr_sd = 0.54,
                  recr_s  = 15 / 74) # nrow(recruits) / sum(all_sdls$n_sdl)

base_nms <- rownames(fixef(surv_mod)) %>%
  gsub(pattern = "_t$", replacement = "", x = .) %>%
  gsub(pattern = "\\:", replacement = "_x_", x = .)

surv_nms <- paste0("surv_", base_nms) %>%
  gsub(pattern = "Intercept", "int", x = .)


surv_list <- setNames(as.list(fixef(surv_mod)[ , 1]), surv_nms)
surv_rans <- setNames(as.list(ranef(surv_mod)[[1]][ , 1, 1]),
                      paste0("surv_int_", rownames(ranef(surv_mod)[[1]])))

surv_pars <- c(surv_list, surv_rans)

base_nms <- rownames(fixef(grow_mod)) %>%
  gsub(pattern = "_t$", replacement = "", x = .) %>%
  gsub(pattern = "\\:", replacement = "_x_", x = .)

grow_nms <- paste0("grow_", base_nms) %>%
  gsub(pattern = "Intercept", "int", x = .)

grow_list <- setNames(as.list(fixef(grow_mod)[ , 1]), grow_nms)
grow_rans <- c(setNames(as.list(ranef(grow_mod)[[1]][ , 1, 1]), 
                        paste0("grow_int_", rownames(ranef(grow_mod)[[1]]))),
               setNames(as.list(ranef(grow_mod)[[1]][ , 1, 2]), 
                        paste0("grow_sigma_int_", rownames(ranef(grow_mod)[[1]]))))

grow_pars <- c(grow_list, grow_rans)

base_nms <- rownames(fixef(repr_mod)) %>%
  gsub(pattern = "_t_1$", replacement = "", x = .) %>%
  gsub(pattern = "\\:", replacement = "_x_", x = .)

repr_nms <- paste0("repr_", base_nms) %>%
  gsub(pattern = "Intercept", "int", x = .)

repr_list <- setNames(as.list(fixef(repr_mod)[ , 1]), repr_nms)
repr_rans <- setNames(as.list(ranef(repr_mod)[[1]][ , 1, 1]), 
                      paste0("repr_int_", rownames(ranef(repr_mod)[[1]])))

repr_pars <- c(repr_list, repr_rans)

base_nms <- rownames(fixef(flow_mod)) %>%
  gsub(pattern = "_t_1$", "", x = .) %>%
  gsub(pattern = "\\:", "_x_", x = .)

flow_nms <- paste0("flow_", base_nms) %>%
  gsub(pattern = "Intercept", "int", x = .)

flow_list <- setNames(as.list(fixef(flow_mod)[ , 1]), flow_nms)
flow_rans <- c(setNames(as.list(ranef(flow_mod)[[1]][ , 1, 1]),
                        paste0("flow_int_", 
                               rownames(ranef(flow_mod)[[1]][, , 1]))),
               setNames(as.list(ranef(flow_mod)[[1]][ , 1, 2]),
                        paste0("flow_log_size_", 
                               rownames(ranef(flow_mod)[[1]][, , 2])))
               
)

all_lin_pars <- c(surv_pars, grow_pars, repr_pars,
                  flow_list, flow_rans, recr_list)

# For the simple determinstic models created on the first pass, we need
# to create a list of the environmental parameters too. These should have the form
# <env_name_site>.

clim_data <- all_ramets %>%
  select(site:seas_sw3_t) %>%
  group_by(site) %>%
  distinct() %>%
  ungroup() 

clim_list <- list()

for(i in seq_len(nrow(clim_data))) {
  site <- clim_data$site[i]
  temp <- as.list(clim_data[i, -1]) %>%
    setNames(paste0(names(.), "_", site))
  
  clim_list <- c(clim_list, temp)
  
}

# Next, get sigma growth intercept parameters from model. These aren't
# smooths, so we can just use exp(intercept + slope * z_1 + site_int)


# Data list for deterministic IPMs


# Set up domain info
all_z <- c(all_ramets$log_size, all_ramets$log_size_next, all_sdls$log_size_next)
L_z <- min(all_z, na.rm = TRUE) * 1.2
U_z <- max(all_z, na.rm = TRUE) * 1.2

p_e_vals <- c(seq(1e-3, 1e-1, by = 5e-3), 1e-1)
p_e_seq <- expand.grid(p_e  = p_e_vals,
                       site = unique(all_ramets$site),
                       stringsAsFactors = FALSE)

s_sb_vals <- c(seq(1e-2, 1, by = 5e-2), 1)
s_sb_seq <- expand.grid(s_sb = s_sb_vals,
                        site = unique(all_ramets$site),
                        stringsAsFactors = FALSE)

p_e_out <- cbind(p_e_seq, 
                 lam = NA)
s_sb_out <- cbind(s_sb_seq,
                  lam = NA)

for(i in seq_along(p_e_vals)) {
  
  det_data_list <- c(clim_list,
                     list(clim_data = use_vr_model(clim_data)),  
                     all_lin_pars,
                     list(p_e = p_e_vals[i],
                          s_sb = 1e-1),
                     p_m_pars,
                     seed_pars,
                     recr_list)
  
  carp_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name    = "P_population",
      formula = s_population * G_population * d_z,
      family  = "CC",
      s_population = inv_logit(
        # Fixed Effects
        surv_int + 
          surv_log_size * z_1 +
          surv_temp_dry * temp_dry_t_population +
          surv_temp_wet * temp_wet_t_population +
          surv_prec_dry * prec_dry_t_population +
          surv_prec_wet * prec_wet_t_population +
          surv_sw3_dry  * sw3_dry_t_population  +
          surv_sw3_wet  * sw3_wet_t_population  +
          # Interaction effects
          surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          surv_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
          surv_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
          # Nativity
          surv_native * is_native(population) +
          # Random site effect
          surv_int_population 
      ),
      
      g_sd_population = exp(grow_sigma_int + grow_sigma_int_population + 
                              grow_sigma_log_size * z_1 + 
                              grow_sigma_temp_dry * temp_dry_t_population +
                              grow_sigma_temp_wet * temp_wet_t_population +
                              grow_sigma_prec_dry * prec_dry_t_population + 
                              grow_sigma_prec_wet * prec_wet_t_population +
                              grow_sigma_sw3_dry *  sw3_dry_t_population +
                              grow_sigma_sw3_wet *  sw3_wet_t_population +
                              grow_sigma_native  *  is_native(population)),
      g_mu_population = grow_int + 
        grow_log_size * z_1 + 
        grow_temp_dry * temp_dry_t_population +
        grow_temp_wet * temp_wet_t_population +
        grow_prec_dry * prec_dry_t_population + 
        grow_prec_wet * prec_wet_t_population +
        grow_sw3_dry *  sw3_dry_t_population +
        grow_sw3_wet *  sw3_wet_t_population +
        # Interactions
        grow_log_size_x_temp_dry * temp_dry_t_population * z_1 +
        grow_log_size_x_temp_wet * temp_wet_t_population * z_1 +
        grow_log_size_x_prec_dry * prec_dry_t_population * z_1 +
        grow_log_size_x_prec_wet * prec_wet_t_population * z_1 +
        grow_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
        grow_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
        
        grow_native  *  is_native(population) +
        grow_int_population,
      G_population    = dnorm(z_2, g_mu_population, g_sd_population),
      # Define the rest of the sub-kernel
      states          = list("z"),
      data_list       = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor       = TRUE,
      evict_fun       = truncated_distributions("norm", "G_population")
    ) %>%
    define_kernel(
      name = "P_to_M_population",
      formula = p_repr_population * s_population * 
        n_flow_population * p_m_population * d_z,
      family = "CD",
      p_repr_population = inv_logit(
        repr_int + 
          repr_log_size * z_1 +
          repr_temp_dry * temp_dry_t_population +
          repr_temp_wet * temp_wet_t_population +
          repr_prec_dry * prec_dry_t_population +
          repr_prec_wet * prec_wet_t_population +
          repr_sw1_dry  * sw1_dry_t_population  +
          repr_sw1_wet  * sw1_wet_t_population  +
          # Interaction effects
          repr_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          repr_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          repr_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          repr_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          repr_log_size_x_sw1_dry  * sw1_dry_t_population  * z_1 +
          repr_log_size_x_sw1_wet  * sw1_wet_t_population  * z_1 +
          # Nativity
          repr_native * is_native(population) +
          repr_log_size_x_native * is_native(population) * z_1 +
          # Random site effect
          repr_int_population
      ),
      s_population = inv_logit(
        # Fixed Effects
        surv_int + 
          surv_log_size * z_1 +
          surv_temp_dry * temp_dry_t_population +
          surv_temp_wet * temp_wet_t_population +
          surv_prec_dry * prec_dry_t_population +
          surv_prec_wet * prec_wet_t_population +
          surv_sw3_dry  * sw3_dry_t_population  +
          surv_sw3_wet  * sw3_wet_t_population  +
          # Interaction effects
          surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          surv_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
          surv_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
          # Nativity
          surv_native * is_native(population) +
          # Random site effect
          surv_int_population
      ),
      n_flow_population = exp(
        # Fixed effects
        flow_int +
          (flow_log_size + flow_log_size_population) * z_1 +
          flow_temp_wet * temp_wet_t_1_population +
          flow_temp_dry * temp_dry_t_1_population +
          flow_prec_wet * prec_wet_t_1_population +
          flow_prec_dry * prec_dry_t_1_population  +
          flow_sw3_wet  * sw3_wet_t_1_population +
          flow_sw3_dry  * sw3_dry_t_1_population +
          # Interactions 
          flow_log_size_x_temp_wet  * temp_wet_t_1_population * z_1 +
          flow_log_size_x_temp_dry  * temp_dry_t_1_population * z_1 +
          flow_log_size_x_prec_wet * prec_wet_t_1_population * z_1 +
          flow_log_size_x_prec_dry  * prec_dry_t_1_population * z_1 +
          flow_log_size_x_sw3_wet   * sw3_wet_t_1_population  * z_1 +
          flow_log_size_x_sw3_dry   * sw3_dry_t_1_population  * z_1 +
          # Nativity
          flow_native * is_native(population) +
          # Random site intercept
          flow_int_population
      ),
      recr_size       = dnorm(z_2, recr_mu, recr_sd),
      states          = list(c("z", "mf")),
      data_list       = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name      = "M_to_SDL",
      formula   = r_f * p_e * g_i * v_s,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("mf", "sdl"))
    ) %>%
    define_kernel(
      name      = "M_to_SB",
      formula   = r_f * (1 - g_i) * v_s,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("mf", "sb"))
    ) %>%
    define_kernel(
      name      = "SB_to_SB",
      formula   = s_sb * (1 - g_sb),
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("sb"))
    ) %>%
    define_kernel(
      name      = "SB_to_SDL",
      formula   = s_sb * g_sb * p_e,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("sb", "sdl"))
    ) %>%
    define_kernel(
      name      = "SDL_to_P",
      formula   = recr_s * r_d,
      family    = "DC",
      r_d       = dnorm(z_2, recr_mu, recr_sd),
      data_list = det_data_list,
      states    = list(c("sdl", "z")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm",
                                          "r_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        c("P_population", "P_to_M_population",
          "M_to_SDL", "M_to_SB",
          "SB_to_SB", "SB_to_SDL",
          "SDL_to_P"),
        rep("midpoint", 7),
        state_start = c("z", "z",  "mf",  "mf", "sb", "sb", "sdl"),
        state_end   =  c("z", "mf", "sdl", "sb", "sb", "sdl", "z")
      )
    ) %>%
    define_domains(z = c(L_z, U_z, 100)) %>%
    define_pop_state(n_z_population   = runif(100),
                     n_mf_population  = rpois(1, 20),
                     n_sdl_population = rpois(1, 10),
                     n_sb_population  = rpois(1, 10)) %>%
    make_ipm(
      usr_funs = list(
        inv_logit     = inv_logit,
        pred_grow_gam = pred_grow_gam,
        pred_repr_gam = pred_repr_gam,
        is_native     = is_native
      ),
      iterations = 200
    )
  
  all_lams <- lambda(carp_ipm) %>%
    setNames(gsub("lambda_", "", names(.)))
  
  conv <- is_conv_to_asymptotic(carp_ipm)
  
  if(!all(conv)) {
    
    nms <- gsub("lambda_", "", names(conv[!conv])) %>% paste(collapse = ", ")
    message("p_e = ", p_e_vals[i], ": ", nms, " did not converge!\n")
  }
  
  for(j in seq_along(all_lams)) {
    
    p_e_out[p_e_out$site == names(all_lams)[j] & 
              p_e_out$p_e == p_e_vals[i], "lam"] <- round(all_lams[j], 5)
    
  }
  
  
}



for(i in seq_along(s_sb_vals)) {
  
  det_data_list <- c(clim_list,
                     list(clim_data = use_vr_model(clim_data)),  
                     all_lin_pars,
                     list(p_e = 1e-2,
                          s_sb = s_sb_vals[i]),
                     p_m_pars,
                     seed_pars,
                     recr_list)
  
  carp_ipm <- init_ipm("general", "di", "det") %>%
    define_kernel(
      name    = "P_population",
      formula = s_population * G_population * d_z,
      family  = "CC",
      s_population = inv_logit(
        # Fixed Effects
        surv_int + 
          surv_log_size * z_1 +
          surv_temp_dry * temp_dry_t_population +
          surv_temp_wet * temp_wet_t_population +
          surv_prec_dry * prec_dry_t_population +
          surv_prec_wet * prec_wet_t_population +
          surv_sw3_dry  * sw3_dry_t_population  +
          surv_sw3_wet  * sw3_wet_t_population  +
          # Interaction effects
          surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          surv_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
          surv_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
          # Nativity
          surv_native * is_native(population) +
          # Random site effect
          surv_int_population 
      ),
      
      g_sd_population = exp(grow_sigma_int + grow_sigma_int_population + 
                              grow_sigma_log_size * z_1 + 
                              grow_sigma_temp_dry * temp_dry_t_population +
                              grow_sigma_temp_wet * temp_wet_t_population +
                              grow_sigma_prec_dry * prec_dry_t_population + 
                              grow_sigma_prec_wet * prec_wet_t_population +
                              grow_sigma_sw3_dry *  sw3_dry_t_population +
                              grow_sigma_sw3_wet *  sw3_wet_t_population +
                              grow_sigma_native  *  is_native(population)),
      g_mu_population = grow_int + 
        grow_log_size * z_1 + 
        grow_temp_dry * temp_dry_t_population +
        grow_temp_wet * temp_wet_t_population +
        grow_prec_dry * prec_dry_t_population + 
        grow_prec_wet * prec_wet_t_population +
        grow_sw3_dry *  sw3_dry_t_population +
        grow_sw3_wet *  sw3_wet_t_population +
        # Interactions
        grow_log_size_x_temp_dry * temp_dry_t_population * z_1 +
        grow_log_size_x_temp_wet * temp_wet_t_population * z_1 +
        grow_log_size_x_prec_dry * prec_dry_t_population * z_1 +
        grow_log_size_x_prec_wet * prec_wet_t_population * z_1 +
        grow_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
        grow_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
        
        grow_native  *  is_native(population) +
        grow_int_population,
      G_population    = dnorm(z_2, g_mu_population, g_sd_population),
      # Define the rest of the sub-kernel
      states          = list("z"),
      data_list       = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor       = TRUE,
      evict_fun       = truncated_distributions("norm", "G_population")
    ) %>%
    define_kernel(
      name = "P_to_M_population",
      formula = p_repr_population * s_population * 
        n_flow_population * p_m_population * d_z,
      family = "CD",
      p_repr_population = inv_logit(
        repr_int + 
          repr_log_size * z_1 +
          repr_temp_dry * temp_dry_t_population +
          repr_temp_wet * temp_wet_t_population +
          repr_prec_dry * prec_dry_t_population +
          repr_prec_wet * prec_wet_t_population +
          repr_sw1_dry  * sw1_dry_t_population  +
          repr_sw1_wet  * sw1_wet_t_population  +
          # Interaction effects
          repr_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          repr_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          repr_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          repr_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          repr_log_size_x_sw1_dry  * sw1_dry_t_population  * z_1 +
          repr_log_size_x_sw1_wet  * sw1_wet_t_population  * z_1 +
          # Nativity
          repr_native * is_native(population) +
          repr_log_size_x_native * is_native(population) * z_1 +
          # Random site effect
          repr_int_population
      ),
      s_population = inv_logit(
        # Fixed Effects
        surv_int + 
          surv_log_size * z_1 +
          surv_temp_dry * temp_dry_t_population +
          surv_temp_wet * temp_wet_t_population +
          surv_prec_dry * prec_dry_t_population +
          surv_prec_wet * prec_wet_t_population +
          surv_sw3_dry  * sw3_dry_t_population  +
          surv_sw3_wet  * sw3_wet_t_population  +
          # Interaction effects
          surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
          surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
          surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
          surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
          surv_log_size_x_sw3_dry  * sw3_dry_t_population  * z_1 +
          surv_log_size_x_sw3_wet  * sw3_wet_t_population  * z_1 +
          # Nativity
          surv_native * is_native(population) +
          # Random site effect
          surv_int_population
      ),
      n_flow_population = exp(
        # Fixed effects
        flow_int +
          (flow_log_size + flow_log_size_population) * z_1 +
          flow_temp_wet * temp_wet_t_1_population +
          flow_temp_dry * temp_dry_t_1_population +
          flow_prec_wet * prec_wet_t_1_population +
          flow_prec_dry * prec_dry_t_1_population  +
          flow_sw3_wet  * sw3_wet_t_1_population +
          flow_sw3_dry  * sw3_dry_t_1_population +
          # Interactions 
          flow_log_size_x_temp_wet  * temp_wet_t_1_population * z_1 +
          flow_log_size_x_temp_dry  * temp_dry_t_1_population * z_1 +
          flow_log_size_x_prec_wet * prec_wet_t_1_population * z_1 +
          flow_log_size_x_prec_dry  * prec_dry_t_1_population * z_1 +
          flow_log_size_x_sw3_wet   * sw3_wet_t_1_population  * z_1 +
          flow_log_size_x_sw3_dry   * sw3_dry_t_1_population  * z_1 +
          # Nativity
          flow_native * is_native(population) +
          # Random site intercept
          flow_int_population
      ),
      recr_size       = dnorm(z_2, recr_mu, recr_sd),
      states          = list(c("z", "mf")),
      data_list       = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name      = "M_to_SDL",
      formula   = r_f * p_e * g_i * v_s,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("mf", "sdl"))
    ) %>%
    define_kernel(
      name      = "M_to_SB",
      formula   = r_f * (1 - g_i) * v_s,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("mf", "sb"))
    ) %>%
    define_kernel(
      name      = "SB_to_SB",
      formula   = s_sb * (1 - g_sb),
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("sb"))
    ) %>%
    define_kernel(
      name      = "SB_to_SDL",
      formula   = s_sb * g_sb * p_e,
      family    = "DD",
      data_list = det_data_list,
      states    = list(c("sb", "sdl"))
    ) %>%
    define_kernel(
      name      = "SDL_to_P",
      formula   = recr_s * r_d,
      family    = "DC",
      r_d       = dnorm(z_2, recr_mu, recr_sd),
      data_list = det_data_list,
      states    = list(c("sdl", "z")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm",
                                          "r_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        c("P_population", "P_to_M_population",
          "M_to_SDL", "M_to_SB",
          "SB_to_SB", "SB_to_SDL",
          "SDL_to_P"),
        rep("midpoint", 7),
        state_start = c("z", "z",  "mf",  "mf", "sb", "sb", "sdl"),
        state_end   =  c("z", "mf", "sdl", "sb", "sb", "sdl", "z")
      )
    ) %>%
    define_domains(z = c(L_z, U_z, 100)) %>%
    define_pop_state(n_z_population   = runif(100),
                     n_mf_population  = rpois(1, 20),
                     n_sdl_population = rpois(1, 10),
                     n_sb_population  = rpois(1, 10)) %>%
    make_ipm(
      usr_funs = list(
        inv_logit     = inv_logit,
        pred_grow_gam = pred_grow_gam,
        pred_repr_gam = pred_repr_gam,
        is_native     = is_native
      ),
      iterations = 200
    )
  
  all_lams <- lambda(carp_ipm) %>%
    setNames(gsub("lambda_", "", names(.)))
  
  conv <- is_conv_to_asymptotic(carp_ipm)
  
  if(!all(conv)) {
    
    nms <- gsub("lambda_", "", names(conv[!conv])) %>% paste(collapse = ", ")
    message("s_sb = ", s_sb_vals[i], ": ", nms, " did not converge!\n")
  }
  
  for(j in seq_along(all_lams)) {
    
    s_sb_out[s_sb_out$site == names(all_lams)[j] & 
              s_sb_out$s_sb == s_sb_vals[i], "lam"] <- round(all_lams[j], 5)
    
  }
  
  
}

p_e_plt <- ggplot(p_e_out, 
                  aes(x = p_e, y = lam)) + 
  geom_line() +
  facet_wrap(~site,
             scales = "free_y") +
  theme_bw() + 
  geom_vline(xintercept = 1e-2, color = "red", linetype = "dotted", size = 1.3) +
  scale_x_log10("Establishment probability")  +
  ylab(expression("Per-capita Growth Rate (" ~ lambda ~ ")"))

s_sb_plt <- ggplot(s_sb_out, 
                  aes(x = s_sb, y = lam)) + 
  geom_line() +
  facet_wrap(~site,
             scales = "free_y") +
  geom_vline(xintercept = 1e-1, color = "red", linetype = "dotted", size = 1.3) +
  theme_bw() +
  scale_x_log10("Seedbank Survival rate") +
  ylab(expression("Per-capita Growth Rate (" ~ lambda ~ ")"))

pdf("ipms/Figures/ipm_output/p_e_sensitivity_plot_prob_resamp.pdf")

  print(p_e_plt)

dev.off()


pdf("ipms/Figures/ipm_output/s_sb_sensitivity_plot_prob_resamp.pdf")

  print(s_sb_plt)

dev.off()


  