# Create IPMs from GAMs and linear models

# Load all the selected models and get them into the right list format

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")
all_sdls   <- readRDS("ipms/Data/seedlings.rds")
all_surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
all_flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig.rds")

surv_mod <- all_surv_mods[[1]][["times_sw2_seas"]]
flow_mod <- all_flow_mods[[1]][["times_sw2_ann"]]

grow_mod <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
repr_mod <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")
recr_mod <- readRDS("ipms/Model_Fits/recr_size_brm.rds")

seed_pars <- readRDS("ipms/Model_Fits/discrete_pars/seed_pars.rds")
p_m_pars  <- readRDS("ipms/Model_Fits/discrete_pars/pod_maturation_pars.rds")

data_list <- list(
  grow_mod       = grow_mod, # These can use pred_* functions from util-ipm.R
  repr_mod       = repr_mod  # These can use pred_* functions from util-ipm.R
)

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
                      paste0("surv_intercept_", rownames(ranef(surv_mod)[[1]])))

surv_pars <- c(surv_list, surv_rans)

base_nms <- rownames(fixef(flow_mod)) %>%
  gsub(pattern = "_t_1$", "", x = .) %>%
  gsub(pattern = "\\:", "_x_", x = .)

flow_nms <- paste0("flow_", base_nms) %>%
  gsub(pattern = "Intercept", "int", x = .)

flow_list <- setNames(as.list(fixef(flow_mod)[ , 1]), flow_nms)
flow_rans <- setNames(as.list(ranef(flow_mod)[[1]][ , 1, 1]),
                      paste0("flow_intercept_", rownames(ranef(flow_mod)[[1]])))

all_lin_pars <- c(surv_pars, flow_list, flow_rans, recr_list)

# For the simple determinstic models created on the first pass, we need
# to create a list of the environmental parameters too. These should have the form
# <env_name_site>.

clim_data <- all_ramets %>%
  select(site:seas_sw3_t) %>%
  group_by(site) %>%
  distinct() %>%
  ungroup() 

write.csv(clim_data, 
          file = "ipms/Model_Fits/ipms/site_climate_values.csv",
          row.names = FALSE,
          quote = FALSE)

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
det_data_list <- c(data_list, 
                   clim_list, # For surv, flow_n
                   list(clim_data = use_vr_model(clim_data)),  # for repr, grow_mod
                   all_lin_pars,
                   list(p_e = 1e-5,
                        s_sb = 1e-2),
                   p_m_pars,
                   seed_pars,
                   recr_list)

# Set up domain info
all_z <- c(all_ramets$log_size, all_ramets$log_size_next, all_sdls$log_size_next)
L_z <- min(all_z, na.rm = TRUE) * 1.2
U_z <- max(all_z, na.rm = TRUE) * 1.2

start <- Sys.time()

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
        surv_sw2_dry  * sw2_dry_t_population  +
        surv_sw2_wet  * sw2_wet_t_population  +
        # Interaction effects
        surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
        surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
        surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
        surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
        surv_log_size_x_sw2_dry  * sw2_dry_t_population  * z_1 +
        surv_log_size_x_sw2_wet  * sw2_wet_t_population  * z_1 +
        # Nativity
        surv_native * is_native(population) +
        # Random site effect
        surv_intercept_population 
    ),
    all_g = pred_grow_gam(grow_mod,
                          z_1,
                          clim_data, 
                          site = population),
    g_sd_population = all_g$sigma,
    g_mu_population = all_g$mean_est,
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
    p_repr_population = pred_repr_gam(repr_mod,
                                      z_1,
                                      clim_data, 
                                      site = population),
    s_population = inv_logit(
      # Fixed Effects
      surv_int + 
        surv_log_size * z_1 +
        surv_temp_dry * temp_dry_t_population +
        surv_temp_wet * temp_wet_t_population +
        surv_prec_dry * prec_dry_t_population +
        surv_prec_wet * prec_wet_t_population +
        surv_sw2_dry  * sw2_dry_t_population  +
        surv_sw2_wet  * sw2_wet_t_population  +
        # Interaction effects
        surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
        surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
        surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
        surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
        surv_log_size_x_sw2_dry  * sw2_dry_t_population  * z_1 +
        surv_log_size_x_sw2_wet  * sw2_wet_t_population  * z_1 +
        # Nativity
        surv_native * is_native(population) +
        # Random site effect
        surv_intercept_population
      ),
    n_flow_population = exp(
      # Fixed effects
      flow_int +
        flow_log_size   * z_1 +
        flow_mean_temp  * mean_temp_t_1_population +
        flow_seas_temp  * seas_temp_t_1_population +
        flow_total_prec * total_prec_t_1_population +
        flow_seas_prec  * seas_prec_t_1_population  +
        flow_mean_sw2   * mean_sw2_t_1_population +
        flow_seas_sw2   * seas_sw2_t_1_population +
        # Interactions 
        flow_log_size_x_mean_temp  * mean_temp_t_1_population  * z_1 +
        flow_log_size_x_seas_temp  * seas_temp_t_1_population  * z_1 +
        flow_log_size_x_total_prec * total_prec_t_1_population * z_1 +
        flow_log_size_x_seas_prec  * seas_prec_t_1_population  * z_1 +
        flow_log_size_x_mean_sw2   * mean_sw2_t_1_population   * z_1 +
        flow_log_size_x_seas_sw2   * seas_sw2_t_1_population   * z_1 +
        # Nativity
        flow_native * is_native(population) +
        # Random site intercept
        flow_intercept_population
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
    iterations = 75
  )

print(Sys.time() - start)

obs_lams <- lambda(carp_ipm)
obs_lams
is_conv_to_asymptotic(carp_ipm, tolerance = 1e-4)

# Resampling G(L/A)MM Posteriors ------

boot_lams <- as_tibble(t(obs_lams)) %>%
  mutate(obs = "yes")

# Using the complete posterior to compute uncertainty in the IPM. 
# The vital rate models aren't conditional on each other, so each
# can be resampled separately. 

surv_draws <- as.data.frame(surv_mod) %>%
  select(c(1:16, 19:31)) %>%
  rename_draws("surv") 
names(surv_draws) <- gsub(pattern = "^surv_z$", 
                          "surv_log_size", 
                          x = names(surv_draws))

flow_draws <- as.data.frame(flow_mod) %>%
  select(c(1:16, 20:32)) %>%
  rename_draws("flow") 
names(flow_draws) <- gsub(pattern = "^flow_z$", 
                          "flow_log_size", 
                          x = names(flow_draws))

temp_lams <- vector("list", 4000L)

bad_its <- integer()
i_bad   <- 1L

start <- Sys.time()

for(i in seq_len(4000)) {
  
  if(i %% 100 == 0) {
    message(i, " of 4000 iterations done.")
  }
  
  boot_surv_pars <- as.list(surv_draws[i, ])
  boot_flow_pars <- as.list(flow_draws[i, ])

  det_data_list <- c(boot_surv_pars, 
                     boot_flow_pars,
                     clim_list, 
                     list(grow_mod = grow_mod,
                          repr_mod = repr_mod,
                          iter = i),
                     list(clim_data = use_vr_model(clim_data)),
                     list(p_e = 1e-5,
                          s_sb = 1e-2),
                     p_m_pars,
                     seed_pars,
                     recr_list)
  
  # Some iterations apparently produce negative/NAs when built? 
  # I'm not totally sure how that's possible, but wrapup in tryCatch() and save
  # to inspect later
  
  boot_ipm <- tryCatch({
    init_ipm("general", "di", "det") %>%
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
            surv_sw2_dry  * sw2_dry_t_population  +
            surv_sw2_wet  * sw2_wet_t_population  +
            # Interaction effects
            surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
            surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
            surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
            surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
            surv_log_size_x_sw2_dry  * sw2_dry_t_population  * z_1 +
            surv_log_size_x_sw2_wet  * sw2_wet_t_population  * z_1 +
            # Nativity
            surv_native * is_native(population) +
            # Random site effect
            surv_intercept_population 
        ),
        all_g = pred_grow_gam(grow_mod,
                              z_1,
                              clim_data, 
                              site = population),
        g_sd_population = all_g$sigma,
        g_mu_population = all_g$mean_est,
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
        p_repr_population = pred_repr_gam(repr_mod,
                                          z_1,
                                          clim_data, 
                                          site = population),
        s_population = inv_logit(
          # Fixed Effects
          surv_int + 
            surv_log_size * z_1 +
            surv_temp_dry * temp_dry_t_population +
            surv_temp_wet * temp_wet_t_population +
            surv_prec_dry * prec_dry_t_population +
            surv_prec_wet * prec_wet_t_population +
            surv_sw2_dry  * sw2_dry_t_population  +
            surv_sw2_wet  * sw2_wet_t_population  +
            # Interaction effects
            surv_log_size_x_temp_dry * temp_dry_t_population * z_1 +
            surv_log_size_x_temp_wet * temp_wet_t_population * z_1 +
            surv_log_size_x_prec_dry * prec_dry_t_population * z_1 +
            surv_log_size_x_prec_wet * prec_wet_t_population * z_1 +
            surv_log_size_x_sw2_dry  * sw2_dry_t_population  * z_1 +
            surv_log_size_x_sw2_wet  * sw2_wet_t_population  * z_1 +
            # Nativity
            surv_native * is_native(population) +
            # Random site effect
            surv_intercept_population
        ),
        n_flow_population = exp(
          # Fixed effects
          flow_int +
            flow_log_size   * z_1 +
            flow_mean_temp  * mean_temp_t_1_population +
            flow_seas_temp  * seas_temp_t_1_population +
            flow_total_prec * total_prec_t_1_population +
            flow_seas_prec  * seas_prec_t_1_population  +
            flow_mean_sw2   * mean_sw2_t_1_population +
            flow_seas_sw2   * seas_sw2_t_1_population +
            # Interactions 
            flow_log_size_x_mean_temp  * mean_temp_t_1_population  * z_1 +
            flow_log_size_x_seas_temp  * seas_temp_t_1_population  * z_1 +
            flow_log_size_x_total_prec * total_prec_t_1_population * z_1 +
            flow_log_size_x_seas_prec  * seas_prec_t_1_population  * z_1 +
            flow_log_size_x_mean_sw2   * mean_sw2_t_1_population   * z_1 +
            flow_log_size_x_seas_sw2   * seas_sw2_t_1_population   * z_1 +
            # Nativity
            flow_native * is_native(population) +
            # Random site intercept
            flow_intercept_population
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
        iterations = 75
      )
  },
  error = function(e) {
    i       <- get("i", envir = .GlobalEnv, inherits = FALSE)
    bad_its <- get("bad_its", envir = .GlobalEnv, inherits = FALSE)
    bad_its <<- c(bad_its, i)
    
    return(e)
  })
  
  if(inherits(boot_ipm, "error")) next
  
  temp_lams[[i]] <- as_tibble(t(lambda(boot_ipm))) %>%
    mutate(obs = "no")
  
  if(!all(is_conv_to_asymptotic(boot_ipm, tolerance = 1e-4))) {
    message("Iteration ", i, " has at least one unconverged model!")
  }
  
}

print(Sys.time() - start) # Simple IPM: 2.48h, general IPM: 

all_lams <- bind_rows(boot_lams, temp_lams) %>%
  as.data.frame()

write.csv(all_lams,
          file = "ipms/Model_Fits/ipms/gam_mod_general_ipm_site_lambdas.csv",
          row.names = FALSE,
          quote = FALSE)

