# Create IPMs from GAMs and linear models

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
det_data_list <- c(clim_list, # For surv, flow_n
                   list(clim_data = use_vr_model(clim_data)),  # for repr, grow_mod
                   all_lin_pars,
                   list(p_e = 1e-2,
                        s_sb = 1e-1),
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
    iterations = 150
  )


print(Sys.time() - start)

obs_lams <- lambda(carp_ipm)
obs_lams
is_conv_to_asymptotic(carp_ipm, tolerance = 1e-4) %>% all()

# Resampling G(L/A)MM Posteriors ------

boot_lams <- as_tibble(t(obs_lams)) %>%
  mutate(obs = "yes")

# Using the complete posterior to compute uncertainty in the IPM. 
# The vital rate models aren't conditional on each other, so each
# can be resampled separately. However, I will resample based on 
# the posterior probability of each parameter set, so that the results
# reflect the relative plausibility of each one.

surv_draws <- as.data.frame(surv_mod) 
surv_prob  <- surv_draws$lp__ %>%
  normalize_lp()
surv_draws <- surv_draws %>%
  select(c(1:16, 18:30)) %>%
  rename_draws("surv") %>%
  setNames(gsub("intercept", "int", names(.)))
names(surv_draws) <- gsub(pattern = "^surv_z$", 
                          "surv_log_size", 
                          x = names(surv_draws))

grow_draws <- as.data.frame(grow_mod)
grow_prob <- grow_draws$lp__ %>%
  normalize_lp()

grow_draws <- grow_draws %>%
  select(1:24, 29:54) %>%
  rename_draws("grow") %>%
  setNames(gsub("intercept", "int", names(.)))

names(grow_draws) <- gsub(pattern = "^grow_z$", 
                          "grow_log_size", 
                          x = names(grow_draws))



repr_draws <- as.data.frame(repr_mod) 
repr_prob  <- repr_draws$lp__ %>%
  normalize_lp()

repr_draws <- repr_draws %>%
  select(1:16, 19:31) %>%
  rename_draws("repr") %>%
  setNames(gsub("intercept", "int", names(.)))

names(repr_draws) <- gsub(pattern = "^repr_z$", 
                          "repr_log_size", 
                          x = names(repr_draws))

flow_draws <- as.data.frame(flow_mod)
flow_prob  <- flow_draws$lp__ %>%
  normalize_lp()

flow_draws <- flow_draws %>%
  select(c(1:16, 47:72)) %>%
  rename_draws("flow") %>%
  setNames(gsub("intercept", "int", names(.)))

names(flow_draws) <- gsub(pattern = "^flow_z$", 
                          "flow_log_size", 
                          x = names(flow_draws))

names(flow_draws)[grepl(",log_size", names(flow_draws))] <- gsub(
  "int", "log_size", names(flow_draws)[grepl(",log_size", names(flow_draws))]
)

names(flow_draws)[grepl(",log_size", names(flow_draws))] <- gsub(
  ",log_size", "", names(flow_draws)[grepl(",log_size", names(flow_draws))]
)

temp_lams <- vector("list", 4000L)

bad_its <- integer()
i_bad   <- 1L

surv_its <- grow_its <- repr_its <- flow_its <- integer(4000L)

start <- Sys.time()

for(i in seq_len(4000)) {
  
  if(i %% 100 == 0) {
    message(i, " of 4000 iterations done.")
  }
  surv_its[i] <- surv_i <- sample(1:nrow(surv_draws), 1, prob = surv_prob)
  grow_its[i] <- grow_i <- sample(1:nrow(grow_draws), 1, prob = grow_prob)
  repr_its[i] <- repr_i <- sample(1:nrow(repr_draws), 1, prob = repr_prob)
  flow_its[i] <- flow_i <- sample(1:nrow(flow_draws), 1, prob = flow_prob)
  
  boot_surv_pars <- as.list(surv_draws[surv_i, ])
  boot_flow_pars <- as.list(flow_draws[flow_i, ])
  boot_grow_pars <- as.list(grow_draws[grow_i, ])
  boot_repr_pars <- as.list(repr_draws[repr_i, ])

  det_data_list <- c(boot_surv_pars, 
                     boot_flow_pars,
                     boot_grow_pars,
                     boot_repr_pars,
                     clim_list, 
                     list(clim_data = use_vr_model(clim_data)),
                     list(p_e = 1e-2,
                          s_sb = 1e-1),
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
        iterations = 150
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

print(Sys.time() - start) # Simple IPM: 2.48h, GAM IPM: 2 days, GLM IPM: 4.4h

all_lams <- bind_rows(boot_lams, temp_lams) %>%
  as.data.frame()

surv_draws[surv_its, ] %>%
  write.csv(., file = "ipms/Model_Fits/posterior_lambda_surv_pars.csv",
            row.names = FALSE,
            quote = FALSE)

grow_draws[grow_its, ] %>%
  write.csv(., file = "ipms/Model_Fits/posterior_lambda_grow_pars.csv",
            row.names = FALSE,
            quote = FALSE)

repr_draws[repr_its, ] %>%
  write.csv(., file = "ipms/Model_Fits/posterior_lambda_repr_pars.csv",
            row.names = FALSE,
            quote = FALSE)

flow_draws[surv_its, ] %>%
  write.csv(., file = "ipms/Model_Fits/posterior_lambda_flow_pars.csv",
            row.names = FALSE,
            quote = FALSE)

write.csv(all_lams,
          file = "ipms/Model_Fits/ipms/lin_mod_general_ipm_site_lambdas_prob_resamp.csv",
          row.names = FALSE,
          quote = FALSE)

