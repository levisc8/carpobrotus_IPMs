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

data_list <- list(
  grow_mod       = grow_mod, # These can use pred_* functions from util-ipm.R
  repr_mod       = repr_mod  # These can use pred_* functions from util-ipm.R
)

recr_list <- list(recr_mu = -8.62,
                  recr_sd = 0.54,
                  p_germ  = 0.01)

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
                   list(clim_data = clim_data),  # for repr, grow_mod
                   all_lin_pars)

# Set up domain info
all_z <- c(all_ramets$log_size, all_ramets$log_size_next, all_sdls$log_size_next)
L_z <- min(all_z, na.rm = TRUE) * 1.2
U_z <- max(all_z, na.rm = TRUE) * 1.2

start <- Sys.time()

carp_ipm <- init_ipm("simple", "di", "det") %>%
  define_kernel(
    name    = "P_population",
    formula = s_population * G_population,
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
    data_list = det_data_list,
    uses_par_sets   = TRUE,
    par_set_indices = list(population = unique(all_ramets$site)),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "G_population")
  ) %>%
  define_kernel(
    name = "F_population",
    formula = p_repr_population * n_flow_population * recr_size * p_germ,
    family = "CC",
    p_repr_population = pred_repr_gam(repr_mod,
                                      z_1,
                                      clim_data, 
                                      site = population),
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
    states          = list("z"),
    data_list = det_data_list,
    uses_par_sets   = TRUE,
    par_set_indices = list(population = unique(all_ramets$site)),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "recr_size")
  ) %>%
  define_impl(
    make_impl_args_list(c("P_population", "F_population"),
                        rep("midpoint", 2),
                        rep("z", 2),
                        rep("z", 2))
  ) %>%
  define_domains(z = c(L_z, U_z, 100)) %>%
  define_pop_state(n_z_population = runif(100)) %>%
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
is_conv_to_asymptotic(carp_ipm, tolerance = 1e-4)

# Resampling GAM Posteriors ------

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

start <- Sys.time()

bad_its <- integer()
i_bad   <- 1L

for(i in seq_len(4000)) {
  
  if(i %% 100 == 0) {
    message(i, " of 4000 iterations done.")
  }
  
  boot_surv_pars <- as.list(surv_draws[i, ])
  boot_flow_pars <- as.list(flow_draws[i, ])

  det_data_list <- c(boot_surv_pars, 
                     boot_flow_pars, 
                     list(grow_mod = grow_mod,
                          repr_mod = repr_mod,
                          iter = i,
                          clim_data = clim_data),
                     recr_list,
                     clim_list)
  
  # Some iterations apparently produce negative/NAs when built? 
  # I'm not totally sure how that's possible, but wrapup in tryCatch() and save
  # to inspect later
  
  boot_ipm <- tryCatch({
    init_ipm("simple", "di", "det") %>%
    define_kernel(
      name    = "P_population",
      formula = s_population * G_population,
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
                            site = population,
                            draw_id = iter),
      g_sd_population = all_g$sigma,
      g_mu_population = all_g$mean_est,
      G_population    = dnorm(z_2, g_mu_population, g_sd_population),
      # Define the rest of the sub-kernel
      states          = list("z"),
      data_list = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G_population")
    ) %>%
    define_kernel(
      name = "F_population",
      formula = p_repr_population * n_flow_population * recr_size * p_germ,
      family = "CC",
      p_repr_population = pred_repr_gam(repr_mod,
                                        z_1,
                                        clim_data, 
                                        site = population,
                                        draw_id = iter),
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
      states          = list("z"),
      data_list = det_data_list,
      uses_par_sets   = TRUE,
      par_set_indices = list(population = unique(all_ramets$site)),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "recr_size")
    ) %>%
    define_impl(
      make_impl_args_list(c("P_population", "F_population"),
                          rep("midpoint", 2),
                          rep("z", 2),
                          rep("z", 2))
    ) %>%
    define_domains(z = c(L_z, U_z, 100)) %>%
    define_pop_state(n_z_population = runif(100)) %>%
    make_ipm(
      usr_funs = list(
        inv_logit     = inv_logit,
        pred_grow_gam = pred_grow_gam,
        pred_repr_gam = pred_repr_gam,
        is_native     = is_native
      ),
      iterations = 100
    )
  },
  error = function(e) {
    i       <- get("i", envir = caller_env())
    bad_its <- get("bad_its", envir = caller_env())
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

print(Sys.time() - start) # ROund 1 was 2.48 hours

all_lams <- bind_rows(boot_lams, temp_lams) %>%
  as.data.frame()

write.csv(all_lams,
          file = "ipms/Model_Fits/ipms/gam_mod_site_lambdas.csv",
          row.names = FALSE,
          quote = FALSE)





# IPM w/ only linear models --------- 

# Trying this for comparison w/ outputs from GAM/linear-based IPM, as the
# lambdas for that are unusually low. This will help work out if it's a weird
# gam thing, or if the continuous stage suffers from some unknown flaw.
# 
# all_flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig.rds")
# all_surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
# all_repr_mods <- readRDS("ipms/Model_Fits/ramet_repro_list_krig.rds")
# all_grow_mods <- readRDS("ipms/Model_Fits/ramet_growth_list_krig.rds")
# 
# flow_mod <- all_flow_mods[[1]][["times_sw3_seas"]]
# surv_mod <- all_surv_mods[[1]][["times_sw3_seas"]]
# repr_mod <- all_repr_mods[[1]][["times_sw1_seas"]]
# grow_mod <- all_grow_mods[[1]][["times_sw3_seas"]]
# 
# 
# grow_pars <- c(name_fixed_pars(grow_mod, "grow", 3),
#                name_ran_pars(grow_mod, "grow"))
# 
# surv_pars <- c(name_fixed_pars(surv_mod, "surv", 3),
#                name_ran_pars(surv_mod, "surv"))
# flow_pars <- c(name_fixed_pars(flow_mod, "flow", 3),
#                name_ran_pars(flow_mod, "flow"))
# 
# repr_pars <- c(name_fixed_pars(repr_mod, "repr", 1),
#                name_ran_pars(repr_mod, "repr"))
# 
# recr_list <- list(
#   
#   recr_mu        = -8.62,
#   recr_sd        = 0.54,
#   p_germ         = 0.01
# )
# 
# inv_logit <- function(x) {
#   1 / (1 + exp(-(x)))
# }
# 
# det_data_list <- c(grow_pars, surv_pars, 
#                    repr_pars, flow_pars, 
#                    recr_list, clim_list)
# 
# start <- Sys.time()
# 
# carp_ipm <- init_ipm("simple", "di", "det") %>%
#   define_kernel(
#     name            = "P_population",
#     formula         = s_population * G_population,
#     family          = "CC",
#     s_population    = inv_logit(
#       # Linear terms
#       surv_int + 
#         surv_z * z_1 + 
#         surv_temp_dry * temp_dry_t_population + 
#         surv_temp_wet * temp_wet_t_population +
#         surv_prec_dry * prec_dry_t_population +
#         surv_prec_wet * prec_wet_t_population +
#         surv_sw3_dry  * sw3_dry_t_population +
#         surv_sw3_wet  * sw3_wet_t_population +
#         # Interactions
#         surv_log_size_x_temp_dry * z_1 * temp_dry_t_population + 
#         surv_log_size_x_temp_wet * z_1 * temp_wet_t_population +
#         surv_log_size_x_prec_dry * z_1 * prec_dry_t_population +
#         surv_log_size_x_prec_wet * z_1 * prec_wet_t_population +
#         surv_log_size_x_sw3_dry  * z_1 * sw3_dry_t_population +
#         surv_log_size_x_sw3_wet  * z_1 * sw3_wet_t_population +
#         # Random site intercept and nativity coefficients
#         surv_intercept_population + 
#         surv_native * is_native(population) +
#         surv_log_size_x_native * z_1 * is_native(population)
#     ),
#     g_mu_population = grow_int + 
#       grow_z * z_1 + 
#       grow_temp_dry * temp_dry_t_population + 
#       grow_temp_wet * temp_wet_t_population +
#       grow_prec_dry * prec_dry_t_population +
#       grow_prec_wet * prec_wet_t_population +
#       grow_sw3_dry  * sw3_dry_t_population +
#       grow_sw3_wet  * sw3_wet_t_population +
#       # Interactions
#       grow_log_size_x_temp_dry * z_1 * temp_dry_t_population + 
#       grow_log_size_x_temp_wet * z_1 * temp_wet_t_population +
#       grow_log_size_x_prec_dry * z_1 * prec_dry_t_population +
#       grow_log_size_x_prec_wet * z_1 * prec_wet_t_population +
#       grow_log_size_x_sw3_dry  * z_1 * sw3_dry_t_population +
#       grow_log_size_x_sw3_wet  * z_1 * sw3_wet_t_population +
#       # Random site intercept and nativity coefficients
#       grow_intercept_population + 
#       grow_native * is_native(population) +
#       grow_log_size_x_native * z_1 * is_native(population),
#     g_sd_population = exp(
#       grow_sigma_int + 
#         grow_sigma_log_size * z_1 +
#         grow_sigma_temp_dry * temp_dry_t_population + 
#         grow_sigma_temp_wet * temp_wet_t_population +
#         grow_sigma_prec_dry * prec_dry_t_population +
#         grow_sigma_prec_wet * prec_wet_t_population +
#         grow_sigma_sw3_dry  * sw3_dry_t_population +
#         grow_sigma_sw3_wet  * sw3_wet_t_population +
#         grow_sigma_native   * is_native(population) +
#         grow_sigma_intercept_population
#         ),
#     G_population    = dnorm(z_2, g_mu_population, g_sd_population),
#     # Define the rest of the sub-kernel
#     states          = list("z"),
#     data_list       = det_data_list,
#     uses_par_sets   = TRUE,
#     par_set_indices = list(population = unique(all_ramets$site)),
#     evict_cor       = TRUE,
#     evict_fun       = truncated_distributions("norm", "G_population")
#   ) %>%
#   define_kernel(
#     name              = "F_population",
#     formula           = p_repr_population * n_flow_population * recr_size * p_germ,
#     family            = "CC",
#     p_repr_population = inv_logit(
#       # Linear terms
#       repr_int + 
#         repr_z * z_1 + 
#         repr_temp_dry * temp_dry_t_1_population + 
#         repr_temp_wet * temp_wet_t_1_population +
#         repr_prec_dry * prec_dry_t_1_population +
#         repr_prec_wet * prec_wet_t_1_population +
#         repr_sw1_dry  * sw3_dry_t_1_population +
#         repr_sw1_wet  * sw3_wet_t_1_population +
#         # Interactions
#         repr_log_size_x_temp_dry * z_1 * temp_dry_t_1_population + 
#         repr_log_size_x_temp_wet * z_1 * temp_wet_t_1_population +
#         repr_log_size_x_prec_dry * z_1 * prec_dry_t_1_population +
#         repr_log_size_x_prec_wet * z_1 * prec_wet_t_1_population +
#         repr_log_size_x_sw1_dry  * z_1 * sw1_dry_t_1_population +
#         repr_log_size_x_sw1_wet  * z_1 * sw1_wet_t_1_population +
#         # Random site intercept and nativity coefficients
#         repr_intercept_population + 
#         repr_native * is_native(population) +
#         repr_log_size_x_native * z_1 * is_native(population)
#     ),
#     n_flow_population = exp(
#       flow_int + 
#         flow_z * z_1 + 
#         flow_temp_dry * temp_dry_t_1_population + 
#         flow_temp_wet * temp_wet_t_1_population +
#         flow_prec_dry * prec_dry_t_1_population +
#         flow_prec_wet * prec_wet_t_1_population +
#         flow_sw3_dry  * sw3_dry_t_1_population +
#         flow_sw3_wet  * sw3_wet_t_1_population +
#         # Interactions
#         flow_log_size_x_temp_dry * z_1 * temp_dry_t_1_population + 
#         flow_log_size_x_temp_wet * z_1 * temp_wet_t_1_population +
#         flow_log_size_x_prec_dry * z_1 * prec_dry_t_1_population +
#         flow_log_size_x_prec_wet * z_1 * prec_wet_t_1_population +
#         flow_log_size_x_sw3_dry  * z_1 * sw3_dry_t_1_population +
#         flow_log_size_x_sw3_wet  * z_1 * sw3_wet_t_1_population +
#         # Random site intercept and nativity coefficients
#         flow_intercept_population + 
#         flow_native * is_native(population) +
#         flow_log_size_x_native * z_1 * is_native(population)
#     ),
#     recr_size         = dnorm(z_2, recr_mu, recr_sd),
#     states            = list("z"),
#     data_list         = det_data_list,
#     uses_par_sets     = TRUE,
#     par_set_indices   = list(population = unique(all_ramets$site)),
#     evict_cor         = TRUE,
#     evict_fun         = truncated_distributions("norm", "recr_size")    
#   ) %>%
#   define_impl(
#     make_impl_args_list(c("P_population", "F_population"),
#                         rep("midpoint", 2),
#                         rep("z", 2),
#                         rep("z", 2))
#   ) %>%
#   define_domains(z = c(L_z, U_z, 150)) %>%
#   define_pop_state(n_z_population = runif(150)) %>%
#   make_ipm(
#     usr_funs = list(
#      inv_logit = inv_logit,
#      is_native = is_native
#     ),
#     iterations = 150,
#     return_all_envs = TRUE
#   )
# 
# Sys.time() - start
# 
# obs_lams <- lambda(carp_ipm)

# Set up resampling uncertainty computations
# 
# boot_lams <- as_tibble(t(obs_lams)) %>%
#   mutate(obs = "yes")
# 
# # Using the complete posterior to compute uncertainty in the IPM. 
# # The vital rate models aren't conditional on each other, so each
# # can be resampled separately. 
# 
# grow_draws <- as.data.frame(grow_mod) %>%
#   select(c(1:25, 30:55)) %>%
#   rename_draws("grow")
# 
# 
# repr_draws <- as.data.frame(repr_mod) %>%
#   select(c(1:16, 19:31)) %>%
#   rename_draws("repr")

# surv_draws <- as.data.frame(surv_mod) %>%
#   select(c(1:16, 19:31)) %>%
#   rename_draws("surv") 
# flow_draws <- as.data.frame(flow_mod) %>%
#   select(c(1:16, 20:32)) %>%
#   rename_draws("flow")
# 
# temp_lams <- vector("list", 4000L)
# 
# start <- Sys.time()
# 
# for(i in seq_len(4000)) {
#   
#   if(i %% 100 == 0) {
#     message(i, " of 4000 iterations done.")
#   }
#   
#   boot_grow_pars <- as.list(grow_draws[i, ])
#   boot_surv_pars <- as.list(surv_draws[i, ])
#   boot_repr_pars <- as.list(repr_draws[i, ])
#   boot_flow_pars <- as.list(flow_draws[i, ])
#   
#   det_data_list <- c(boot_grow_pars, boot_surv_pars, 
#                      boot_repr_pars, boot_flow_pars, 
#                      recr_list, clim_list)
#   
# 
#   boot_ipm <- init_ipm("simple", "di", "det") %>%
#     define_kernel(
#       name            = "P_population",
#       formula         = s_population * G_population,
#       family          = "CC",
#       s_population    = inv_logit(
#         # Linear terms
#         surv_int + 
#           surv_z * z_1 + 
#           surv_temp_dry * temp_dry_t_population + 
#           surv_temp_wet * temp_wet_t_population +
#           surv_prec_dry * prec_dry_t_population +
#           surv_prec_wet * prec_wet_t_population +
#           surv_sw3_dry  * sw3_dry_t_population +
#           surv_sw3_wet  * sw3_wet_t_population +
#           # Interactions
#           surv_log_size_x_temp_dry * z_1 * temp_dry_t_population + 
#           surv_log_size_x_temp_wet * z_1 * temp_wet_t_population +
#           surv_log_size_x_prec_dry * z_1 * prec_dry_t_population +
#           surv_log_size_x_prec_wet * z_1 * prec_wet_t_population +
#           surv_log_size_x_sw3_dry  * z_1 * sw3_dry_t_population +
#           surv_log_size_x_sw3_wet  * z_1 * sw3_wet_t_population +
#           # Random site intercept and nativity coefficients
#           surv_intercept_population + 
#           surv_native * is_native(population) +
#           surv_log_size_x_native * z_1 * is_native(population)
#       ),
#       g_mu_population = grow_int + 
#         grow_z * z_1 + 
#         grow_temp_dry * temp_dry_t_population + 
#         grow_temp_wet * temp_wet_t_population +
#         grow_prec_dry * prec_dry_t_population +
#         grow_prec_wet * prec_wet_t_population +
#         grow_sw3_dry  * sw3_dry_t_population +
#         grow_sw3_wet  * sw3_wet_t_population +
#         # Interactions
#         grow_log_size_x_temp_dry * z_1 * temp_dry_t_population + 
#         grow_log_size_x_temp_wet * z_1 * temp_wet_t_population +
#         grow_log_size_x_prec_dry * z_1 * prec_dry_t_population +
#         grow_log_size_x_prec_wet * z_1 * prec_wet_t_population +
#         grow_log_size_x_sw3_dry  * z_1 * sw3_dry_t_population +
#         grow_log_size_x_sw3_wet  * z_1 * sw3_wet_t_population +
#         # Random site intercept and nativity coefficients
#         grow_intercept_population + 
#         grow_native * is_native(population) +
#         grow_log_size_x_native * z_1 * is_native(population),
#       g_sd_population = exp(
#         grow_sigma_int + 
#           grow_sigma_log_size * z_1 +
#           grow_sigma_temp_dry * temp_dry_t_population + 
#           grow_sigma_temp_wet * temp_wet_t_population +
#           grow_sigma_prec_dry * prec_dry_t_population +
#           grow_sigma_prec_wet * prec_wet_t_population +
#           grow_sigma_sw3_dry  * sw3_dry_t_population +
#           grow_sigma_sw3_wet  * sw3_wet_t_population +
#           grow_sigma_native   * is_native(population) +
#           grow_sigma_intercept_population
#       ),
#       G_population    = dnorm(z_2, g_mu_population, g_sd_population),
#       # Define the rest of the sub-kernel
#       states          = list("z"),
#       data_list       = det_data_list,
#       uses_par_sets   = TRUE,
#       par_set_indices = list(population = unique(all_ramets$site)),
#       evict_cor       = TRUE,
#       evict_fun       = truncated_distributions("norm", "G_population")
#     ) %>%
#     define_kernel(
#       name              = "F_population",
#       formula           = p_repr_population * n_flow_population * recr_size * p_germ,
#       family            = "CC",
#       p_repr_population = inv_logit(
#         # Linear terms
#         repr_int + 
#           repr_z * z_1 + 
#           repr_temp_dry * temp_dry_t_1_population + 
#           repr_temp_wet * temp_wet_t_1_population +
#           repr_prec_dry * prec_dry_t_1_population +
#           repr_prec_wet * prec_wet_t_1_population +
#           repr_sw1_dry  * sw3_dry_t_1_population +
#           repr_sw1_wet  * sw3_wet_t_1_population +
#           # Interactions
#           repr_log_size_x_temp_dry * z_1 * temp_dry_t_1_population + 
#           repr_log_size_x_temp_wet * z_1 * temp_wet_t_1_population +
#           repr_log_size_x_prec_dry * z_1 * prec_dry_t_1_population +
#           repr_log_size_x_prec_wet * z_1 * prec_wet_t_1_population +
#           repr_log_size_x_sw1_dry  * z_1 * sw1_dry_t_1_population +
#           repr_log_size_x_sw1_wet  * z_1 * sw1_wet_t_1_population +
#           # Random site intercept and nativity coefficients
#           repr_intercept_population + 
#           repr_native * is_native(population) +
#           repr_log_size_x_native * z_1 * is_native(population)
#       ),
#       n_flow_population = exp(
#         flow_int + 
#           flow_z * z_1 + 
#           flow_temp_dry * temp_dry_t_1_population + 
#           flow_temp_wet * temp_wet_t_1_population +
#           flow_prec_dry * prec_dry_t_1_population +
#           flow_prec_wet * prec_wet_t_1_population +
#           flow_sw3_dry  * sw3_dry_t_1_population +
#           flow_sw3_wet  * sw3_wet_t_1_population +
#           # Interactions
#           flow_log_size_x_temp_dry * z_1 * temp_dry_t_1_population + 
#           flow_log_size_x_temp_wet * z_1 * temp_wet_t_1_population +
#           flow_log_size_x_prec_dry * z_1 * prec_dry_t_1_population +
#           flow_log_size_x_prec_wet * z_1 * prec_wet_t_1_population +
#           flow_log_size_x_sw3_dry  * z_1 * sw3_dry_t_1_population +
#           flow_log_size_x_sw3_wet  * z_1 * sw3_wet_t_1_population +
#           # Random site intercept and nativity coefficients
#           flow_intercept_population + 
#           flow_native * is_native(population) +
#           flow_log_size_x_native * z_1 * is_native(population)
#       ),
#       recr_size         = dnorm(z_2, recr_mu, recr_sd),
#       states            = list("z"),
#       data_list         = det_data_list,
#       uses_par_sets     = TRUE,
#       par_set_indices   = list(population = unique(all_ramets$site)),
#       evict_cor         = TRUE,
#       evict_fun         = truncated_distributions("norm", "recr_size")    
#     ) %>%
#     define_impl(
#       make_impl_args_list(c("P_population", "F_population"),
#                           rep("midpoint", 2),
#                           rep("z", 2),
#                           rep("z", 2))
#     ) %>%
#     define_domains(z = c(L_z, U_z, 150)) %>%
#     define_pop_state(n_z_population = runif(150)) %>%
#     make_ipm(
#       usr_funs = list(
#         inv_logit = inv_logit,
#         is_native = is_native),
#       iterations = 150
#     )
#   
#   temp_lams[[i]] <- as_tibble(t(lambda(boot_ipm))) %>%
#     mutate(obs = "no")
#   
#   if(!all(is_conv_to_asymptotic(boot_ipm, tolerance = 1e-5))) {
#     message("Iteration ", i, " has at least one unconverged model!")
#   }
# 
# }
# 
# Sys.time() - start # ROund 1 was 2.48 hours
# 
# all_lams <- bind_rows(boot_lams, temp_lams) %>%
#   as.data.frame()
# 
# write.csv(all_lams,
#           file = "ipms/Model_Fits/ipms/lin_mod_site_lambdas.csv",
#           row.names = FALSE,
#           quote = FALSE)
# 
# 
# pdf("ipms/Figures/ipm_output/lin_mod_site_lambdas.pdf",
#     height = 8,
#     width = 8)
#   
#   par(mfrow = c(2,2))
#   for(i in seq_len(ncol(all_lams) - 1)) {
#     
#     use_lam <- all_lams[ , i, drop = TRUE]
#     max_p <- ifelse(max(use_lam) > 3, 3,
#                    max(use_lam))
#     min_p <- ifelse(min(use_lam) < 0.05, 0, min(use_lam))
#     
#     
#     hist(use_lam,
#               main = names(all_lams)[i],
#               breaks = 200,
#               xlim = c(min_p, max_p))
#     abline(v = mean(use_lam), col = "red")
#     abline(v = quantile(use_lam, prob = 0.025), 
#            col = "blue", 
#            lty = 2)
#     abline(v = quantile(use_lam, prob = 0.975), 
#            col = "blue", 
#            lty = 2)
#   }
# 
# 
# dev.off()