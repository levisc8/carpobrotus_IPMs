# Prototype survival and growth hierarchical models

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds") %>%
  mutate(
    native = case_when(
      site %in% c("Melkboss",      "Vogelgat", 
                  "St_Francis",    "Struisbaai",
                  "Springfontein", "Rooisand")   ~ 1,
      TRUE ~ 0
    )
  )

# These are ramets that we did not observe at t+1 due to differences in flight
# pattern of the drone. Thus, we do not know one way or another whether they
# survived/grew. We can use their data for flowering because that only requires
# observation at time t. However, they must be removed from the data set
# to correctly model survival + growth.

rm_ramets  <- read.csv('ipms/Data/t_2_omissions.csv',
                       stringsAsFactors = FALSE) %>%
  setNames(c("site", "id"))

surv_data  <- anti_join(all_ramets, rm_ramets, by = c("site", "id")) %>%
  filter(!is.na(alive) & !is.na(log_size))

grow_data <- readRDS("ipms/Data/growth_data.rds") %>%
  anti_join(rm_ramets, by = c("site", "id")) %>%
  filter(!is.na(size) & !is.na(size_next) & alive == 1) %>% 
  mutate(
    native = case_when(
      site %in% c("Melkboss",      "Vogelgat", 
                  "St_Francis",    "Struisbaai",
                  "Springfontein", "Rooisand")   ~ 1,
      TRUE ~ 0
    )
  )

# Survival -----

# brm_survival_gams <- list(clim_krig = fit_vr_model(surv_data, "alive",
#                                                    gam = TRUE))
# brm_survival_lins <- list(clim_krig = fit_vr_model(surv_data, "alive", 
#                                                    gam = FALSE))

# 
# saveRDS(brm_survival_gams, file = 'ipms/Model_Fits/ramet_survival_list_krig_gam.rds')
# saveRDS(brm_survival_lins, file = 'ipms/Model_Fits/ramet_survival_list_krig.rds')
brm_survival_gams <- readRDS(file = 'ipms/Model_Fits/ramet_survival_list_krig_gam.rds')
brm_survival_lins <- readRDS(file = 'ipms/Model_Fits/ramet_survival_list_krig.rds')

# plot_models(brm_survival_gams, "survival", gam = TRUE) 
# plot_models(brm_survival_lins, "survival", gam = FALSE)
# plot_preds(brm_survival_list, "survival_native", native = "yes")

all_surv_waic <- waic(brm_survival_gams[[1]][[1]],
                      brm_survival_gams[[1]][[2]],
                      brm_survival_gams[[1]][[3]],
                      brm_survival_gams[[1]][[4]],
                      brm_survival_gams[[1]][[5]],
                      brm_survival_gams[[1]][[6]],
                      brm_survival_gams[[1]][[7]],
                      brm_survival_gams[[1]][[8]],
                      brm_survival_gams[[1]][[9]],
                      brm_survival_gams[[1]][[10]],
                      brm_survival_gams[[1]][[11]],
                      brm_survival_gams[[1]][[12]],
                      brm_survival_gams[[1]][[13]],
                      brm_survival_lins[[1]][[1]],
                      brm_survival_lins[[1]][[2]],
                      brm_survival_lins[[1]][[3]],
                      brm_survival_lins[[1]][[4]],
                      brm_survival_lins[[1]][[5]],
                      brm_survival_lins[[1]][[6]],
                      brm_survival_lins[[1]][[7]],
                      brm_survival_lins[[1]][[8]],
                      brm_survival_lins[[1]][[9]],
                      brm_survival_lins[[1]][[10]],
                      brm_survival_lins[[1]][[11]],
                      brm_survival_lins[[1]][[12]],
                      brm_survival_lins[[1]][[13]],
                      model_names = c(paste0(names(brm_survival_gams[[1]][1:13]),
                                             "_gam"),
                                      paste0(names(brm_survival_lins[[1]][1:13]),
                                             "_lin")))

# Growth ----

# brm_grow_gams <- list(clim_krig = fit_vr_model(grow_data, "log_size_next", 
#                                                gam = TRUE))
# brm_grow_lins <- list(clim_krig = fit_vr_model(grow_data, "log_size_next", 
#                                                gam = FALSE))
# 
# saveRDS(brm_grow_gams, file = 'ipms/Model_Fits/ramet_growth_list_krig_gam.rds')
# saveRDS(brm_grow_lins, file = 'ipms/Model_Fits/ramet_growth_list_krig.rds')
brm_grow_list <- readRDS('ipms/Model_Fits/ramet_growth_list_krig_gam.rds')
brm_grow_list <- readRDS('ipms/Model_Fits/ramet_growth_list_krig.rds')
# plot_models(brm_grow_lins, "log_size_next", gam = FALSE)
# plot_models(brm_grow_gams, "log_size_next", gam = TRUE)


all_grow_waic <- waic(brm_grow_gams[[1]][[1]],
                      brm_grow_gams[[1]][[2]],
                      brm_grow_gams[[1]][[3]],
                      brm_grow_gams[[1]][[4]],
                      brm_grow_gams[[1]][[5]],
                      brm_grow_gams[[1]][[6]],
                      brm_grow_gams[[1]][[7]],
                      brm_grow_gams[[1]][[8]],
                      brm_grow_gams[[1]][[9]],
                      brm_grow_gams[[1]][[10]],
                      brm_grow_gams[[1]][[11]],
                      brm_grow_gams[[1]][[12]],
                      brm_grow_gams[[1]][[13]],
                      brm_grow_lins[[1]][[1]],
                      brm_grow_lins[[1]][[2]],
                      brm_grow_lins[[1]][[3]],
                      brm_grow_lins[[1]][[4]],
                      brm_grow_lins[[1]][[5]],
                      brm_grow_lins[[1]][[6]],
                      brm_grow_lins[[1]][[7]],
                      brm_grow_lins[[1]][[8]],
                      brm_grow_lins[[1]][[9]],
                      brm_grow_lins[[1]][[10]],
                      brm_grow_lins[[1]][[11]],
                      brm_grow_lins[[1]][[12]],
                      brm_grow_lins[[1]][[13]],
                      model_names = c(paste0(names(brm_grow_gams[[1]][1:13]),
                                             "_gam"),
                                      paste0(names(brm_grow_lins[[1]][1:13]),
                                             "_lin")))

all_grow_waic

# surv_gam_list <- readRDS('ipms/Model_Fits/ramet_survival_list_krig_gam.rds')
# surv_lin_list <- readRDS('ipms/Model_Fits/ramet_survival_list_krig.rds')

# surv_test_model <- brm(
#   alive ~ log_size + 
#     s(temp_dry_t, bs = "cs", k = 4) + 
#     temp_wet_t * log_size + 
#     prec_dry_t * log_size + 
#     prec_wet_t * log_size + 
#     sw2_dry_t * log_size + 
#     sw2_wet_t * log_size + 
#     sw1_dry_t * log_size + 
#     sw1_wet_t * log_size +
#     (1 | site),
#   data = surv_data,
#   family = bernoulli(),
#   chains     = 4,    
#   backend    = "cmdstanr",
#   inits      = "random",
#   cores      = getOption("mc.cores", 4L),
#   save_model = "ipms/Stan/lin_gam_mix_survival.stan",
#   control    = list(adapt_delta = 0.999,
#                     max_treedepth = 15))
#
# saveRDS(surv_test_model, "ipms/Model_Fits/surv_lin_gam_mix.rds")
# surv_test_model <- readRDS("ipms/Model_Fits/surv_lin_gam_mix.rds")
#
# surv_waic <- waic(surv_test_model,
#                   surv_lin_list$clim_krig$times_sw2_seas,
#                   surv_lin_list$clim_krig$times_sw1_seas)
# surv_waic
# 
# 
# grow_3 <- bf(
#   log_size_next ~ log_size +
#     t2(temp_wet_t, log_size, bs = "cs", k = 4) + 
#     t2(prec_dry_t, log_size, bs = "cs", k = 4) +
#     t2(sw3_dry_t, log_size, bs = "cs", k = 4) +
#     temp_dry_t * log_size +
#     prec_wet_t * log_size +
#     sw3_wet_t * log_size +
#     (1 | site),
#   sigma ~ log_size + native + temp_dry_t + prec_wet_t + (1|site)
# )
# 
# grow_2_1 <- bf(
#   log_size_next ~ log_size +
#     t2(mean_temp_t, log_size, bs = "cs", k = 4) +
#     seas_temp_t * log_size +
#     t2(total_prec_t, log_size, bs = "cs", k = 4) +
#     t2(seas_prec_t, log_size, bs = "cs", k = 4) +
#     s(mean_sw2_t, bs = "cs", k = 4) +
#     seas_sw2_t * log_size +
#     s(mean_sw1_t, bs = "cs", k = 4) +
#     seas_sw1_t * log_size +
#     (1 | site),
#   sigma ~ log_size + 
#     mean_temp_t + seas_temp_t + 
#     seas_prec_t + 
#     seas_sw1_t + seas_sw2_t +
#     (1 | site)
# )
# 
# grow_3_test_model <- brm(
#   grow_3,
#   data = grow_data,
#   family = gaussian(),
#   chains     = 4,
#   backend    = "cmdstanr",
#   init       = "random",
#   cores      = getOption("mc.cores", 4L),
#   save_model = "ipms/Stan/lin_gam_mix_growth_3.stan",
#   control    = list(adapt_delta = 0.999,
#                     max_treedepth = 15))
# 
# grow_2_1_test_model <- brm(
#   grow_2_1,
#   data       = grow_data,
#   family     = gaussian(),
#   chains     = 4,
#   backend    = "cmdstanr",
#   init       = "random",
#   cores      = getOption("mc.cores", 4L),
#   save_model = "ipms/Stan/lin_gam_mix_growth_2_1.stan",
#   control    = list(adapt_delta = 0.999,
#                     max_treedepth = 15))
# 
# saveRDS(grow_2_1_test_model, "ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
# saveRDS(grow_3_test_model, "ipms/Model_Fits/grow_3_lin_gam_mix.rds")

grow_gam_list <- readRDS('ipms/Model_Fits/ramet_growth_list_krig_gam.rds')
grow_2_1_test_model <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
grow_3_test_model <- readRDS("ipms/Model_Fits/grow_3_lin_gam_mix.rds")

grow_waic <- waic(grow_3_test_model,
                  grow_2_1_test_model,
                  grow_gam_list$clim_krig$times_sw1_ann,
                  grow_gam_list$clim_krig$times_sw2_ann,
                  grow_gam_list$clim_krig$times_sw3_seas)
 
grow_waic

