# Fecundity models
# Starting with hierarchical models using populations as the sole random effect.
# initially went with populations nested within countries, but that proved difficult
# as estimating effects for Portugal/Spain + Israel became difficult with so few
# levels. I also don't think that political boundaries are making too much difference,
# as control efforts seem at least somewhat comparable across sites within each
# country. 
#
# We need to do a bit of model selection. Following the recommendations of 
# Ellner, Rees, and Childs (2016), fit a series of progressively more complicated
# models and then compare them - choosing the simplest useful one.
#
# I've switch over to a Bayesian approach and using WAIC. I think
# this is the final alteration before calling the models usable.
#
# These models are density-INdependent (for now). I am still trying to come up
# with something better than a mean-field approximation. Additionally, these 
# polygons aren't necessarily individuals - some probably are all genetically identical,
# but some are almost certainly not (despite being a single polygon).

# Optionally, read in transformed data if skipping steps 2 + 3
# all_ramets <- readRDS('Data/all_ramets_di.rds')

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds") %>%
  mutate(
    native = case_when(
      site %in% c("Melkboss",      "Vogelgat", 
                  "St_Francis",    "Struisbaai",
                  "Springfontein", "Rooisand")   ~ 1,
      TRUE ~ 0
    )
  )

repro_data <- filter(all_ramets, !is.na(repro) & !is.na(size)) #%>%
  # filter(site != "Havatselet")
flower_data <- filter(all_ramets, !is.na(flower_n) & !is.na(size)) #%>%
  # filter(site != "Havatselet")

# BRMS fits ------------------

# brm_repro_list <- list(krig_clim = fit_vr_model(repro_data, "repro"))

# saveRDS(brm_repro_list, file = 'ipms/Model_Fits/ramet_repro_list_krig_gam.rds')
# brm_repro_list <- readRDS('ipms/Model_Fits/ramet_repro_list_krig_gam.rds')

# plot_models(brm_repro_list, "repro")
# plot_preds(brm_repro_list, "repro_native", native = "yes")

# Flower production models -----------------

# brm_flower_n_list <- list(krig_clim = fit_vr_model(flower_data, "flower_n"))
# 
# saveRDS(brm_flower_n_list, file = 'ipms/Model_Fits/ramet_flower_n_list_krig_gam.rds')
# brm_flower_n_list <- readRDS('ipms/Model_Fits/ramet_flower_n_list_krig_gam.rds')
# 
# plot_models(brm_flower_n_list, "flower_n")
# plot_preds(brm_flower_n_list, "flower_n_native", native = "yes")

# # Recruit size distribution model-----------
# recruits <- readRDS("ipms/Data/seedlings.rds")
# # 
# recr_size_model <- brm(log_size_next ~ 1,
#                        data = recruits,
#                        family = gaussian(),
#                        chains = 4L,
#                        backend = "cmdstanr",
#                        cores = getOption("mc.cores", 4L),
#                        save_model = 'ipms/Stan/recr_size_model.stan',
#                        control = list(adapt_delta = 0.99))
# 
# saveRDS(recr_size_model, file = "ipms/Model_Fits/recr_size_brm.rds")
# 
# sink(file = 'ipms/Model_Summaries/recr_size_model.txt')
# cat('Intercept only, varies across populations\n\n *********************\n\n')
# print(summary(recr_size_model))
# cat('\n\nEnd output')
# sink()
# 
# pars <- summary(recr_size_model)
# 
# xx <- seq(min(recruits$log_size_next) - 1, 
#           max(recruits$log_size_next) + 1, 
#           length.out = 100)
# 
# yy <- dnorm(xx, 
#             mean = unlist(pars$fixed[1]),
#             sd   = unlist(pars$spec_pars[1]))
# 
# 
# 
# pdf("ipms/Model_Summaries/recr_size_model.pdf")
# 
# # Create histogram w/ normal density overlay
#   hist(recruits$log_size_next, freq = FALSE, ylim = c(0, 0.8))
#   lines(xx, yy, col = "red", lty = 2)
# 
# 
# dev.off()

# Test against final models from GAM/linear fits 

repr_gam_list <- readRDS('ipms/Model_Fits/ramet_repro_list_krig_gam.rds')
repr_lin_list <- readRDS('ipms/Model_Fits/ramet_repro_list_krig.rds')
# 
# repr_test_model <- brm(
#   repro ~ log_size + 
#     s(mean_temp_t_1, bs = "cs", k = 4) +
#     seas_temp_t_1 * log_size +
#     s(total_prec_t_1, bs = "cs", k = 4) + 
#     seas_prec_t_1 * log_size +
#     mean_sw2_t_1 * log_size + 
#     seas_sw2_t_1 * log_size +
#     (1 | site),
#   data = repro_data,
#   family = bernoulli(),
#   chains     = 4,    
#   backend    = "cmdstanr",
#   inits      = "random",
#   cores      = getOption("mc.cores", 4L),
#   save_model = "ipms/Stan/lin_gam_mix_repro.stan",
#   control    = list(adapt_delta = 0.999,
#                     max_treedepth = 15))

# saveRDS(repr_test_model, "ipms/Model_Fits/repro_lin_gam_mix.rds")
repr_test_model <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")


flow_gam_list <- readRDS('ipms/Model_Fits/ramet_flower_n_list_krig_gam.rds')
flow_lin_list <- readRDS('ipms/Model_Fits/ramet_flower_n_list_krig.rds')
# 
# flow_test_model <- brm(
#   flower_n ~ log_size + 
#     temp_dry_t_1 * log_size +
#     temp_wet_t_1 * log_size +
#     t2(prec_dry_t_1, log_size, bs = "cs", k = 4) + 
#     t2(prec_wet_t_1, log_size, bs = "cs", k = 4) + 
#     sw2_dry_t_1 * log_size + 
#     sw2_wet_t_1 * log_size + 
#     (1 | site),
#   data = flower_data,
#   family = negbinomial(),
#   chains     = 4,    
#   backend    = "cmdstanr",
#   inits      = "0",
#   cores      = getOption("mc.cores", 4L),
#   save_model = "ipms/Stan/lin_gam_mix_flower_n.stan",
#   control    = list(adapt_delta = 0.999,
#                     max_treedepth = 15))

# saveRDS(flow_test_model, "ipms/Model_Fits/flower_n_lin_gam_mix.rds")
flow_test_model <- readRDS("ipms/Model_Fits/flower_n_lin_gam_mix.rds")

repr_waic <- waic(repr_lin_list$krig_clim$times_sw2_ann, repr_test_model)
flow_waic <- waic(flow_test_model, flow_gam_list$krig_clim$times_sw2_seas)

repr_waic
flow_waic

# # Finish fecundity model fitting
