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

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds") 

repro_data <- filter(all_ramets, !is.na(repro) & !is.na(size)) %>%
  filter(site != "Havatselet")
flower_data <- filter(all_ramets, !is.na(flower_n) & !is.na(size)) %>%
  filter(site != "Havatselet")

# BRMS fits ------------------

repro_map_list <- fit_vr_model(repro_data, "repro", "map_rec")
repro_mat_list <- fit_vr_model(repro_data, "repro", "mat_rec")
repro_p_seas_list <- fit_vr_model(repro_data, "repro", "p_seas_rec")
repro_t_co_q_list <- fit_vr_model(repro_data, "repro", "t_co_qu_rec")
repro_t_seas_list <- fit_vr_model(repro_data, "repro", "t_seas_rec")

brm_repro_list <- list(map = repro_map_list,
                       mat = repro_mat_list,
                       p_seas = repro_p_seas_list,
                       t_seas = repro_t_seas_list,
                       co_qu  = repro_t_co_q_list)

saveRDS(brm_repro_list, file = 'ipms/Model_Fits/ramet_repro_list_no_is.rds')
# brm_repro_list <- readRDS('ipms/Model_Fits/ramet_repro_list.rds')

plot_models(brm_repro_list, "repro")
plot_preds(brm_repro_list, "repro")

# Flower production models -----------------
flower_n_map_list <- fit_vr_model(flower_data, "flower_n", "map_rec")
flower_n_mat_list <- fit_vr_model(flower_data, "flower_n", "mat_rec")
flower_n_p_seas_list <- fit_vr_model(flower_data, "flower_n", "p_seas_rec")
flower_n_t_co_q_list <- fit_vr_model(flower_data, "flower_n", "t_co_qu_rec")
flower_n_t_seas_list <- fit_vr_model(flower_data, "flower_n", "t_seas_rec")

brm_flower_n_list <- list(map = flower_n_map_list,
                          mat = flower_n_mat_list,
                          p_seas = flower_n_p_seas_list,
                          t_seas = flower_n_t_seas_list,
                          co_qu  = flower_n_t_co_q_list)

saveRDS(brm_flower_n_list, file = 'ipms/Model_Fits/ramet_flower_n_list_no_is.rds')
# brm_flower_n_list <- readRDS('ipms/Model_Fits/ramet_flower_n_list.rds')

plot_models(brm_flower_n_list, "flower_n")
plot_preds(brm_flower_n_list, "flower_n")

# # Recruit size distribution model
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
# 
# 
# # Finish fecundity model fitting
