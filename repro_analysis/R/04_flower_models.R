# Fit the vital rate models using brms

all_ramets <- readRDS("repro_analysis/Data/demography/all_ramets_di.rds")

clim_data <- select(all_data, site, Sampled:p_seas_rec)

clim_ramets <- left_join(all_ramets, clim_data, by = c("population" = "site"))

nat_sites <- c("Melkboss", "Rooisand", "Springfontein",
               "St_Francis", "Struisbaai", "Vogelgat")

clim_ramets$native <- ifelse(clim_ramets$population %in% nat_sites,
                             1,
                             0)

saveRDS(clim_ramets, file = "repro_analysis/Data/demography/all_ramets_clim.rds")

# # Temp Coldest Quarter ------------- 
# 
# repro_mod_t_co <- brm(repro ~ log_size + t_co_qu_rec * native,
#                       data = clim_ramets,
#                       family = bernoulli(),
#                       chains = 4,
#                       iter = 6000,
#                       warmup = 2000,
#                       cores = getOption("mc.cores", 4L),
#                       save_model = 'stan/pr_flower_mod_t_co_qu.stan',
#                       control = list(adapt_delta = 0.99,
#                                      max_treedepth = 15))
# 
# flower_mod_t_co <- brm(flower_n ~ log_size + t_co_qu_rec * native,
#                        data = clim_ramets,
#                        family = negbinomial(),
#                        chains = 4,
#                        iter = 6000,
#                        warmup = 2000,
#                        cores = getOption("mc.cores", 4L),
#                        save_model = 'stan/flower_n_mod_t_co_qu.stan',
#                        control = list(adapt_delta = 0.99,
#                                       max_treedepth = 15))
# 
# 
# ## Mean Annual Temperature ------------------
# 
# repro_mod_mat <- brm(repro ~ log_size + mat_rec * native,
#                      data = clim_ramets,
#                      family = bernoulli(),
#                      chains = 4,
#                      iter = 6000,
#                      warmup = 2000,
#                      cores = getOption("mc.cores", 4L),
#                      save_model = 'stan/pr_flower_mod_mat.stan',
#                      control = list(adapt_delta = 0.99,
#                                     max_treedepth = 15))
# 
# flower_mod_mat <- brm(flower_n ~ log_size + mat_rec * native,
#                       data = clim_ramets,
#                       family = negbinomial(),
#                       chains = 4,
#                       iter = 6000,
#                       warmup = 2000,
#                       cores = getOption("mc.cores", 4L),
#                       save_model = 'stan/flower_n_mod_mat.stan',
#                       control = list(adapt_delta = 0.99,
#                                      max_treedepth = 15))
# 
# # Mean annual precip ----------------
# 
# repro_mod_map <- brm(repro ~ log_size + map_rec * native,
#                      data = clim_ramets,
#                      family = bernoulli(),
#                      chains = 4,
#                      iter = 6000,
#                      warmup = 2000,
#                      cores = getOption("mc.cores", 4L),
#                      save_model = 'stan/pr_flower_mod_map.stan',
#                      control = list(adapt_delta = 0.99,
#                                     max_treedepth = 15))
# 
# flower_mod_map <- brm(flower_n ~ log_size + map_rec * native,
#                       data = clim_ramets,
#                       family = negbinomial(),
#                       chains = 4,
#                       iter = 6000,
#                       warmup = 2000,
#                       cores = getOption("mc.cores", 4L),
#                       save_model = 'stan/flower_n_mod_map.stan',
#                       control = list(adapt_delta = 0.99,
#                                      max_treedepth = 15))
# 
# 
# # Save everything -------------
# 
# pdf("repro_analysis/Manuscript/Figures/repro_mod_diagnostics.pdf")
# 
# plot(repro_mod_t_co,
#      ask = FALSE)
# 
# 
# pp_check(repro_mod_t_co,
#          type     = 'bars_grouped',
#          group    = 'native',
#          nsamples = 100L,
#          freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_t_co,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# plot(repro_mod_map,
#      ask = FALSE)
# 
# p <- pp_check(repro_mod_map,
#               type     = 'bars_grouped',
#               group    = 'native',
#               nsamples = 100L,
#               freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_map,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# 
# plot(repro_mod_mat,
#      ask = FALSE)
# 
# p <- pp_check(repro_mod_mat,
#               type     = 'bars_grouped',
#               group    = 'native',
#               nsamples = 100L,
#               freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_mat,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# dev.off()
# 
# sink(file   = "repro_analysis/Manuscript/Figures/repro_mod_diagnostics.txt",
#      append = FALSE)
# 
# cat("\n\n****************** Temp Coldest Quarter ************\n\n")
# print(summary(repro_mod_t_co))
# cat("\n\n************** End Temp Coldest Quarter ************\n\n")
# 
# cat("\n\n****************** Mean Precip   *******************\n\n")
# print(summary(repro_mod_map))
# cat("\n\n************** End Mean Precip *********************\n\n")
# 
# 
# cat("\n\n****************** Mean Temp ***********************\n\n")
# print(summary(repro_mod_mat))
# cat("\n\n************** End Mean Temp ***********************\n\n")
# 
# 
# sink()
# 
# pdf("repro_analysis/Manuscript/Figures/flower_mod_diagnostics.pdf")
# 
# plot(flower_mod_t_co,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_t_co,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# plot(flower_mod_map,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_map,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# 
# plot(flower_mod_mat,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_mat,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# dev.off()
# 
# sink(file   = "repro_analysis/Manuscript/Figures/flower_mod_diagnostics.txt",
#      append = FALSE)
# 
# cat("\n\n****************** Temp Coldest Quarter ************\n\n")
# print(summary(flower_mod_t_co))
# cat("\n\n************** End Temp Coldest Quarter ************\n\n")
# 
# cat("\n\n****************** Mean Precip   *******************\n\n")
# print(summary(flower_mod_map))
# cat("\n\n************** End Mean Precip *********************\n\n")
# 
# 
# cat("\n\n****************** Mean Temp ***********************\n\n")
# print(summary(flower_mod_mat))
# cat("\n\n************** End Mean Temp ***********************\n\n")
# 
# sink()
# 
# vr_mod_list <- list(repro_mod_t_co  = repro_mod_t_co,
#                     repro_mod_map   = repro_mod_map,
#                     repro_mod_mat   = repro_mod_mat,
#                     flower_mod_t_co = flower_mod_t_co,
#                     flower_mod_map  = flower_mod_map,
#                     flower_mod_mat  = flower_mod_mat)
# 
# saveRDS(vr_mod_list, file = "repro_analysis/model_fits/vr_mod_list.rds")