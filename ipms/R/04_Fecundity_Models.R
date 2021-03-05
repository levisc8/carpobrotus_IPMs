# Fecundity models
# Starting with hierarchical models using populations as the sole random effect.
# initially went with populations nested within countries, but that proved difficult
# as estimating effects for Portugal/Spain + Israel became difficult with so few
# levels. I also don't think that political boundaries are making too much difference,
# as control efforts seem at least somewhat comparable across sites within each
# country. 
#
# I am skipping fitting models for genets because the more I work on this,
# the less certain I feel about our ability to actually create true genets
# from this data set. Hence, I will only fit models for ramets for now.
# Including genets shouldn't take too much extra effort if we decide to attack
# that later.
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

repro_data <- filter(all_ramets, !is.na(repro) & !is.na(log_size))
flower_data <- filter(all_ramets, !is.na(flower_n) & !is.na(log_size))

# BRMS fits ------------------
# intercept only, varies across sites

ramet_repro_int_only_brm <- brm(repro ~ 1 + (1|population),
                                data = repro_data,
                                family = bernoulli(),
                                chains = 4,
                                cores = getOption("mc.cores", 4L),
                                save_model = 'Stan/repro_int_only.stan',
                                save_dso = TRUE,
                                control = list(adapt_delta = 0.99))

# # size fixed, intercept varies across sites

ramet_repro_size_int_r_brm <- brm(repro ~ log_size + (1|population),
                                  data = repro_data,
                                  family = bernoulli(),
                                  chains = 4,
                                  cores = getOption("mc.cores", 4L),
                                  save_model = 'Stan/repro_size_fixed_int_varies.stan',
                                  save_dso = TRUE,
                                  control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites, but are not correlated

ramet_repro_slope_int_uncor_brm <- brm(repro ~ log_size +
                                         (log_size || population),
                                       data = repro_data,
                                       family = bernoulli(),
                                       chains = 4,
                                       cores = getOption("mc.cores", 4L),
                                       save_model = 'Stan/repro_size_int_varies.stan',
                                       save_dso = TRUE,
                                       control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites and are correlated.

ramet_repro_slope_int_cor_brm <- brm(repro ~ log_size + 
                                       (log_size|population),
                                     data = repro_data,
                                     family = bernoulli(),
                                     chains = 4,
                                     cores = getOption("mc.cores", 4L),
                                     save_model = 'Stan/repro_size_int_varies_cor.stan',
                                     save_dso = TRUE,
                                     control = list(adapt_delta = 0.99))

brm_repro_list <- list(int_only = ramet_repro_int_only_brm,
                       size_int_r = ramet_repro_size_int_r_brm,
                       slope_int_uncor = ramet_repro_slope_int_uncor_brm,
                       slope_int_cor = ramet_repro_slope_int_cor_brm)

brm_repro_waic <- waic(ramet_repro_int_only_brm,
                       ramet_repro_size_int_r_brm,
                       ramet_repro_slope_int_uncor_brm,
                       ramet_repro_slope_int_cor_brm)

brm_repro_list$waic <- brm_repro_waic

saveRDS(brm_repro_list, file = 'Model_Fits/ramet_repro_list_brms.rds')
# models <- readRDS('Model_Fits/ramet_repro_list_brms.rds')
# 
# 
# ramet_repro_int_only_brm <- models$int_only
# ramet_repro_size_int_r_brm <- models$size_int_r
# ramet_repro_slope_int_uncor_brm <- models$slope_int_uncor
# ramet_repro_slope_int_cor_brm <- models$slope_int_cor
# 
# brm_repro_waic <- models$waic

pdf("Model_Summaries/Trace_Plots/ramet_int_only_repro.pdf")

plot(ramet_repro_int_only_brm,
     ask = FALSE)

p <- pp_check(ramet_repro_int_only_brm,
              type     = 'bars_grouped',
              group    = 'population',
              nsamples = 100L,
              freq     = FALSE)

print(p)

dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_slope_int_random_repro.pdf')
plot(ramet_repro_size_int_r_brm,
     ask = FALSE)

p <- pp_check(ramet_repro_size_int_r_brm,
              type     = 'bars_grouped',
              group    = 'population',
              nsamples = 100L,
              freq     = FALSE)

print(p)
dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_uncor_slope_int_repro.pdf',
    width = 8,
    height = 8)

plot(ramet_repro_slope_int_uncor_brm,
     ask = FALSE)

p <- pp_check(ramet_repro_slope_int_uncor_brm,
              type     = 'bars_grouped',
              group    = 'population',
              nsamples = 100L,
              freq     = FALSE)
print(p)

dev.off()

pdf("Model_Summaries/Trace_Plots/ramet_cor_slope_int_repro.pdf",
    width = 8,
    height = 8)

plot(ramet_repro_slope_int_cor_brm,
     ask = FALSE)

p <- pp_check(ramet_repro_slope_int_cor_brm,
              type     = 'bars_grouped',
              group    = 'population',
              nsamples = 100L,
              freq     = FALSE)
print(p)
dev.off()


sink(file = 'Model_Summaries/ramet_repro_brm.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_repro_int_only_brm))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_repro_size_int_r_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(ramet_repro_slope_int_uncor_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_repro_slope_int_cor_brm))

cat('\n\n*********************\n\nWAIC Results\n\n')
print(
  brm_repro_waic
)

cat('\n\nEnd output')
sink()


# Flower production models -----------------

# # intercept only, varies across sites
ramet_flower_n_int_only_brm <- brm(flower_n ~ flower_col + (1|population),
                                   data = flower_data,
                                   family = negbinomial(),
                                   cores = getOption("mc.cores", 4L),
                                   save_model = 'Stan/flower_n_int_only.stan',
                                   save_dso = TRUE,
                                   control = list(adapt_delta = 0.99,
                                                  max_treedepth = 12))

# # size fixed, intercept varies across sites
ramet_flower_n_size_int_r_brm <- brm(flower_n ~ log_size * flower_col +
                                       (1|population),
                                     data = flower_data,
                                     family = negbinomial(),
                                     cores = getOption("mc.cores", 4L),
                                     save_model = 'Stan/flower_n_size_fixed_int_varies.stan',
                                     save_dso = TRUE,
                                     control = list(adapt_delta = 0.99,
                                                    max_treedepth = 12))

# # slope and intercepts vary across sites, but are not correlated 
ramet_flower_n_slope_int_uncor_brm <- brm(flower_n ~ log_size * flower_col +
                                            (log_size || population),
                                          data = flower_data,
                                          family = negbinomial(),
                                          cores = getOption("mc.cores", 4L),
                                          save_model = 'Stan/flower_n_size_int_varies.stan',
                                          save_dso = TRUE,
                                          control = list(adapt_delta = 0.99,
                                                         max_treedepth = 12))

# # slope and intercepts vary across sites and are correlated.
ramet_flower_n_slope_int_cor_brm <- brm(flower_n ~ log_size * flower_col + 
                                          (log_size|population),
                                        data = flower_data,
                                        family = negbinomial(),
                                        cores = getOption("mc.cores", 4L),
                                        save_model = 'Stan/flower_n_size_int_varies_cor.stan',
                                        save_dso = TRUE,
                                        control = list(adapt_delta = 0.99,
                                                       max_treedepth = 12))
brm_flower_n_list <- list(int_only = ramet_flower_n_int_only_brm,
                          size_int_r = ramet_flower_n_size_int_r_brm,
                          slope_int_uncor = ramet_flower_n_slope_int_uncor_brm,
                          slope_int_cor = ramet_flower_n_slope_int_cor_brm)

brm_flower_n_waic <- waic(ramet_flower_n_int_only_brm,
                          ramet_flower_n_size_int_r_brm,
                          ramet_flower_n_slope_int_uncor_brm,
                          ramet_flower_n_slope_int_cor_brm)

brm_flower_n_list$waic <- brm_flower_n_waic

saveRDS(brm_flower_n_list, file = 'Model_Fits/ramet_flower_list_brms.rds')

# models <- readRDS('Model_Fits/ramet_flower_list_brms.rds')
# 
# ramet_flower_n_int_only_brm <- models$int_only
# ramet_flower_n_size_int_r_brm <- models$size_int_r
# ramet_flower_n_slope_int_uncor_brm <- models$slope_int_uncor
# ramet_flower_n_slope_int_cor_brm <- models$slope_int_cor
# brm_flower_n_waic <- models$waic

pdf("Model_Summaries/Trace_Plots/ramet_int_only_flower_n.pdf")

plot(ramet_flower_n_int_only_brm,
     ask = FALSE)

p <- pp_check(ramet_flower_n_int_only_brm,
              type     = 'scatter_avg_grouped',
              group    = 'population',
              nsamples = 100L)

print(p + geom_abline(slope = 1, intercept = 0))

dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_slope_int_random_flower_n.pdf')

plot(ramet_flower_n_size_int_r_brm,
     ask = FALSE)

p <- pp_check(ramet_flower_n_size_int_r_brm,
              type     = 'scatter_avg_grouped',
              group    = 'population',
              nsamples = 100L)

print(p + geom_abline(slope = 1, intercept = 0))

dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_uncor_slope_int_flower_n.pdf',
    width = 8,
    height = 8)

plot(ramet_flower_n_slope_int_uncor_brm,
     ask = FALSE)

p <- pp_check(ramet_flower_n_slope_int_uncor_brm,
              type     = 'scatter_avg_grouped',
              group    = 'population',
              nsamples = 100L)

print(p + geom_abline(slope = 1, intercept = 0))


dev.off()

pdf("Model_Summaries/Trace_Plots/ramet_cor_slope_int_flower_n.pdf",
    width = 8,
    height = 8)

plot(ramet_flower_n_slope_int_cor_brm,
     ask = FALSE)

p <- pp_check(ramet_flower_n_slope_int_cor_brm,
              type     = 'scatter_avg_grouped',
              group    = 'population',
              nsamples = 100L)

print(p + geom_abline(slope = 1, intercept = 0))



dev.off()

sink(file = 'Model_Summaries/ramet_flower_brm.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_flower_n_int_only_brm))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_flower_n_size_int_r_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(ramet_flower_n_slope_int_uncor_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_flower_n_slope_int_cor_brm))
cat('\n\n*********************\n\nCross Validation Results\n\n')
print(
  brm_flower_n_waic
)

cat('\n\nEnd output')
sink()


# Recruit size distribution model
recruits <- filter(all_ramets, is.na(log_size) & id > 7999)

recr_size_model <- brm(log_size_next ~ 1 + (1 | population),
                       data = recruits,
                       family = gaussian(),
                       chains = 4L,
                       cores = getOption("mc.cores", 4L),
                       save_model = 'Stan/recr_size_model.stan',
                       save_dso = TRUE,
                       control = list(adapt_delta = 0.99))

saveRDS(recr_size_model, file = "Model_Fits/recr_size_brm.rds")

sink(file = 'Model_Summaries/recr_size_model.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(recr_size_model))
cat('\n\nEnd output')
sink()


# Finish fecundity model fitting
