# Prototype survival and growth hierarchical models

surv_data <- filter(all_ramets, !is.na(alive) & !is.na(log_size))

# Correct a typo from data entry from Rough Island, survival is 0 or 1, cannot be greater than 1!

surv_data$alive[surv_data$alive > 1] <- 1

grow_data <- filter(surv_data, alive == 1)


# Survival -----
# Intercept only! 

ramet_surv_int_only_brm <- brm(alive ~ 1 + (1 | population),
                               data = surv_data,
                               family = bernoulli(),
                               chains = 4,
                               cores = getOption("mc.cores", 4L),
                               save_model = 'Stan/surv_int_only.stan',
                               save_dso = TRUE,
                               control = list(adapt_delta = 0.99))

# # size fixed, intercept varies across sites

ramet_surv_size_int_r_brm <- brm(alive ~ log_size + (1|population),
                                 data = surv_data,
                                 family = bernoulli(),
                                 chains = 4,
                                 cores = getOption("mc.cores", 4L),
                                 save_model = 'Stan/surv_size_fixed_int_varies.stan',
                                 save_dso = TRUE,
                                 control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites, but are not correlated 

ramet_surv_slope_int_uncor_brm <- brm(alive ~ log_size +
                                        (log_size || population),
                                      data = surv_data,
                                      family = bernoulli(),
                                      chains = 4,
                                      cores = getOption("mc.cores", 4L),
                                      save_model = 'Stan/surv_size_int_varies.stan',
                                      save_dso = TRUE,
                                      control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites and are correlated.

ramet_surv_slope_int_cor_brm <- brm(alive ~ log_size + 
                                      (log_size|population),
                                    data = surv_data,
                                    family = bernoulli(),
                                    chains = 4,
                                    cores = getOption("mc.cores", 4L),
                                    save_model = 'Stan/surv_size_int_varies_cor.stan',
                                    save_dso = TRUE,
                                    control = list(adapt_delta = 0.99))

brm_surv_list <- list(int_only = ramet_surv_int_only_brm,
                      size_int_r = ramet_surv_size_int_r_brm,
                      slope_int_uncor = ramet_surv_slope_int_uncor_brm,
                      slope_int_cor = ramet_surv_slope_int_cor_brm)

brm_surv_waic <- waic(ramet_surv_int_only_brm,
                      ramet_surv_size_int_r_brm,
                      ramet_surv_slope_int_uncor_brm,
                      ramet_surv_slope_int_cor_brm)

brm_surv_list$waic <- brm_surv_waic

saveRDS(brm_surv_list, file = 'Model_Fits/ramet_surv_list_brms.rds')

pdf("Model_Summaries/Trace_Plots/ramet_int_only_surv.pdf")
plot(ramet_surv_int_only_brm,
     ask = FALSE)
p <- pp_check(ramet_surv_int_only_brm,
              type = 'bars_grouped',
              group = 'population',
              nsamples = 100L,
              freq = FALSE)
print(p)
dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_slope_int_random_surv.pdf')
plot(ramet_surv_size_int_r_brm,
     ask = FALSE)
p <- pp_check(ramet_surv_size_int_r_brm,
              type = 'bars_grouped',
              group = 'population',
              nsamples = 100L,
              freq = FALSE)
print(p)
dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_uncor_slope_int_surv.pdf',
    width = 8,
    height = 8)
plot(ramet_surv_slope_int_uncor_brm,
     ask = FALSE)
p <- pp_check(ramet_surv_slope_int_uncor_brm,
              type = 'bars_grouped',
              group = 'population',
              nsamples = 100L,
              freq = FALSE)
print(p)
dev.off()

pdf("Model_Summaries/Trace_Plots/ramet_cor_slope_int_surv.pdf",
    width = 8,
    height = 8)
plot(ramet_surv_slope_int_cor_brm,
     ask = FALSE)
p <- pp_check(ramet_surv_slope_int_cor_brm, 
              type = 'bars_grouped',
              group = 'population',
              nsamples = 100L,
              freq = FALSE)
print(p)

dev.off()


sink(file = 'Model_Summaries/ramet_surv_brm.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_surv_int_only_brm))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_surv_size_int_r_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(ramet_surv_slope_int_uncor_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_surv_slope_int_cor_brm))

cat('\n\n*********************\n\nWAIC Results\n\n')
print(
  brm_surv_waic
)

cat('\n\nEnd output')
sink()



# Growth ----

# Prototype growth hierarchical models 

# Intercept only! 

ramet_grow_int_only_brm <- brm(log_size_next ~ 1 + (1 | population),
                               data = grow_data,
                               family = gaussian(),
                               chains = 4,
                               cores = getOption("mc.cores", 4L),
                               save_model = 'Stan/grow_int_only.stan',
                               save_dso = TRUE,
                               control = list(adapt_delta = 0.99))

# # size fixed, intercept varies across sites

ramet_grow_size_int_r_brm <- brm(log_size_next ~ log_size + 
                                   (1|population),
                                 data = grow_data,
                                 family = gaussian(),
                                 chains = 4,
                                 cores = getOption("mc.cores", 4L),
                                 save_model = 'Stan/grow_size_fixed_int_varies.stan',
                                 save_dso = TRUE,
                                 control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites, but are not correlated 

ramet_grow_slope_int_uncor_brm <- brm(log_size_next ~ log_size +
                                        (log_size || population),
                                      data = grow_data,
                                      family = gaussian(),
                                      chains = 4,
                                      cores = getOption("mc.cores", 4L),
                                      save_model = 'Stan/grow_size_int_varies.stan',
                                      save_dso = TRUE,
                                      control = list(adapt_delta = 0.99))

# # slope and intercepts vary across sites and are correlated.

ramet_grow_slope_int_cor_brm <- brm(log_size_next ~ log_size + 
                                      (log_size|population),
                                    data = grow_data,
                                    family = gaussian(),
                                    chains = 4,
                                    cores = getOption("mc.cores", 4L),
                                    save_model = 'Stan/grow_size_int_varies_cor.stan',
                                    save_dso = TRUE,
                                    control = list(adapt_delta = 0.99))

brm_grow_list <- list(int_only = ramet_grow_int_only_brm,
                      size_int_r = ramet_grow_size_int_r_brm,
                      slope_int_uncor = ramet_grow_slope_int_uncor_brm,
                      slope_int_cor = ramet_grow_slope_int_cor_brm)

brm_grow_waic <- waic(ramet_grow_int_only_brm,
                      ramet_grow_size_int_r_brm,
                      ramet_grow_slope_int_uncor_brm,
                      ramet_grow_slope_int_cor_brm)

brm_grow_list$waic <- brm_grow_waic

saveRDS(brm_grow_list, file = 'Model_Fits/ramet_grow_list_brms.rds')


pdf("Model_Summaries/Trace_Plots/ramet_int_only_grow.pdf")
plot(ramet_grow_int_only_brm,
     ask = FALSE)
p <- pp_check(ramet_grow_int_only_brm, 
              type = 'scatter_avg_grouped',
              group = 'population',
              nsamples = 100L)
print(p + geom_abline(slope = 1, intercept = 0))
dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_slope_int_random_grow.pdf')
plot(ramet_grow_size_int_r_brm,
     ask = FALSE)
p <- pp_check(ramet_grow_size_int_r_brm, 
              type = 'scatter_avg_grouped',
              group = 'population',
              nsamples = 100L)
print(p + geom_abline(slope = 1, intercept = 0))
dev.off()

pdf('Model_Summaries/Trace_Plots/ramet_uncor_slope_int_grow.pdf',
    width = 8,
    height = 8)
plot(ramet_grow_slope_int_uncor_brm,
     ask = FALSE)
p <- pp_check(ramet_grow_slope_int_uncor_brm, 
              type = 'scatter_avg_grouped',
              group = 'population',
              nsamples = 100L)
print(p + geom_abline(slope = 1, intercept = 0))
dev.off()

pdf("Model_Summaries/Trace_Plots/ramet_cor_slope_int_grow.pdf",
    width = 8,
    height = 8)
plot(ramet_grow_slope_int_cor_brm,
     ask = FALSE)
p <- pp_check(ramet_grow_slope_int_cor_brm, 
              type = 'scatter_avg_grouped',
              group = 'population',
              nsamples = 100L)
print(p + geom_abline(slope = 1, intercept = 0))
dev.off()


sink(file = 'Model_Summaries/ramet_grow_brm.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_grow_int_only_brm))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_grow_size_int_r_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(ramet_grow_slope_int_uncor_brm))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_grow_slope_int_cor_brm))

cat('\n\n*********************\n\nWAIC Results\n\n')
print(
  brm_grow_waic
)

cat('\n\nEnd output')
sink()
