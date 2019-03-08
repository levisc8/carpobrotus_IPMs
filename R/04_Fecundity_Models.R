# Fecundity models
# Starting with hierarchical models using populations nested within countries as random effects. 
# We also need to do a bit of model selection. Following the recommendations of 
# Ellner, Rees, and Childs (2016), fit a series of progressively more complicated
# models and then use a likelihood ratio test to see which is worth keeping

# At the bottom, I switch over to a Bayesian approach. The complicated random
# effects structure IS appropriate, but currently having issues with convergence/
# singular fits in lme4. 
# https://stats.stackexchange.com/questions/378939/dealing-with-singular-fit-in-mixed-models
# suggests switching over to a Bayesian model and then inspecting trace plots for 
# parameters to see why the lme fit may not be converging. 


# Start with ramet level

# intercept only, varies across sites
ramet_repro_int_only <- glmer(repro ~ 1 + (1|population/country),
                              data = all_ramets,
                              family = binomial())

# size fixed, intercept varies across sites
ramet_repro_size_int_r <- glmer(repro ~ log_size + (1|population/country),
                                data = all_ramets, 
                                family = binomial())

# slope and intercepts vary across sites, but are not correlated - REFIT starting
# with old output
ramet_repro_slope_int_uncor <- glmer(repro ~ log_size + 
                                       (1|population/country) + 
                                       (0 + log_size|population/country),
                                     data = all_ramets,
                                     family = binomial())

new_start <- getME(ramet_repro_slope_int_uncor, c('theta', 'fixef'))

ramet_repro_slope_int_uncor_2 <- update(ramet_repro_slope_int_uncor,
                                        start = new_start)

# slope and intercepts vary across sites and are correlated. I doubt this will
# yield additional information, but would be interesting if there was spatial 
# autocorrelation in parameters...
ramet_repro_slope_int_cor <- glmer(repro ~ log_size + (log_size|population/country),
                                   data = all_ramets,
                                   family = binomial())

ramet_LRT <- anova(ramet_repro_int_only,
                   ramet_repro_size_int_r,
                   ramet_repro_slope_int_uncor_2,
                   ramet_repro_slope_int_cor)

# Now genets
# intercept only, varies across sites
genet_repro_int_only <- glmer(repro ~ 1 + (1|population/country),
                              data = all_genets,
                              family = binomial())

# size fixed, intercept varies across sites
genet_repro_size_int_r <- glmer(repro ~ log_size + (1|population/country),
                                data = all_genets, 
                                family = binomial())

# slope and intercepts vary across sites, but are not correlated 
genet_repro_slope_int_uncor <- glmer(repro ~ log_size + 
                                       (1|population/country) + 
                                       (0 + log_size|population/country),
                                     data = all_genets,
                                     family = binomial())

new_start <- getME(genet_repro_slope_int_uncor, c('theta', 'fixef'))

genet_repro_slope_int_uncor_2 <- update(genet_repro_slope_int_uncor,
                                        start = new_start)


# All correlated as above - REFIT starting with output from unconverged model
genet_repro_slope_int_cor <- glmer(repro ~ log_size + (log_size|population/country),
                                   data = all_genets,
                                   family = binomial())

new_start <- getME(genet_repro_slope_int_cor, c('theta', 'fixef'))

genet_repro_slope_int_cor_2 <- update(genet_repro_slope_int_cor,
                                      start = new_start)



genet_LRT <- anova(genet_repro_int_only,
                   genet_repro_size_int_r,
                   genet_repro_slope_int_uncor_2,
                   genet_repro_slope_int_cor_2)


sink(file = 'Model_Summaries/ramet_pr_repro_lme.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_repro_int_only))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_repro_size_int_r))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(ramet_repro_slope_int_uncor_2))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_repro_slope_int_cor))

cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')

print(ramet_LRT)

cat('\n\nEnd output')
sink()

sink(file = 'Model_Summaries/genet_pr_repro_lme.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(genet_repro_int_only))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(genet_repro_size_int_r))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(genet_repro_slope_int_uncor))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(genet_repro_slope_int_cor_2))
cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')
print(genet_LRT)
cat('\n\nEnd output')
sink()
# Fit some flower production models

# Start with ramet level

# intercept only, varies across sites
ramet_flower_n_int_only <- glmer.nb(flower_n ~ 1 + (1|population/country),
                                    data = all_ramets)

# size fixed, intercept varies across sites
ramet_flower_n_size_int_r <- glmer.nb(flower_n ~ log_size + (1|population/country),
                                      data = all_ramets)

# slope and intercepts vary across sites, but are not correlated. This model
# will not converge with the new zealand data added in, so I'm leaving it out
# for now. I don't think it will make much difference anyway. 

# ramet_flower_n_slope_int_uncor <- glmer.nb(flower_n ~ log_size + 
#                                        (1|population/country) + 
#                                        (0 + log_size|population/country),
#                                      data = all_ramets,
#                                      nb.control = glmerControl(optCtrl = list(maxfun = 1e6)))

# slope and intercepts vary across sites and are correlated. I doubt this will
# yield additional information, but would be interesting if there was spatial 
# autocorrelation in parameters...
ramet_flower_n_slope_int_cor <- glmer.nb(flower_n ~ log_size + (log_size|population/country),
                                         data = all_ramets)

ramet_flower_LRT <- anova(ramet_flower_n_int_only,
                          ramet_flower_n_size_int_r,
                          # ramet_flower_n_slope_int_uncor,
                          ramet_flower_n_slope_int_cor)

# Now genets
# intercept only, varies across sites
genet_flower_n_int_only <- glmer.nb(flower_n ~ 1 + (1|population/country),
                                    data = all_genets)

# size fixed, intercept varies across sites
genet_flower_n_size_int_r <- glmer.nb(flower_n ~ log_size + (1|population/country),
                                      data = all_genets)

# slope and intercepts vary across sites, but are not correlated
genet_flower_n_slope_int_uncor <- glmer.nb(flower_n ~ log_size + 
                                             (1|population/country) + 
                                             (0 + log_size|population/country),
                                           data = all_genets)

# All correlated as above
genet_flower_n_slope_int_cor <- glmer.nb(flower_n ~ log_size + (log_size|population/country),
                                         data = all_genets)


genet_flower_LRT <- anova(genet_flower_n_int_only,
                          genet_flower_n_size_int_r,
                          genet_flower_n_slope_int_uncor,
                          genet_flower_n_slope_int_cor)

sink(file = 'Model_Summaries/ramet_flower_lme.txt') 
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(ramet_flower_n_int_only))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(ramet_flower_n_size_int_r))
# cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
#     ' but are not correlated\n\n*********************\n\n')
# print(summary(ramet_flower_n_slope_int_uncor))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(ramet_flower_n_slope_int_cor))
cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')
print(ramet_flower_LRT)
cat('\n\nEnd output')
sink()


sink(file = 'Model_Summaries/genet_flower_lme.txt')
cat('Intercept only, varies across populations\n\n *********************\n\n')
print(summary(genet_flower_n_int_only))
cat('\n\n*********************\n\nSlope fixed, intercept varies across',
    'populations\n\n*********************\n\n')
print(summary(genet_flower_n_size_int_r))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' but are not correlated\n\n*********************\n\n')
print(summary(genet_flower_n_slope_int_uncor))
cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
    ' and are correlated\n\n*********************\n\n')
print(summary(genet_flower_n_slope_int_cor))
cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')
print(genet_flower_LRT)
cat('\n\nEnd output')
sink()

# Next, add predicted values so that predicted means can be plotted. Only using
# the best models based on likelihood ratio tests

all_ramets$repro_pred <- predict(ramet_repro_size_int_r,
                                 type = 'response')

all_ramets$flower_pred <- predict(ramet_flower_n_slope_int_cor,
                                  newdata = all_ramets,
                                  type = 'response')

all_ramets$flower_pred[is.na(all_ramets$flower_n)] <- NA_real_


all_genets$repro_pred <- predict(genet_repro_size_int_r,
                                 type = 'response')

all_genets$flower_pred <- predict(genet_flower_n_size_int_r,
                                  newdata = all_genets,
                                  type = 'response')

all_genets$flower_pred[is.na(all_genets$flower_n)] <- NA_real_


# BRMS fits ------------------
# intercept only, varies across sites
ramet_repro_int_only_brm <- brm(repro ~ 1 + (1|population/country),
                                data = all_ramets,
                                family = bernoulli(),
                                cores = getOption("mc.cores", 4L),
                                save_model = 'Stan/int_only.stan',
                                save_dso = TRUE)

# size fixed, intercept varies across sites
ramet_repro_size_int_r_brm <- brm(repro ~ log_size + (1|population/country),
                                  data = all_ramets,
                                  family = bernoulli(),
                                  cores = getOption("mc.cores", 4L),
                                  save_model = 'Stan/size_fixed_int_varies.stan',
                                  save_dso = TRUE)

# slope and intercepts vary across sites, but are not correlated 
ramet_repro_slope_int_uncor_brm <- brm(repro ~ log_size + 
                                         (1|population/country) + 
                                         (0 + log_size|population/country),
                                       data = all_ramets,
                                       family = bernoulli(),
                                       cores = getOption("mc.cores", 4L),
                                       save_model = 'Stan/size_int_varies.stan',
                                       save_dso = TRUE)

# slope and intercepts vary across sites and are correlated.
ramet_repro_slope_int_cor_brm <- brm(repro ~ log_size + (log_size|population/country),
                                     data = all_ramets,
                                     family = bernoulli(),
                                     cores = getOption("mc.cores", 4L),
                                     save_model = 'Stan/size_int_varies_cor.stan',
                                     save_dso = TRUE)


plot(ramet_repro_int_only_brm)
ggsave('ramet_int_only_repro.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_repro_size_int_r_brm)
ggsave('ramet_slope_int_random_repro.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_repro_slope_int_uncor_brm )
ggsave('ramet_uncor_slope_int_repro.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_repro_slope_int_cor_brm )
ggsave('ramet_cor_slope_int_repro.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)


# intercept only, varies across sites
ramet_flower_n_int_only_brm <- brm(flower_n ~ 1 + (1|population/country),
                                   data = all_ramets,
                                   family = negbinomial(),
                                   cores = getOption("mc.cores", 4L),
                                   save_model = 'Stan/int_only.stan',
                                   save_dso = TRUE)

# size fixed, intercept varies across sites
ramet_flower_n_size_int_r_brm <- brm(flower_n ~ log_size + (1|population/country),
                                     data = all_ramets,
                                     family = negbinomial(),
                                     cores = getOption("mc.cores", 4L),
                                     save_model = 'Stan/size_fixed_int_varies.stan',
                                     save_dso = TRUE)

# slope and intercepts vary across sites, but are not correlated 
ramet_flower_n_slope_int_uncor_brm <- brm(flower_n ~ log_size + 
                                            (1|population/country) + 
                                            (0 + log_size|population/country),
                                          data = all_ramets,
                                          family = negbinomial(),
                                          cores = getOption("mc.cores", 4L),
                                          save_model = 'Stan/size_int_varies.stan',
                                          save_dso = TRUE)

# slope and intercepts vary across sites and are correlated.
ramet_flower_n_slope_int_cor_brm <- brm(flower_n ~ log_size + (log_size|population/country),
                                        data = all_ramets,
                                        family = negbinomial(),
                                        cores = getOption("mc.cores", 4L),
                                        save_model = 'Stan/size_int_varies_cor.stan',
                                        save_dso = TRUE)


plot(ramet_flower_n_int_only_brm)
ggsave('ramet_int_only_flower_n.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_flower_n_size_int_r_brm)
ggsave('ramet_slope_int_random_flower_n.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_flower_n_slope_int_uncor_brm )
ggsave('ramet_uncor_slope_int_flower_n.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)

plot(ramet_flower_n_slope_int_cor_brm )
ggsave('ramet_cor_slope_int_flower_n.pdf',
       path = "Model_Summaries/Trace_Plots",
       height = 8,
       width = 8,
       unit = "in",
       dpi = 400)


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
cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')
cat('\n\nEnd output')
sink()


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
cat('\n\n*********************\n\nLikelihood ratio tests',
    '\n\n*********************\n\n')
cat('\n\nEnd output')
sink()

# Finish fecundity model fitting