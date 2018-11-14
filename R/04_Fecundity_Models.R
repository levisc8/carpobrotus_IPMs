# Fecundity models
# Starting with hierarchical models using population as random effect. 
# We also need to do a bit of model selection. Following the recommendations of 
# Ellner, Rees, and Childs (2016), fit a series of progressively more complicated
# models and then use a likelihood ratio test to see which is worth keeping

# This is mostly for prototyping purposes. I will likely move to a Bayesian model
# once we get a better idea of what the right form of said model should be.


# Start with ramet level

# intercept only, varies across sites
ramet_repro_int_only <- glmer(repro ~ 1 + (1|population),
                              data = all_ramets,
                              family = binomial())

# size fixed, intercept varies across sites
ramet_repro_size_int_r <- glmer(repro ~ log_size + (1|population),
                                data = all_ramets, 
                                family = binomial())

# slope and intercepts vary across sites, but are not correlated - REFIT starting
# with old output
ramet_repro_slope_int_uncor <- glmer(repro ~ log_size + 
                                             (1|population) + 
                                             (0 + log_size|population),
                                     data = all_ramets,
                                     family = binomial())

new_start <- getME(ramet_repro_slope_int_uncor, c('theta', 'fixef'))

ramet_repro_slope_int_uncor_2 <- update(ramet_repro_slope_int_uncor,
                                        start = new_start)

# slope and intercepts vary across sites and are correlated. I doubt this will
# yield additional information, but would be interesting if there was spatial 
# autocorrelation in parameters...
ramet_repro_slope_int_cor <- glmer(repro ~ log_size + (log_size|population),
                                   data = all_ramets,
                                   family = binomial())

ramet_LRT <- anova(ramet_repro_int_only,
                   ramet_repro_size_int_r,
                   ramet_repro_slope_int_uncor_2,
                   ramet_repro_slope_int_cor)

# Now genets
# intercept only, varies across sites
genet_repro_int_only <- glmer(repro ~ 1 + (1|population),
                              data = all_genets,
                              family = binomial())

# size fixed, intercept varies across sites
genet_repro_size_int_r <- glmer(repro ~ log_size + (1|population),
                                data = all_genets, 
                                family = binomial())

# slope and intercepts vary across sites, but are not correlated 
genet_repro_slope_int_uncor <- glmer(repro ~ log_size + 
                                       (1|population) + 
                                       (0 + log_size|population),
                                     data = all_genets,
                                     family = binomial())


# All correlated as above - REFIT starting with output from unconverged model
genet_repro_slope_int_cor <- glmer(repro ~ log_size + (log_size|population),
                                   data = all_genets,
                                   family = binomial())

new_start <- getME(genet_repro_slope_int_cor, c('theta', 'fixef'))

genet_repro_slope_int_cor_2 <- update(genet_repro_slope_int_cor,
                                        start = new_start)



genet_LRT <- anova(genet_repro_int_only,
                   genet_repro_size_int_r,
                   genet_repro_slope_int_uncor,
                   genet_repro_slope_int_cor_2)


sink(file = 'Model_Summaries/ramet_pr_repro_models.txt') 
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

sink(file = 'Model_Summaries/genet_pr_repro_models.txt') 
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
ramet_flower_n_int_only <- glmer.nb(flower_n ~ 1 + (1|population),
                              data = all_ramets)

# size fixed, intercept varies across sites
ramet_flower_n_size_int_r <- glmer.nb(flower_n ~ log_size + (1|population),
                                data = all_ramets)

# slope and intercepts vary across sites, but are not correlated
ramet_flower_n_slope_int_uncor <- glmer.nb(flower_n ~ log_size + 
                                       (1|population) + 
                                       (0 + log_size|population),
                                     data = all_ramets)

# slope and intercepts vary across sites and are correlated. I doubt this will
# yield additional information, but would be interesting if there was spatial 
# autocorrelation in parameters...
ramet_flower_n_slope_int_cor <- glmer.nb(flower_n ~ log_size + (log_size|population),
                                   data = all_ramets)

ramet_flower_LRT <- anova(ramet_flower_n_int_only,
                          ramet_flower_n_size_int_r,
                          ramet_flower_n_slope_int_uncor,
                          ramet_flower_n_slope_int_cor)

# Now genets
# intercept only, varies across sites
genet_flower_n_int_only <- glmer.nb(flower_n ~ 1 + (1|population),
                              data = all_genets)

# size fixed, intercept varies across sites
genet_flower_n_size_int_r <- glmer.nb(flower_n ~ log_size + (1|population),
                                data = all_genets)

# slope and intercepts vary across sites, but are not correlated
genet_flower_n_slope_int_uncor <- glmer.nb(flower_n ~ log_size + 
                                       (1|population) + 
                                       (0 + log_size|population),
                                     data = all_genets)

# All correlated as above
genet_flower_n_slope_int_cor <- glmer.nb(flower_n ~ log_size + (log_size|population),
                                   data = all_genets)


genet_flower_LRT <- anova(genet_flower_n_int_only,
                          genet_flower_n_size_int_r,
                          genet_flower_n_slope_int_uncor,
                          genet_flower_n_slope_int_cor)

sink(file = 'Model_Summaries/ramet_flower_models.txt') 
  cat('Intercept only, varies across populations\n\n *********************\n\n')
  print(summary(ramet_flower_n_int_only))
  cat('\n\n*********************\n\nSlope fixed, intercept varies across',
      'populations\n\n*********************\n\n')
  print(summary(ramet_flower_n_size_int_r))
  cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
      ' but are not correlated\n\n*********************\n\n')
  print(summary(ramet_flower_n_slope_int_uncor))
  cat('\n\n*********************\n\nSlope and intercepts vary across sites,',
      ' and are correlated\n\n*********************\n\n')
  print(summary(ramet_flower_n_slope_int_cor))
  cat('\n\n*********************\n\nLikelihood ratio tests',
      '\n\n*********************\n\n')
  print(ramet_flower_LRT)
  cat('\n\nEnd output')
sink()


sink(file = 'Model_Summaries/genet_flower_models.txt')
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

all_ramets$flower_pred <- add_predictions(ramet_flower_n_slope_int_uncor,
                                          all_ramets)


all_genets$repro_pred <- predict(genet_repro_size_int_r,
                                 type = 'response')
all_genets$flower_pred <- add_predictions(genet_flower_n_size_int_r,
                                          all_genets)

# Finish fecundity model fitting