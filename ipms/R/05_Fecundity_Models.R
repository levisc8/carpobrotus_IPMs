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
recruits <- readRDS("ipms/Data/seedlings.rds")
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

## flower -> fruit -> seed transition parameters

flow_maturation <- read.csv("ipms/Data/seed_pods.csv",
                            stringsAsFactors = FALSE)

# TODO: re-check plant IDs that aren't actually in all_ramets. Use t_1 polygon
# files to get sizes + flower_n. We can't use merged ramets because sizes aren't
# preserved. What we really want is pods/m^2. conversion is 
# (n_mat / (tot_pods  / (0.16m^2 * n_quad))) * flower_n.
# 0.16 m^2 is the size of the quadrat I used to count seed pods in the field
# (40 cm X 40 cm)

all_sites <- read.csv('ipms/Data/All_Field_Sites.csv', 
                      stringsAsFactors = FALSE)

pops      <- all_sites[all_sites$Demography == 1 &
                         all_sites$Polygon_checked == 1 &
                         all_sites$Polygon_2019_checked == 1 &
                         all_sites$Site != "Lake_Ellsmere" &
                         all_sites$Notes == 'ready', # &
                       #all_sites$Country != 'Israel', 
                       c('Country', 'Site')] %>%
  arrange(Site)

t_1_temp <- list()

sf_use_s2(FALSE)

for(i in pops$Site) {

  pop <- i
  country <- pops$Country[pops$Site == i]
  
  t_1_temp[[i]] <- t_1_data_for_seed_pods(
        glue('ipms/Polygons/{country}/{pop}/Polygon_Final.shp'),
        pop,
        # Variables to keep
        population, 
        id,
        size, flower_n)
  
  st_geometry(t_1_temp[[i]]) <- NULL
}

t_1_data <- bind_rows(t_1_temp)

t_1_data$f_p_z <- t_1_data$flower_n / exp(t_1_data$log_size)

mat_data <- flow_maturation %>%
  left_join(t_1_data, by = c("population" = "population", "plant" = "id")) %>%
  filter(population != "Lake_Ellsmere") %>%
  select(population:n_ab, flower_n, log_size, f_p_z) %>%
  group_by(population, plant) %>%
  summarise(
    # number of quads per plant
    n_quads = max(quad),
    
    # total area sampled per plant
    tot_a   = n_quads * 0.16,
    
    # Total mature pods sampled
    n_p_p   = sum(n_mat),  
    
    # Pods/m^2 sampled
    est_p   = n_p_p / (tot_a),
    
    # Flowers/m^2 observed
    f_p_z   = unique(f_p_z),
    
    # pods per unit area / flowers per unit area (gives expected value of total
    # pods per plant). This can be greater than 1 because we did not necessarily
    # observe *all* flowers with the drone - some could have emerged later or 
    # earlier in the season. NB: Maybe use as a measure for how early/late
    # we were in the season to capture flowering??
    est_mat = est_p / f_p_z) %>%
  ungroup() %>%
  group_by(population) %>%
  summarise(pop_mat_factor = mean(est_mat, na.rm = TRUE)) %>%
  bind_rows(data.frame(population = "Havatselet",
                       pop_mat_factor = mean(.$pop_mat_factor, na.rm = TRUE)))

p_m_pars <- as.list(mat_data$pop_mat_factor) %>%
  setNames(nm = paste0("p_m_", mat_data$population)) 

saveRDS(p_m_pars, "ipms/Model_Fits/discrete_pars/pod_maturation_pars.rds")

# D'Antonio (1990) finds the following germination rates : Seeds in intact fruits
# extracted < 1 year after maturation: 0.24 Seeds in intact fruits extracted > 1
# year after maturation: 0.73 Talk to TMK + RSG about how to include animal
# feces data from this paper.
# https://www.jstor.org/stable/pdf/2404312.pdf?casa_token=EJiz83s5YREAAAAA:NvwT98rakddeZTiXki8r90zIHWabrJUqPh1M94UeYDQITufnLq4PbS0cGXcZNKYcl00xgPA1dpTxfwkxIU766O4oo2ke_0Zsib5I-7Yy7ZpQtoouYjU
# Vila + D'Antonio find variable viability across different sites, so I'll take
# mean of those values (NB: this is NOT the mean of values across "species". I 
# don't think the species delination they used back then still applies (Novoa 
# et al. forthcoming)).

seed_pars <- list(
  v_s  = mean(c(43.33, 80.66, 75.33)) / 100,
  g_i  = 0.24,
  g_sb = 0.73
)

# Vila & D'Antonio (1998) find seeds/fruit varies across species/hybrids. They 
# are: 
# C. edulis:    1573 +/- 65
# Hybrid:       704  +/- 43
# C. chilensis: 483  +/- 27
# Vila + D'Antonio find variable viability across different "species", so I'll
# take mean of those values (NB: I don't think the species delination they used
# back then still applies (Novoa et al. forthcoming)).
# source for table above: https://www.jstor.org/stable/pdf/176600.pdf?casa_token=axvw7Gbp014AAAAA:1dL3VbuAI5vP9zMaElpBpX48BhzHtcDpkS8Etvl-1mCbPzTrEZRbvixkfZjtgh3HLYsr5zFWPUJqCBf2QmD9ZQByiJZOul9HsCnlfPW7-JO3mpQnl1w

seed_pars <- c(seed_pars, list(r_f = mean(c(1573, 704, 483))))

saveRDS(seed_pars, "ipms/Model_Fits/discrete_pars/seed_pars.rds")
# Seedling survival:
# This is taken from two plots at Rooisand. I counted seedlings in 40x40 cm
# plots and then came back a year later and used the drone to determine number
# of survivors and how big they were. Recruit size is modeled above with an
# intercept only regression, so here is the survival. It'll just be 
# n_alive / n_tot. I couldn't find plot 2, so we need to exclude that from
# the calculation of n_tot.

all_sdls <- read.csv("ipms/Data/sdl_plots.csv",
                     stringsAsFactors = FALSE) %>%
  .[-2, ]

sdl_s    <- nrow(recruits) / sum(all_sdls$n_sdl)
  
# # Finish fecundity model fitting
