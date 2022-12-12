# LTRE Calculations: Random Forest with complete posterior + lambdas

# all_surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
# all_flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_lin_slope_ran.rds")
# all_grow_mods <- readRDS("ipms/Model_Fits/ramet_growth_list_krig.rds")
# all_repr_mods <- readRDS("ipms/Model_Fits/ramet_repro_list_krig.rds")
# 
# 
# surv_mod <- all_surv_mods[[1]][["times_sw3_seas"]]
# flow_mod <- all_flow_mods[[1]][["times_sw3_seas"]]
# grow_mod <- all_grow_mods[[1]][["times_sw3_seas"]]
# repr_mod <- all_repr_mods[[1]][["times_sw1_seas"]]

recr_mod <- readRDS("ipms/Model_Fits/recr_size_brm.rds")
recr_draws <- as.data.frame(recr_mod) %>%
  select(1:2) %>%
  setNames(c("recr_size_int", "recr_size_sigma")) %>%
  replicate(n = 13, expr = list(.)) %>%
  do.call(what = "rbind", args = .)

surv_draws <- read.csv("ipms/Model_Fits/posterior_lambda_surv_pars.csv") %>%
  rename_draws("surv") %>%
  setNames(gsub("intercept", "int", names(.)))
names(surv_draws) <- gsub(pattern = "^surv_z$",
                          "surv_log_size",
                          x = names(surv_draws))

all_surv <- gather(data = surv_draws, surv_int_Colares:surv_int_Whirinaki, 
                     value = "surv_intercept_ran", key = "site") %>%
  select(-site)

grow_draws <- read.csv("ipms/Model_Fits/posterior_lambda_grow_pars.csv") %>%
  rename_draws("grow") %>%
  setNames(gsub("intercept", "int", names(.)))

names(grow_draws) <- gsub(pattern = "^grow_z$", 
                          "grow_log_size", 
                          x = names(grow_draws))

all_grow_1 <- gather(data = grow_draws, grow_int_Colares:grow_int_Whirinaki, 
                     value = "grow_intercept_ran", key = "site") %>%
  select(-site)

all_grow_2 <- gather(data = grow_draws, 
                     grow_sigma_int_Colares:grow_sigma_int_Whirinaki, 
                     value = "grow_sigma_intercept_ran", key = "site") %>%
  select(-site)

all_grow <- cbind(all_grow_1, 
                  grow_sigma_intercept_ran = all_grow_2$grow_sigma_intercept_ran) %>%
  select(-c(grow_sigma_int_Colares:grow_sigma_int_Whirinaki))


repr_draws <- read.csv("ipms/Model_Fits/posterior_lambda_repr_pars.csv") %>%
  rename_draws("repr") %>%
  setNames(gsub("intercept", "int", names(.)))

names(repr_draws) <- gsub(pattern = "^repr_z$", 
                          "repr_log_size", 
                          x = names(repr_draws))

all_repr <- gather(data = repr_draws, repr_int_Colares:repr_int_Whirinaki, 
                     value = "repr_intercept_ran", key = "site") %>%
  select(-site)


flow_draws <- read.csv("ipms/Model_Fits/posterior_lambda_flow_pars.csv") %>%
  rename_draws("flow") %>%
  setNames(gsub("intercept", "int", names(.)))

names(flow_draws) <- gsub(pattern = "^flow_z$", 
                          "flow_log_size", 
                          x = names(flow_draws))

names(flow_draws)[grepl(",log_size", names(flow_draws))] <- gsub(
  "int", "log_size", names(flow_draws)[grepl(",log_size", names(flow_draws))]
)

names(flow_draws)[grepl(",log_size", names(flow_draws))] <- gsub(
  ",log_size", "", names(flow_draws)[grepl(",log_size", names(flow_draws))]
)
flow_draws_1 <- gather(data = flow_draws, flow_int_Colares:flow_int_Whirinaki, 
                       value = "flow_intercept_ran", key = "site")
flow_draws_2 <- gather(data = flow_draws, flow_log_size_Colares:flow_log_size_Whirinaki, 
                       value = "flow_slope_ran", key = "site")

lambdas <- read.csv("ipms/Model_Fits/ipms/lin_mod_general_ipm_site_lambdas_prob_resamp.csv",
                    stringsAsFactors = FALSE) %>%
  filter(obs == "no") %>%
  select(-obs) %>%
  gather(key = "site", value = "lambda") %>%
  mutate(site = gsub("lambda_", "", site))


all_flow <- cbind(lambda = lambdas$lambda,
                  flow_draws_1, 
                  flow_slope_ran = flow_draws_2$flow_slope_ran) %>%
  select(-c(flow_log_size_Colares:flow_log_size_Whirinaki))

all_data <- cbind(all_flow, all_repr, all_grow, all_surv, recr_draws)
lams <- all_data$lambda
pop_draws_site <- select(all_data, -c(lambda))

split_ind <- sample(1:2, 
                    nrow(pop_draws_site), 
                    replace = TRUE,
                    prob = c(0.7, 0.3))

lams_test  <- lams[split_ind == 2]
pop_test   <- pop_draws_site[split_ind == 2, ]
lams_train <- lams[split_ind == 1]
pop_train  <- pop_draws_site[split_ind == 1, ]

# temp_rf_site <- randomForest(x          = pop_train,
#                              y          = lams_train,
#                              xtest      = pop_test,
#                              ytest      = lams_test,
#                              ntree      = 751,
#                              importance = TRUE)
# 
# saveRDS(temp_rf_site,
#         file = "ipms/Model_Fits/ltre/all_site_rf_ltre_with_site_term.rds")

all_rf_site <- readRDS("ipms/Model_Fits/ltre/all_site_rf_ltre_with_site_term.rds")

imps <- all_rf_site$importance %>%
  as.data.frame() %>%
  mutate(var = rownames(.)) %>%
  select(c(1, 3))%>% 
  cbind(imp_sd = all_rf_site$importanceSD) %>%
  setNames(c("inc_mse", "var", "imp_sd")) %>%
  arrange(desc(inc_mse)) %>%
  # slice_head(n = 15) %>%x
  mutate(var = factor(var, levels = as.character(var)),
         var = fct_reorder(var, inc_mse, .fun = max))



r2_dat <- data.frame(n_trees = 1:751,
                     var_exp = all_rf_site$rsq * 100)

r2_plt <- ggplot(r2_dat, aes(x = n_trees, y = var_exp)) + 
  geom_line(color = "red", size = 1.5) +
  theme_bw() +
  ylab("% Variance Explained") +
  xlab("# of Trees")

vi_plt <- ggplot(imps[1:20, ], aes(x = var, y = inc_mse)) + 
  geom_point(size = 2) + 
  geom_linerange(aes(ymin = pmax(0, inc_mse - 2 * imp_sd), ymax = inc_mse + 2 * imp_sd),
                 size = 1.2) +
  coord_flip() +
  theme_bw() +
  ylab("Importance") +
  xlab("Vital Rate Regression Parameter") +
  ggtitle("Including 'site' term")


vi_plt_minus_site <- ggplot(filter(imps, var != "site")[1:19, ], 
                            aes(x = var, y = inc_mse)) + 
  geom_point(size = 2) + 
  geom_linerange(aes(ymin = pmax(0, inc_mse - 2 * imp_sd), ymax = inc_mse + 2 * imp_sd),
                 size = 1.2) +
  coord_flip() +
  theme_bw() +
  ylab("Importance") +
  xlab("Vital Rate Regression Parameter") +
  ggtitle("Excluding 'site' term")

pdf("ipms/Figures/ipm_output/all_site_ltre_rf_w_site_term.pdf")

print(vi_plt / vi_plt_minus_site)

dev.off()

pdf("ipms/Figures/ipm_output/all_site_ltre_r2_rf_w_site_term.pdf")

print(r2_plt)

dev.off()
