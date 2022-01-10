# Prototype survival and growth hierarchical models

all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")

surv_data <- filter(all_ramets, !is.na(alive) & !is.na(log_size))

grow_data <- readRDS("ipms/Data/growth_data.rds") %>%
  filter(!is.na(size) & !is.na(size_next) & alive == 1)


# Survival -----
# Intercept only!

survival_map_list <- fit_vr_model(surv_data, "alive", "map_rec")
survival_mat_list <- fit_vr_model(surv_data, "alive", "mat_rec")
survival_p_seas_list <- fit_vr_model(surv_data, "alive", "p_seas_rec")
survival_t_co_q_list <- fit_vr_model(surv_data, "alive", "t_co_qu_rec")
survival_t_seas_list <- fit_vr_model(surv_data, "alive", "t_seas_rec")

brm_survival_list <- list(map = survival_map_list,
                          mat = survival_mat_list,
                          p_seas = survival_p_seas_list,
                          t_seas = survival_t_seas_list,
                          co_qu  = survival_t_co_q_list)

saveRDS(brm_survival_list, file = 'ipms/Model_Fits/ramet_survival_list.rds')
# brm_survival_list <- readRDS('Model_Fits/ramet_survival_list.rds')

plot_models(brm_survival_list, "survival")

# Growth ----

grow_map_list <- fit_vr_model(grow_data, "log_size_next", "map_rec")
grow_mat_list <- fit_vr_model(grow_data, "log_size_next", "mat_rec")
grow_p_seas_list <- fit_vr_model(grow_data, "log_size_next", "p_seas_rec")
grow_t_co_q_list <- fit_vr_model(grow_data, "log_size_next", "t_co_qu_rec")
grow_t_seas_list <- fit_vr_model(grow_data, "log_size_next", "t_seas_rec")

brm_grow_list <- list(map = grow_map_list,
                       mat = grow_mat_list,
                       p_seas = grow_p_seas_list,
                       t_seas = grow_t_seas_list,
                       co_qu  = grow_t_co_q_list)

saveRDS(brm_grow_list, file = 'ipms/Model_Fits/ramet_grow_list.rds')
# brm_grow_list <- readRDS('Model_Fits/ramet_grow_list.rds')

plot_models(brm_grow_list, "log_size_next")
