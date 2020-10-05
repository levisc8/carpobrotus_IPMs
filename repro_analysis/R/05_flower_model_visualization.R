# Load fitted model objects

vr_mod_list <- readRDS("repro_analysis/model_fits/vr_mod_list.rds")

clim_ramets <- readRDS("repro_analysis/Data/demography/all_ramets_clim.rds")

size_seq <- seq(min(clim_ramets$log_size),
                max(clim_ramets$log_size),
                length.out = 50) 
# 
# map_seq <- seq(min(clim_ramets$map_rec),
#                max(clim_ramets$map_rec),
#                length.out = 50)
# 
# mat_seq <- seq(min(clim_ramets$mat_rec),
#                max(clim_ramets$mat_rec),
#                length.out = 50) 
# 
# t_c_seq <- seq(min(clim_ramets$t_co_qu_rec),
#                max(clim_ramets$t_co_qu_rec),
#                length.out = 50) 
# 
# t_co_pred_dat <- expand.grid(log_size    = size_seq,
#                              t_co_qu_rec = t_c_seq,
#                              native      = c(1, 0))
# 
# mat_pred_dat <- expand.grid(log_size    = size_seq,
#                             mat_rec     = mat_seq,
#                             native      = c(1, 0))
# 
# map_pred_dat <- expand.grid(log_size    = size_seq,
#                             map_rec     = map_seq,
#                             native      = c(1, 0))
# 
# pred_repro_t_co_qu  <- predict(vr_mod_list$repro_mod_t_co,
#                                newdata = t_co_pred_dat) %>%
#   as_tibble() %>%
#   cbind(t_co_pred_dat) %>% 
#   mutate(clim_var = "t_co_qu_rec",
#          clim_val = t_co_qu_rec) %>%
#   select(-t_co_qu_rec)
# 
# pred_repro_mat      <- predict(vr_mod_list$repro_mod_mat,
#                                newdata = mat_pred_dat) %>%
#   as_tibble() %>%
#   cbind(mat_pred_dat) %>% 
#   mutate(clim_var = "mat_rec",
#          clim_val = mat_rec) %>%
#   select(-mat_rec)
# 
# all_repro_pred      <- predict(vr_mod_list$repro_mod_map,
#                                newdata = map_pred_dat) %>%
#   as_tibble() %>%
#   cbind(map_pred_dat)%>% 
#   mutate(clim_var = "map_rec",
#          clim_val = map_rec) %>%
#   select(-map_rec) %>%
#   rbind(pred_repro_mat, pred_repro_t_co_qu)%>%
#   mutate(native = as.factor(native)) 
# 
# 
# pred_flower_t_co_qu  <- predict(vr_mod_list$flower_mod_t_co,
#                                newdata = t_co_pred_dat) %>%
#   as_tibble() %>%
#   cbind(t_co_pred_dat) %>% 
#   mutate(clim_var = "t_co_qu_rec",
#          clim_val = t_co_qu_rec) %>%
#   select(-t_co_qu_rec)
# 
# pred_flower_mat      <- predict(vr_mod_list$flower_mod_mat,
#                                newdata = mat_pred_dat) %>%
#   as_tibble() %>%
#   cbind(mat_pred_dat) %>% 
#   mutate(clim_var = "mat_rec",
#          clim_val = mat_rec) %>%
#   select(-mat_rec)
# 
# all_flower_pred      <- predict(vr_mod_list$flower_mod_map,
#                                newdata = map_pred_dat) %>%
#   as_tibble() %>%
#   cbind(map_pred_dat)%>% 
#   mutate(clim_var = "map_rec",
#          clim_val = map_rec) %>%
#   select(-map_rec)  %>%
#   rbind(pred_flower_mat, pred_flower_t_co_qu) %>%
#   mutate(native = as.factor(native))
# 
# model_preds <- list(flower_preds = all_flower_pred,
#                     repro_preds = all_repro_pred)
# 
# saveRDS(model_preds, 
#         file = "repro_analysis/Data/demography/model_predictions.rds")

model_preds <- readRDS("repro_analysis/Data/demography/model_predictions.rds")

# For plotting, we just want to see the effects of climate x nativity, not size.
# we have to choose a target size to look at for this. For now, using max_size,
# but can be adjusted using the target_size variable.

target_size <- max(size_seq)

all_flower_pred <- model_preds$flower_preds %>%
  filter(log_size == target_size)

all_repro_pred  <- model_preds$repro_preds %>%
  filter(log_size == target_size)


for_plot <- clim_ramets %>%
  select(-c(flower_col, Sampled, size,
            t_seas_rec, p_seas_rec, t_co_mo_rec)) %>%
  pivot_longer(-c(id:id_pop, native),
               names_to = "clim_var",
               values_to = "clim_val") %>%
  filter(log_size > mean(log_size))

f_plt <- ggplot(for_plot, aes(y = flower_n, x = clim_val)) +
  geom_point(aes(color = as.factor(native))) +
  geom_line(data = all_flower_pred,
            aes(y = Estimate,
                x = clim_val,
                color = as.factor(native)), 
            show.legend = FALSE,
            size = 1.5,
            linetype = "dashed") + 
  facet_wrap(~clim_var,
             scales = "free_x") +
  theme_bw() +
  ylim(c(0, max(for_plot$flower_n, na.rm = TRUE) + 100))


pdf("repro_analysis/Manuscript/Figures/flower_mod_viz.pdf")

  print(f_plt)

dev.off()

r_plt <- ggplot(for_plot, 
       aes(y = repro, x = clim_val)) +
  geom_jitter(aes(color = as.factor(native)),
              height = 0.05,
              width  = 0) +
  geom_line(data = all_repro_pred,
            aes(y = Estimate,
                x = clim_val,
                color = as.factor(native)), 
            show.legend = FALSE,
            size = 1.5,
            linetype = "dashed") + 
  facet_wrap(~clim_var,
             scales = "free_x") +
  theme_bw() +
  ylim(c(-0.1, 1.1))

pdf("repro_analysis/Manuscript/Figures/repro_mod_viz.pdf")

  print(r_plt)

dev.off()
