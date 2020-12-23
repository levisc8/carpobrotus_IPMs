# Load fitted model objects

repro_mod_list <- readRDS("repro_analysis/model_fits/vr_mod_list_centered.rds") %>%
  .[grepl("repro", names(.))]

flower_mod_list <- readRDS("repro_analysis/model_fits/vr_mod_list_centered_size_x_nat.rds") %>%
  .[grepl("flower", names(.))]

vr_mod_list <- c(repro_mod_list, flower_mod_list)

clim_ramets <- readRDS("repro_analysis/Data/demography/all_ramets_clim.rds")

clim_ramets$map_rec     <- scale(clim_ramets$map_rec)
clim_ramets$t_co_qu_rec <- scale(clim_ramets$t_co_qu_rec)
clim_ramets$mat_rec     <- scale(clim_ramets$mat_rec)

size_seq <- seq(min(clim_ramets$log_size),
                max(clim_ramets$log_size),
                length.out = 50) 

map_seq <- seq(min(clim_ramets$map_rec),
               max(clim_ramets$map_rec),
               length.out = 50)

mat_seq <- seq(min(clim_ramets$mat_rec),
               max(clim_ramets$mat_rec),
               length.out = 50)

t_c_seq <- seq(min(clim_ramets$t_co_qu_rec),
               max(clim_ramets$t_co_qu_rec),
               length.out = 50)

t_co_pred_dat <- expand.grid(log_size    = size_seq,
                             t_co_qu_rec = t_c_seq,
                             native      = c(1, 0))

mat_pred_dat <- expand.grid(log_size    = size_seq,
                            mat_rec     = mat_seq,
                            native      = c(1, 0))

map_pred_dat <- expand.grid(log_size    = size_seq,
                            map_rec     = map_seq,
                            native      = c(1, 0))

pred_repro_t_co_qu  <- predict(vr_mod_list$repro_mod_t_co,
                               newdata = t_co_pred_dat) %>%
  as_tibble() %>%
  cbind(t_co_pred_dat) %>%
  mutate(clim_var = "t_co_qu_rec",
         clim_val = t_co_qu_rec) %>%
  select(-t_co_qu_rec)

pred_repro_mat      <- predict(vr_mod_list$repro_mod_mat,
                               newdata = mat_pred_dat) %>%
  as_tibble() %>%
  cbind(mat_pred_dat) %>%
  mutate(clim_var = "mat_rec",
         clim_val = mat_rec) %>%
  select(-mat_rec)

all_repro_pred      <- predict(vr_mod_list$repro_mod_map,
                               newdata = map_pred_dat) %>%
  as_tibble() %>%
  cbind(map_pred_dat)%>%
  mutate(clim_var = "map_rec",
         clim_val = map_rec) %>%
  select(-map_rec) %>%
  rbind(pred_repro_mat, pred_repro_t_co_qu)%>%
  mutate(native = as.factor(native))


pred_flower_t_co_qu  <- predict(vr_mod_list$flower_mod_t_co,
                               newdata = t_co_pred_dat) %>%
  as_tibble() %>%
  cbind(t_co_pred_dat) %>%
  mutate(clim_var = "t_co_qu_rec",
         clim_val = t_co_qu_rec) %>%
  select(-t_co_qu_rec)

pred_flower_mat      <- predict(vr_mod_list$flower_mod_mat,
                               newdata = mat_pred_dat) %>%
  as_tibble() %>%
  cbind(mat_pred_dat) %>%
  mutate(clim_var = "mat_rec",
         clim_val = mat_rec) %>%
  select(-mat_rec)

all_flower_pred      <- predict(vr_mod_list$flower_mod_map,
                               newdata = map_pred_dat) %>%
  as_tibble() %>%
  cbind(map_pred_dat)%>%
  mutate(clim_var = "map_rec",
         clim_val = map_rec) %>%
  select(-map_rec)  %>%
  rbind(pred_flower_mat, pred_flower_t_co_qu) %>%
  mutate(native = as.factor(native))
#
model_preds <- list(flower_preds = all_flower_pred,
                    repro_preds = all_repro_pred)

saveRDS(model_preds,
        file = "repro_analysis/Data/demography/model_predictions.rds")

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

# Compute means and variances

mean_f_n <- function(p_r_z, f_n_z) {
  
  return(p_r_z * f_n_z)
  
}

var_f_n <- function(p_r_z, f_n_z) {
  
  omega_bar_z <- mean_f_n(p_r_z, f_n_z)
  
  out <- omega_bar_z + ((1 - p_r_z) / p_r_z) * omega_bar_z ^ 2
  
  return(out)
  
}

z <- unique(model_preds$flower_preds$log_size)
d_z <- (z[2] - z[1])

# Make sure we are combining predictions correctly
if(identical(model_preds$flower_preds[ , 5:8],
             model_preds$repro_preds[ , 5:8])) {
  
  p_r_z <- model_preds$repro_preds$Estimate
  f_n_z <- model_preds$flower_preds$Estimate
  
}

omega_z     <- mean_f_n(p_r_z, f_n_z)
var_omega_z <- var_f_n(p_r_z, f_n_z)



sum_plot_data <- cbind(model_preds$flower_preds,
                       omega_z     = omega_z,
                       var_omega_z = var_omega_z) %>%
  cbind(expand.grid(dummy_size = 1:50,
                    dummy_clim = 1:50)) %>%
  filter(omega_z > 0) %>%
  select(log_size, native:var_omega_z, dummy_size, dummy_clim)

sum_plot_data$native <- as.character(sum_plot_data$native)

sum_plot_data$native[sum_plot_data$native == "1"] <- "Native"
sum_plot_data$native[sum_plot_data$native == "0"] <- "Non-Native"

nat_map <- filter(sum_plot_data,
                  native == "Native" & clim_var == "map_rec" & omega_z > 0)
inv_map <- filter(sum_plot_data,
                  native == "Non-Native" & clim_var == "map_rec")

nat_mat <- filter(sum_plot_data,
                  native == "Native" & clim_var == "mat_rec")
inv_mat <- filter(sum_plot_data,
                  native == "Non-Native" & clim_var == "mat_rec")

nat_tco <- filter(sum_plot_data,
                  native == "Native" & clim_var == "t_co_qu_rec")
inv_tco <- filter(sum_plot_data,
                  native == "Non-Native" & clim_var == "t_co_qu_rec")

nat_map_plot <- gg_image_plot(nat_map, 
              x_vals = range(nat_map$log_size),
              y_vals = range(nat_map$clim_val),
              omega_z  ^ (1/3),
              dummy_size,
              dummy_clim) +
  ggtitle(label = "Native", subtitle = "Annual Precipitation") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

nat_mat_plot <- gg_image_plot(nat_mat, 
                              x_vals = range(nat_mat$log_size),
                              y_vals = range(nat_mat$clim_val),
                              omega_z ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Annual Precipitation") +
  theme(legend.position = "none")

nat_tco_plot <- gg_image_plot(nat_tco, 
                              x_vals = range(nat_mat$log_size),
                              y_vals = range(nat_mat$clim_val),
                              omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Temperature Coldest Quarter") +
  theme(legend.position = "none",
        axis.title.y = element_blank())



inv_map_plot <- gg_image_plot(inv_map, 
                              x_vals = range(inv_map$log_size),
                              y_vals = range(inv_map$clim_val),
                              omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "Non-Native", subtitle = "Annual Precipitation") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

inv_mat_plot <- gg_image_plot(inv_mat, 
                              x_vals = range(inv_mat$log_size),
                              y_vals = range(inv_mat$clim_val),
                              omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Annual Precipitation") +
  theme(legend.position = "bottom") + 
  scale_fill_gradient("Mean Flower Production (cube root)",
                      low = "red",
                      high = "yellow") +
  guides(fill = guide_colorbar(title.position = "top"))

inv_tco_plot <- gg_image_plot(inv_tco, 
                              x_vals = range(inv_mat$log_size),
                              y_vals = range(inv_mat$clim_val),
                              omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Temperature Coldest Quarter") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

leg <- get_gg_legend(inv_mat_plot)


pdf(file = "repro_analysis/Manuscript/Figures/mean_flower_n_clim.pdf",
    height = 10,
    width = 10)

  grid.arrange(
    nat_map_plot, nat_mat_plot, nat_tco_plot,
    inv_map_plot, inv_mat_plot + theme(legend.position = "none"), inv_tco_plot,
    leg,
    layout_matrix = matrix(
      c(1:6, NA_real_, 7, NA_real_),
      nrow = 3, ncol = 3, byrow = TRUE
    ),
    heights = c(3, 3, 1)
  )

dev.off()


# Variance in omega_bar

nat_map_plot <- gg_image_plot(nat_map, 
                              x_vals = range(nat_map$log_size),
                              y_vals = range(nat_map$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "Native", subtitle = "Annual Precipitation") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

nat_mat_plot <- gg_image_plot(nat_mat, 
                              x_vals = range(nat_mat$log_size),
                              y_vals = range(nat_mat$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Annual Temperature") +
  theme(legend.position = "none")

nat_tco_plot <- gg_image_plot(nat_tco, 
                              x_vals = range(nat_mat$log_size),
                              y_vals = range(nat_mat$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Temperature Coldest Quarter") +
  theme(legend.position = "none",
        axis.title.y = element_blank())



inv_map_plot <- gg_image_plot(inv_map, 
                              x_vals = range(inv_map$log_size),
                              y_vals = range(inv_map$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "Non-Native", subtitle = "Annual Precipitation") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

inv_mat_plot <- gg_image_plot(inv_mat, 
                              x_vals = range(inv_mat$log_size),
                              y_vals = range(inv_mat$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Annual Temperature") +
  theme(legend.position = "bottom") + 
  scale_fill_gradient("Var(Flower Production) (Cube Root)",
                      low = "red",
                      high = "yellow") +
  guides(fill = guide_colorbar(title.position = "top"))

inv_tco_plot <- gg_image_plot(inv_tco, 
                              x_vals = range(inv_mat$log_size),
                              y_vals = range(inv_mat$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_size,
                              dummy_clim) +
  ggtitle(label = "", subtitle = "Temperature Coldest Quarter") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

leg <- get_gg_legend(inv_mat_plot)

pdf(file = "repro_analysis/Manuscript/Figures/var_flower_n_clim.pdf",
    height = 10,
    width = 10)

  grid.arrange(
    nat_map_plot, nat_mat_plot, nat_tco_plot,
    inv_map_plot, inv_mat_plot + theme(legend.position = "none"), inv_tco_plot,
    leg,
    layout_matrix = matrix(
      c(1:6, NA_real_, 7, NA_real_),
      nrow = 3, ncol = 3, byrow = TRUE
    ),
    heights = c(3, 3, 1)
  )

dev.off()
