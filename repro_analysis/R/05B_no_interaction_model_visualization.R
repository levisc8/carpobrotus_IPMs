# Load fitted model objects

vr_mod_list <- readRDS("repro_analysis/model_fits/vr_no_nat_mod_list_centered.rds") 

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
                             t_co_qu_rec = t_c_seq)

mat_pred_dat <- expand.grid(log_size    = size_seq,
                            mat_rec     = mat_seq)

map_pred_dat <- expand.grid(log_size    = size_seq,
                            map_rec     = map_seq)

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
  rbind(pred_repro_mat, pred_repro_t_co_qu)


# Temporary predict method for stanfit objects. Should only be used 
# for these specific models!

predict.stanfit <- function(mod, new_data) {
  
  # Add intercept column to create a design matrix for the new 
  # data
  
  new_data <- as.matrix(cbind(1, new_data))
  
  # Use the means from posterior distributions to create
  # a coefficient vector
  
  pars     <- rstan::extract(mod, 
                             pars = c("b_Intercept", "b")) %>%
    do.call(what = "cbind", args = .) %>%
    apply(2, mean)
  
  # Prediction vector new data, and wrap it up in a more friend interface 
  # so it'll stack with the pr(repro) predictions
  
  out <- data.frame(Estimate = exp(new_data %*% pars))
  
  return(out)
  
}

pred_flower_t_co_qu  <- predict(vr_mod_list$flower_mod_t_co,
                                new_data = t_co_pred_dat) %>%
  as_tibble() %>%
  cbind(t_co_pred_dat) %>%
  mutate(clim_var = "t_co_qu_rec",
         clim_val = t_co_qu_rec) %>%
  select(-t_co_qu_rec)

pred_flower_mat      <- predict(vr_mod_list$flower_mod_mat,
                                new_data = mat_pred_dat) %>%
  as_tibble() %>%
  cbind(mat_pred_dat) %>%
  mutate(clim_var = "mat_rec",
         clim_val = mat_rec) %>%
  select(-mat_rec)

all_flower_pred      <- predict(vr_mod_list$flower_mod_map,
                                new_data = map_pred_dat) %>%
  as_tibble() %>%
  cbind(map_pred_dat)%>%
  mutate(clim_var = "map_rec",
         clim_val = map_rec) %>%
  select(-map_rec)  %>%
  rbind(pred_flower_mat, pred_flower_t_co_qu)
#
model_preds <- list(flower_preds = all_flower_pred,
                    repro_preds = all_repro_pred)

saveRDS(model_preds,
        file = "repro_analysis/Data/demography/no_nat_model_predictions.rds")

model_preds <- readRDS("repro_analysis/Data/demography/no_nat_model_predictions.rds")

# For plotting, we just want to see the effects of climate, not size.
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
  pivot_longer(-c(id:id_pop),
               names_to = "clim_var",
               values_to = "clim_val") %>%
  filter(log_size > mean(log_size))

f_plt <- ggplot(for_plot, aes(y = flower_n, x = clim_val)) +
  geom_point() +
  geom_line(data = all_flower_pred,
            aes(y = Estimate,
                x = clim_val), 
            show.legend = FALSE,
            size = 1.5,
            linetype = "dashed") + 
  facet_wrap(~clim_var,
             scales = "free_x") +
  theme_bw() +
  ylim(c(0, max(for_plot$flower_n, na.rm = TRUE) + 100))


pdf("repro_analysis/Manuscript/Figures/flower_mod_no_nat_viz.pdf")

print(f_plt)

dev.off()

r_plt <- ggplot(for_plot, 
                aes(y = repro, x = clim_val)) +
  geom_jitter(height = 0.05,
              width  = 0) +
  geom_line(data = all_repro_pred,
            aes(y = Estimate,
                x = clim_val), 
            show.legend = FALSE,
            size = 1.5,
            linetype = "dashed") + 
  facet_wrap(~clim_var,
             scales = "free_x") +
  theme_bw() +
  ylim(c(-0.1, 1.1))

pdf("repro_analysis/Manuscript/Figures/repro_mod_no_nat_viz.pdf")

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
if(all.equal(model_preds$flower_preds[ , 2:4],
             model_preds$repro_preds[ , c(5:7)])) {
  
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
  select(-Estimate)

sum_plot_data$omega_z[sum_plot_data$omega_z < 1]     <- NA
sum_plot_data$var_omega_z[sum_plot_data$omega_z < 1] <- NA


dat_map <- filter(sum_plot_data, clim_var == "map_rec")

dat_mat <- filter(sum_plot_data, clim_var == "mat_rec")

dat_tco <- filter(sum_plot_data, clim_var == "t_co_qu_rec")

nat_map_plot <- gg_image_plot(dat_map, 
                              y_vals = range(dat_map$log_size),
                              x_vals = range(dat_map$clim_val),
                              omega_z  ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  ggtitle(label = "Annual Precipitation") +
  theme(legend.position = "none",
        axis.title.x = element_blank())

nat_mat_plot <- gg_image_plot(dat_mat, 
                              y_vals = range(dat_mat$log_size),
                              x_vals = range(dat_mat$clim_val),
                              omega_z ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  ggtitle(label = "Annual Temperature") +
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

nat_tco_plot <- gg_image_plot(dat_tco, 
                              y_vals = range(dat_mat$log_size),
                              x_vals = range(dat_mat$clim_val),
                              omega_z  ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  ggtitle(label = "Temperature Coldest Quarter") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

m_leg <- get_gg_legend(nat_mat_plot)


# Variance in omega_bar
var_map_plot <- gg_image_plot(dat_map, 
                              y_vals = range(dat_map$log_size),
                              x_vals = range(dat_map$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

var_mat_plot <- gg_image_plot(dat_mat, 
                              y_vals = range(dat_mat$log_size),
                              x_vals = range(dat_mat$clim_val),
                              var_omega_z ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  theme(legend.position = "bottom",
        axis.title.y = element_blank())

var_tco_plot <- gg_image_plot(dat_tco, 
                              y_vals = range(dat_mat$log_size),
                              x_vals = range(dat_mat$clim_val),
                              var_omega_z  ^ (1/3),
                              dummy_clim,
                              dummy_size) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank())


v_leg <- get_gg_legend(var_mat_plot)

pdf(file = "repro_analysis/Manuscript/Figures/flower_n_no_nat_clim.pdf",
    height = 10,
    width = 10)

grid.arrange(
  m_leg,
  nat_map_plot, nat_mat_plot + theme(legend.position = "none"), nat_tco_plot,
  var_map_plot, var_mat_plot + theme(legend.position = "none"), var_tco_plot,
  v_leg,
  layout_matrix = matrix(
    c(NA_real_, 1, NA_real_,
      2:7,
      NA_real_, 8, NA_real_),
    nrow = 4, ncol = 3, byrow = TRUE
  ),
  heights = c(1, 3, 3, 1)
)

dev.off()
