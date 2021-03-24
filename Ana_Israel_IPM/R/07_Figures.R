# Final plots

exp_var_temp <- lapply(exp_var_out,
                       function(x) {
                         temp <- sort(x[2:1001])
                         return(c(x[1], temp[25], temp[975]))
                       }) %>%
  as_tibble() %>%
  mutate(boot_obs = c('Observed', 'Lower_CI', 'Upper_CI'))

surv_pred <- data.frame(log_size = xx, 
                        pred = predict(surv_mod_lin, 
                                       data.frame(log_size = xx),
                                       type = 'response'))

grow_pred <- data.frame(log_size = xx,
                        pred = predict(grow_exp_var,
                                       data.frame(log_size = xx),
                                       type = 'response'))

pr_pred <- data.frame(log_size = xx,
                        pred = predict(p_r_mod_lin,
                                       data.frame(log_size = xx),
                                       type = 'response'))
fs_pred <- data.frame(log_size = xx,
                      pred = predict(f_s_mod,
                                     data.frame(log_size = xx),
                                     type = 'response'))

yy <- seq(min(recruit_t_2$log_size_next, na.rm = TRUE) - 10,
          max(recruit_t_2$log_size_next, na.rm = TRUE) + 10,
          length.out = 400)

recr_pred <- data.frame(log_size_next = yy,
                        density = dnorm(yy, 
                                        f_d_mu,
                                        f_d_sd))


all_vr <- exp_var_temp %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)

lam_pred <- filter(all_vr, vital_rate == 'lambda')

elas_plot <- filter(all_vr, vital_rate %in% c('p_elas', 'f_elas'))


# Set up default themes for contour plots and line plots

theme_contour <- theme(
  panel.background = element_blank(),
  axis.text        = element_text(size   = 8),
  axis.title.x     = element_text(size   = 10,
                                  margin = margin(
                                    t = 5,
                                    r = 0, 
                                    l = 0, 
                                    b = 1
                                  )
  ),
  axis.title.y     = element_text(size   = 8,
                                  margin = margin(
                                    t = 0,
                                    r = 5,
                                    l = 1,
                                    b = 0
                                  )
  )
)

theme_linerange <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_blank(), # Remove x-axis text + title
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 7), # make y-axis text + title bigger
    axis.title.y      = element_text(size = 9,
                                     margin = margin(
                                       t = 0,
                                       l = 2,
                                       r = 5,
                                       b = 0
                                     )),
    strip.text        = element_text(size = 10), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7)
  )

theme_vr <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_text(size = 7), 
    axis.text.y       = element_text(size = 7), # make y-axis text + title bigger
    axis.title.x      = element_text(size   = 8,
                                     margin = margin(
                                       t = 5,
                                       r = 0, 
                                       l = 0, 
                                       b = 7
                                     )
    ),
    axis.title.y     = element_text(size   = 8,
                                    margin = margin(
                                      t = 3,
                                      r = 10,
                                      l = 1,
                                      b = 0
                                    )
    ),
    strip.text        = element_text(size = 8), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7)
  )

# Now, make the figure panel

grow_plot <- ggplot(all_data, aes(x = log_size,
                                  y = log_size_next)) +
  geom_point(color = 'black', 
             size = 1.25) + 
  geom_line(data = grow_pred,
            aes(x = log_size,
                y = pred),
            linetype = 'dashed',
            size = 1.25,
            color = 'black',
            show.legend = FALSE) + 
  geom_abline(intercept = 0,
              slope = 1,
              color = 'grey70',
              show.legend = FALSE,
              size = 1.25) +
  theme_vr +
  scale_x_continuous('ln(Surface Area, t)',
                     limits = c(-6.1, 3.5)) +
  scale_y_continuous('ln(Surface Area, t + 1)',
                     limits = c(-6.1, 3.5))

surv_plot <- ggplot(all_data, aes(x = log_size,
                                  y = survival)) +
  geom_jitter(color = 'black',
              size = 1.75,
              width = 0,
              height = 0.05) + 
  geom_line(data = surv_pred,
            aes(x = log_size,
                y = pred),
            linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'black') +
  theme_vr + 
  scale_x_continuous('ln(Surface Area, t)', 
                     limits = c(-5.5, 3.5)) +
  scale_y_continuous('Survival (t + 1)',
                     limits = c(-0.1, 1.1),
                     breaks = c(0, 1))

pr_plot <- ggplot(all_data, aes(x = log_size, 
                                y = repro)) +
  geom_jitter(color = 'black',
              size = 1.75,
              width = 0,
              height = 0.05) + 
  geom_line(data = pr_pred,
            aes(x = log_size,
                y = pred),
            linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'black') +
  theme_vr + 
  scale_x_continuous('ln(Surface Area, t)', 
                     limits = c(-5.5, 3.5)) +
  scale_y_continuous('Pr(Reproductive, t)',
                     limits = c(-0.1, 1.1),
                     breaks = c(0, 1))

fs_plot <- ggplot(all_data, aes(x = log_size, 
                                y = flower_n)) +
  geom_point(color = 'black',
              size = 1.75) + 
  geom_line(data = fs_pred,
            aes(x = log_size,
                y = pred),
            linetype = 'dashed',
            size = 1.25,
            show.legend = FALSE,
            color = 'black') +
  theme_vr + 
  scale_x_continuous('ln(Surface Area, t)', 
                     limits = c(-5.5, 3.5)) +
  scale_y_continuous('# of Flowers',
                     limits = c(0, 70),
                     breaks = seq(0, 70, by = 15))


recr_plot <- ggplot(recruit_t_2, aes(x = log_size_next)) +
  geom_histogram(aes(y = ..density..),
                 bins = 10,
                 fill = NA,
                 color = 'black') +
  geom_line(data = recr_pred,
                aes(x = log_size_next,
                    y = density)) + 
  theme_vr + 
  scale_x_continuous('ln(Surface Area, t + 1)',
                     limits = c(-7, 2)) + 
  scale_y_continuous('Probability Density', 
                     limits = c(0, 0.5)) + 
  theme(panel.grid = element_blank())

lam_plot <- ggplot(lam_pred, 
                   aes(x = vital_rate,
                       y = Observed)) + 
  geom_point(color = 'black',
             size = 4) + 
  geom_linerange(aes(ymin = Lower_CI,
                     ymax = Upper_CI),
                 color = 'black',
                 size = 1.75) + 
  theme_linerange + 
  scale_y_continuous(parse(text = 'lambda'),
                     limits = c(0.93, 1.02)) +
  geom_hline(yintercept = 1,
             color = 'grey80',
             linetype = 'dashed',
             size = 2)

tiff(filename = 'Ana_Israel_IPM/Figures/Figure_3.tiff',
    height = 5,
    width = 6,
    units = 'in',
    res = 300)

  grid.arrange(grow_plot, surv_plot,
               pr_plot, fs_plot,
               recr_plot, lam_plot,
               nrow = 3, ncol = 2)
  
dev.off()


# Full kernel elasticity/sensitivity contour plots

all_sens <- k_sens_exp

k_sens_plot <- ggplot(all_sens) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 0.4) + 
  scale_fill_gradient("Sensitivity",
                      low = 'red',
                      high = 'yellow') +
  scale_x_continuous(name = '',
                     limits = c(-7.3, 4),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = 'ln(Surface Area, t + 1)',
                     limits = c(-7.3, 4),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour + 
  theme(legend.position = 'left',
        legend.key = element_rect(size = unit(4, 'in')),
        legend.title = element_text(size = 14),
        strip.text = element_blank(),
        strip.background = element_rect(fill = NULL,
                                        color = NULL))


all_elas <- k_elas_exp

k_elas_plot <- ggplot(all_elas) +
  geom_tile(aes(x = y, y = x, fill = value)) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 0.4)  + 
  scale_fill_gradient("Elasticity",
                      low = 'red',
                      high = 'yellow') +
  scale_x_continuous(name = 'ln(Surface Area, t)',
                     limits = c(-7.3, 4),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = '',
                     limits = c(-7.3, 4),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour + 
  theme(legend.position = 'left',
        legend.key = element_rect(size = unit(4, 'in')),
        legend.title = element_text(size = 14),
        strip.text = element_blank(),
        strip.background = element_rect(fill = NULL,
                                        color = NULL))

sub_kern_elas_plot <- ggplot(elas_plot) + 
  geom_point(
    aes(
      x     = vital_rate, 
      y     = Observed,
      color = vital_rate
    ),
    size = 2.5,
    show.legend = FALSE
  ) + 
  geom_linerange(
    aes(
      x     = vital_rate,
      ymin  = Lower_CI,
      ymax  = Upper_CI,
      color = vital_rate
    ),
    size = 1.3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    breaks = c('f_elas',
               'p_elas'),
    values = c('grey70', 'black')
  ) +
  scale_x_discrete("",
                   breaks = c('f_elas', 'p_elas'),
                   labels = c('Sexual Reproduction',
                              'Survival/Growth')) +
  scale_y_continuous('Elasticity') + 
  theme_vr

tiff(filename = 'Ana_Israel_IPM/Figures/Figure_4.tiff',
    height = 6,
    width = 8,
    units = 'in',
    res = 300)

  grid.arrange(k_sens_plot,
               k_elas_plot,
               sub_kern_elas_plot,
               layout_matrix = matrix(c(1, 3,
                                        2, 3),
                                      nrow = 2, 
                                      byrow = TRUE))

dev.off()


# population summary plots 

theme_sum <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
   axis.text = element_text(size = 12),
   axis.title = element_text(size = 14),
   strip.text        = element_text(size = 12), # Increases plot label size
   strip.background  = element_rect(fill = NA,
                                    color = "black"),
   legend.background = element_rect(fill = NA,  # Box for the legend 
                                    color = 'black'),
   legend.text       = element_text(size = 6),
   legend.title      = element_text(size = 7)
  )


f_n_dat <- all_data %>% 
  select(id, flower_n, flower_n_next) %>% 
  gather(-id, key = "time", value = "value")

size_dat <- all_data %>% 
  select(id, log_size, log_size_next) %>% 
  gather(-id, key = "time", value = "value")

tot_size <- all_data %>% 
  select(id, size, size_next) %>%
  gather(-id, key = "time", value = "value") %>%
  group_by(time) %>%
  summarise(tot_size =  paste("Total area: ", 
                              round(sum(value, na.rm = TRUE), 3),
                              sep = ""))


f_n_dat$time <- ifelse(f_n_dat$time == "flower_n", "t", "t + 1")
size_dat$time <- ifelse(size_dat$time == "log_size", "t", "t + 1")
tot_size$time <- ifelse(tot_size$time == "size", "t", "t + 1")


f_n_hist <- ggplot(f_n_dat, aes(x = value)) + 
  geom_histogram(fill = NA, 
                 color = "black", 
                 size = 1.25) + 
  facet_wrap(~ time,
             scales = "free_x") + 
  labs(x = "Flower #", y = "Count") + 
  theme_sum 

size_2 <- size_dat %>%
  group_by(time) %>% 
  summarise(mean_value = mean(value, na.rm = TRUE),
            min_val = paste("Minimum size: ",
                            round(min(value, na.rm = TRUE), 3), 
                            sep = ""),
            max_val = paste("Maximum size: ",
                            round(max(value, na.rm = TRUE), 3),
                            sep = ""),
            mean_txt = paste("Mean size: ", round(mean_value, 3), sep = ""))

size_hist <- ggplot(size_dat, aes(x = value)) + 
  geom_histogram(fill = NA, 
                 color = "black", 
                 size = 1.25,
                 binwidth = 0.25) +
  geom_vline(data = size_2,
             aes(xintercept = mean_value), 
             size = 1.25,
             color = "red",
             linetype = "longdash") + 
  facet_wrap(~time) + 
  labs(x = "Ln(Surface Area)", y = "Count") + 
  theme_sum + 
  # theme(axis.title = element_text(size = 16))  +
  geom_text(aes(x = -6, y = 23, label = tot_size), data = tot_size, hjust = 0) +
  geom_text(aes(x = -6, y = 22, label = mean_txt), data = size_2, hjust = 0) +
  geom_text(aes(x = -6, y = 21, label = min_val), data = size_2, hjust = 0) +
  geom_text(aes(x = -6, y = 20, label = max_val), data = size_2, hjust = 0)

tiff(filename = 'Ana_Israel_IPM/Figures/Figure_2.tiff',
     height = 9,
     width = 9,
     units = 'in',
     res = 300)

  grid.arrange(f_n_hist,
               size_hist,
               nrow = 2, ncol = 1)

dev.off()

png(filename = 'Ana_Israel_IPM/Figures/Figure_2.png',
     height = 9,
     width = 9,
     units = 'in',
     res = 300)

grid.arrange(f_n_hist,
             size_hist,
             nrow = 2, ncol = 1)

dev.off()
