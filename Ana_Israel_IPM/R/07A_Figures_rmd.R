# Final plots

# First, get the outputs into a reasonable format

const_var_temp <- lapply(const_var_out,
                         function(x) {
                           temp <- sort(x[2:1001])
                           return(c(x[1], temp[25], temp[975]))
                         }) %>%
  as_tibble() %>%
  mutate(boot_obs = c('Observed', 'Lower_CI', 'Upper_CI'))

exp_var_temp <- lapply(exp_var_out,
                       function(x) {
                         temp <- sort(x[2:1001])
                         return(c(x[1], temp[25], temp[975]))
                       }) %>%
  as_tibble() %>%
  mutate(boot_obs = c('Observed', 'Lower_CI', 'Upper_CI'))

for_plot_const <- const_var_temp %>%
  select(lambda, sd_g, g_int, g_slope, p_elas, f_elas, boot_obs) %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Constant_Variance')  %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)


for_plot_exp <- exp_var_temp %>%
  select(lambda, g_sigma_par, g_int, g_slope, p_elas, f_elas, boot_obs) %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Exponential_Variance')  %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)


pooled <- const_var_temp %>%
  select(-c(lambda, sd_g, g_int, g_slope, p_elas, f_elas)) %>% 
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Shared') %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)

for_plot <- rbind(for_plot_const, for_plot_exp, pooled)

# OK, data are now together. Next, we substitute in some fancy expressions
# so that figure titles look really nice
for_plot$Model <- gsub("_", " ", for_plot$Model)



for_plot$vital_rate[for_plot$vital_rate == 'lambda'] <- "paste(lambda)"

for_plot$vital_rate[for_plot$vital_rate == 'f_d_mu'] <- "paste(mu[f[d]])"
for_plot$vital_rate[for_plot$vital_rate == 'f_d_sd'] <- "paste(sigma[f[d]])"

for_plot$vital_rate[for_plot$vital_rate == 'p_r_slope'] <- "paste(italic(p[r](z)), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 'p_r_int'] <- "paste(italic(p[r](z)), ' Intercept')"


for_plot$vital_rate[for_plot$vital_rate == 'f_r_slope'] <- "paste(italic(f[r](z)), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 'f_r_int'] <- "paste(italic(f[r](z)), ' Intercept')"

for_plot$vital_rate[for_plot$vital_rate == 'f_s_slope'] <- "paste(italic(f[s](z)), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 'f_s_int'] <- "paste(italic(f[s](z)), ' Intercept')"

for_plot$vital_rate[for_plot$vital_rate == 'g_slope'] <- "paste(italic(mu[g(z^{','},~z)]), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 'g_int'] <- "paste(italic(mu[g(z^{','},~z)]), ' Intercept')"
for_plot$vital_rate[for_plot$vital_rate == 'g_sigma_par'] <- "paste(sigma[g(z^{','},~z)], ' Exponent')"
for_plot$vital_rate[for_plot$vital_rate == 'sd_g'] <- "paste(sigma[g(z^{','},~z)], ' Constant')"

for_plot$vital_rate[for_plot$vital_rate == 's_slope'] <- "paste(italic(s(z)), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 's_int'] <- "paste(italic(s(z)), ' Intercept')"

for_plot$vital_rate[for_plot$vital_rate == 'f_r'] <- "paste(italic(f[r]))"



elas_plot <- filter(for_plot, vital_rate %in% c('p_elas', 'f_elas'))
for_plot  <- filter(for_plot, !vital_rate %in% c('p_elas', 'f_elas'))

elas_plot$vital_rate[elas_plot$vital_rate == 'p_elas'] <- "paste('Elasticity to Survival/Growth')"
elas_plot$vital_rate[elas_plot$vital_rate == 'f_elas'] <- "paste('Elasticity to Fecundity')"

for_plot$Facet <- relevel(as.factor(for_plot$vital_rate), ref = 'paste(lambda)')
for_plot$Model <- factor(for_plot$Model, 
                         levels = c('Constant Variance',
                                    "Shared",
                                    "Exponential Variance"),
                         ordered = TRUE)

elas_plot$Model <- factor(elas_plot$Model, 
                         levels = c('Constant Variance',
                                    "Shared",
                                    "Exponential Variance"),
                         ordered = TRUE)

for_plot$yints <- NA
for_plot$yints[for_plot$vital_rate == 'paste(lambda)'] <- 1

# Set up default themes for contour plots and line plots

theme_contour <- theme(
  panel.background = element_blank(),
  axis.text        = element_text(size   = 14),
  axis.title.x     = element_text(size   = 16,
                                  margin = margin(
                                    t = 20,
                                    r = 0, 
                                    l = 0, 
                                    b = 2
                                  )
  ),
  axis.title.y     = element_text(size   = 16,
                                  margin = margin(
                                    t = 0,
                                    r = 20,
                                    l = 2,
                                    b = 0
                                  )
  )
)

theme_linerange <- theme_bw() + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_blank(), # Remove x-axis text + title
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 14), # make y-axis text + title bigger
    axis.title.y      = element_text(size = 16),
    strip.text        = element_text(size = 20), # Increases plot label size
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 14)
  )

# Now, make the figure panel

vr_plot <- ggplot(for_plot) + 
  geom_point(
    aes(
      x     = Model, 
      y     = Observed,
      color = Model
    ),
    size = 2.5
  ) + 
  geom_linerange(
    aes(
      x     = Model,
      ymin  = Lower_CI,
      ymax  = Upper_CI,
      color = Model
    ),
    size = 1.3
  ) + 
  geom_hline(
    aes(
      yintercept = yints
    ),
    color       = 'grey50',
    alpha       = 0.7,
    linetype    = 'dashed',
    size        = 1.2,
    na.rm       = TRUE,
    show.legend = FALSE
  ) + 
  facet_wrap(
    ~ Facet,
    scales   = 'free_y',
    labeller = label_parsed) +
  scale_color_manual(
    breaks = c('Constant Variance',
               'Shared',
               'Exponential Variance'),
    values = c('grey70', 'grey50', 'black')
  ) +
  theme_linerange +
  theme(legend.position   = c(0.6125, 0.125))


# Full kernel elasticity/sensitivity contour plots

k_sens_cv_plot <- ggplot(k_sens_const) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = 'ln(Surface Area, T + 1)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour


k_elas_cv_plot <- ggplot(k_elas_const) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour


k_sens_ev_plot <- ggplot(k_sens_exp) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = 'ln(Surface Area, T)',
                     limits = c(-5, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = 'ln(Surface Area, T + 1)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour


k_elas_ev_plot <- ggplot(k_elas_exp) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = 'ln(Surface Area, T)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour

# Sub-kernel elasticity plots ------------

P_elas_cv_plot <- ggplot(P_elas_cv) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = 'ln(Surface Area, T + 1)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour

P_elas_ev_plot <- ggplot(P_elas_ev) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = 'ln(Surface Area, T)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = 'ln(Surface Area, T + 1)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour

F_elas_cv_plot <- ggplot(F_elas_cv) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour

F_elas_ev_plot <- ggplot(F_elas_ev) +
  geom_tile(aes(x = x, y = y, fill = value),
            show.legend = FALSE) +
  geom_contour(aes(x = x, y = y, z = value),
               color = 'black',
               size = 1.1) +
  scale_fill_gradient(low = 'red', high = 'yellow') +
  scale_x_continuous(name = 'ln(Surface Area, T)',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  scale_y_continuous(name = '',
                     limits = c(-4.9, 3.59),
                     breaks = seq(min(d1), max(d1), length.out = 5),
                     labels = round(seq(min(d1), max(d1), length.out = 5),
                                    digits = 3)) + 
  theme_contour



sub_kern_elas_plot <- ggplot(elas_plot) + 
  geom_point(
    aes(
      x     = Model, 
      y     = Observed,
      color = Model
    ),
    size = 2.5
  ) + 
  geom_linerange(
    aes(
      x     = Model,
      ymin  = Lower_CI,
      ymax  = Upper_CI,
      color = Model
    ),
    size = 1.3
  ) + 
  facet_wrap(
    ~ vital_rate,
    labeller = label_parsed) +
  scale_color_manual(
    breaks = c('Constant Variance',
               'Shared',
               'Exponential Variance'),
    values = c('grey70', 'grey50', 'black')
  ) +
  theme_linerange




