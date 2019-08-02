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
  select(lambda, sd_g, g_int, g_slope, boot_obs) %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Constant_Variance')  %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)


for_plot_exp <- exp_var_temp %>%
  select(lambda, g_sigma_par, g_int, g_slope, boot_obs) %>%
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Exponential_Variance')  %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)


pooled <- const_var_temp %>%
  select(-c(lambda, sd_g, g_int, g_slope)) %>% 
  gather(key = 'vital_rate', value = 'value', -boot_obs) %>%
  mutate(Model = 'Shared') %>%
  spread(key = 'boot_obs', value = 'value', drop = FALSE)

for_plot <- rbind(for_plot_const, for_plot_exp, pooled)

# OK, data are now together. Next, we substitute in some fancy expressions
# so that figure titles look really nice
for_plot$Model <- gsub("_", " ", for_plot$Model)


for_plot$vital_rate[for_plot$vital == 'sd_g'] <- "paste(italic(sigma[g(z^{','},~z))]))"

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
for_plot$vital_rate[for_plot$vital_rate == 'sd_g'] <- "paste(sigma[g(z^{','},~z)],
' Constant')"


for_plot$vital_rate[for_plot$vital_rate == 's_slope'] <- "paste(italic(s(z)), ' Slope')"
for_plot$vital_rate[for_plot$vital_rate == 's_int'] <- "paste(italic(s(z)), ' Intercept')"

for_plot$vital_rate[for_plot$vital_rate == 'f_r'] <- "paste(italic(f[r]))"

for_plot$Facet <- relevel(as.factor(for_plot$vital_rate), ref = 'paste(lambda)')
for_plot$Model <- factor(for_plot$Model, 
                         levels = c('Constant Variance',
                                    "Shared",
                                    "Exponential Variance"),
                         ordered = TRUE)
for_plot$yints <- NA
for_plot$yints[for_plot$vital_rate == 'paste(lambda)'] <- 1

# Now, make the figure panel

ggplot(for_plot) + 
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
    labeller = label_parsed
  ) +
  theme_bw() + 
  scale_color_manual(
    breaks = c('Constant Variance',
               'Shared',
               'Exponential Variance'),
    values = c('grey70', 'grey50', 'black')
  ) + 
  theme( 
    # Extras to theme_bw()
    axis.text.x       = element_blank(), # Remove x-axis text + title
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 14), # make y-axis text + title bigger
    axis.title.y      = element_text(size = 16),
    strip.text        = element_text(size = 20), # Increases plot label size
    legend.position   = c(0.6125, 0.125),
    legend.background = element_rect(fill = NA,  # Box for the legend 
                                     color = 'black'),
    legend.text       = element_text(size = 12),
    legend.title      = element_text(size = 14)
  )


ggsave(filename = 'VR_Panel.png',
       path = 'Ana_Israel_IPM/Figures/',
       device = 'png',
       height = 10,
       width = 14,
       units = 'in')
