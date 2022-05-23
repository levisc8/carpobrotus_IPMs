# Preliminary data exploration

grow_data  <- readRDS("ipms/Data/growth_data.rds")
all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")
rm_ramets  <- read.csv('ipms/Data/t_2_omissions.csv',
                       stringsAsFactors = FALSE)

scat_theme <- theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  )

for(i in unique(all_ramets$site)) {
  
  temp_rm   <- filter(rm_ramets, population == i)
  temp_dat  <- filter(all_ramets, site == i)
  temp_gro  <- filter(grow_data, 
                      site == i & alive == 1 & !id %in% temp_rm$ramet)
  temp_sur  <- filter(all_ramets, 
                      site == i & !is.na(alive) & !id %in% temp_rm$ramet)
  
  temp_recr <- filter(temp_dat, id >= 8000)
  check_new <- filter(temp_dat, id < 8000 & is.na(size))
  
  s_plot <- ggplot(temp_sur, aes(x = log_size, y = alive)) +
    geom_jitter(width = 0, height = 0.05, size = 1.5) +
    scat_theme +
    ggtitle("Survival")
  
  g_plot <- ggplot(temp_gro, aes(x = log_size, y = log_size_next)) +
    geom_point(size = 1.5) +
    geom_abline(slope = 1, intercept = 0,
                color = "grey70",
                size = 1,
                linetype = "dashed") +
    # geom_smooth(formula = y ~ s(x, bs = "cs"),
    #             method = "gam",
    #             color = "red",
    #             linetype = "dotted",
    #             size = 1.5,
    #             se = FALSE) + 
    stat_smooth(formula = y ~ s(x, bs = "cs", k = 10),
                method = "gam",
                color = "red",
                size = 1.5,
                se = TRUE,
                method.args = list(family = gaussian())) +
    scat_theme +
    ggtitle("growth")
  
  pr_plot <- ggplot(temp_dat, aes(x = log_size, y = repro)) +
    geom_jitter(width = 0, height = 0.05, size = 1.5) +
    scat_theme +
    ggtitle("Pr(flowering)")
  
  fn_plot <- ggplot(temp_dat, aes(x = log_size, y = flower_n)) +
    geom_point(size = 1.5) +
    scat_theme +
    ggtitle("Flower #")
  
  rc_plot <- ggplot(temp_recr, aes(x = log_size_next)) +
    geom_histogram() +
    geom_vline(aes(xintercept = mean(log_size_next)),
               color = "red", 
               linetype = "dashed",
               size = 2) +
    scat_theme +
    ggtitle("Newly Observed Plants")
  
  dummy <- ggplot(check_new, aes(x = log_size_next)) +
    geom_histogram() +
    geom_vline(aes(xintercept = mean(log_size_next)),
               color = "red", 
               linetype = "dashed",
               size = 2) +
    scat_theme +
    ggtitle("Big New Plants")
  
  pdf(glue("ipms/Figures/site/prelim_plots_{i}.pdf"))
  
    grid.arrange(
      s_plot, g_plot,
      pr_plot, fn_plot,
      rc_plot, dummy,
      layout_matrix = matrix(c(1:6), 
                             nrow = 3,
                             ncol = 2,
                             byrow = TRUE)
    )
    
  dev.off()
  
}


# 
# clim_dat <- all_ramets %>%
#   mutate(incr = log_size_next - size) %>%
#   select(site:p_seas_rec, alive, incr, repro, flower_n) %>%
#   pivot_longer(cols      = Mean_Annual_Temp_hist:p_seas_rec, 
#                names_to  = "clim_var", 
#                values_to = "clim_val") %>%
#   pivot_longer(cols      = alive:flower_n,
#                names_to  = "demo_var",
#                values_to = "demo_val") 
# clim_dat$clim_val <- unlist(clim_dat$clim_val)
# 
# 
# for(i in unique(clim_dat$clim_var)) {
#   
#   surv_dat <- filter(clim_dat, clim_var == i & demo_var == "alive")
#   grow_dat <- filter(clim_dat, clim_var == i & demo_var == "incr")
#   repr_dat <- filter(clim_dat, clim_var == i & demo_var == "repro")
#   flow_dat <- filter(clim_dat, clim_var == i & demo_var == "flower_n")
#   
#   s <- ggplot(surv_dat, aes(x = clim_val, y = demo_val)) + 
#     geom_jitter() +
#     geom_smooth(method = "glm",
#                 formula = y ~ x,
#                 method.args = list(family = binomial())) +
#   ggtitle(i, subtitle = "survival")
#   
#   
#   g <- ggplot(grow_dat, aes(x = clim_val, y = demo_val)) + 
#     geom_point() +
#     geom_smooth(method = "glm",
#                 formula = y ~ x,
#                 method.args = list(family = gaussian())) +
#     ggtitle(i, subtitle = "Growth Increment")
#   
#   
#   r <- ggplot(repr_dat, aes(x = clim_val, y = demo_val)) + 
#     geom_jitter() +
#     geom_smooth(method = "glm",
#                 formula = y ~ x,
#                 method.args = list(family = binomial())) +
#     ggtitle(i, subtitle = "Pr(flowering)")
#   
#   
#   f <- ggplot(flow_dat, aes(x = clim_val, y = demo_val)) + 
#     geom_point() +
#     geom_smooth(method = "glm",
#                 formula = y ~ x,
#                 method.args = list(family = poisson())) + 
#     ggtitle(i, subtitle = "Flower #")
#     
#   pdf(glue("ipms/Figures/clim/prelim_plots_{i}.pdf"))
#   
#     grid.arrange(s, g, r, f, nrow = 2, ncol = 2)
#   
#   dev.off()  
# }
