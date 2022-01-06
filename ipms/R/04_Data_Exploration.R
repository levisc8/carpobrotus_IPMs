# Preliminary data exploration

library(ggplot2)
library(gridExtra)

grow_data <- readRDS("ipms/Data/growth_data.rds")
all_ramets <- readRDS("ipms/Data/all_ramets_di.rds")

scat_theme <- theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  )

for(i in unique(all_ramets$population)) {
  
  temp_dat  <- filter(all_ramets, population == i)
  temp_gro  <- filter(grow_data, population == i & alive == 1)
  
  temp_recr <- filter(temp_dat, id >= 8000)
  check_new <- filter(temp_dat, id < 8000 & is.na(size))
  
  s_plot <- ggplot(temp_dat, aes(x = size, y = alive)) +
    geom_jitter(width = 0, height = 0.05, size = 1.5) +
    scat_theme +
    ggtitle("Survival")
  
  g_plot <- ggplot(temp_gro, aes(x = size, y = size_next)) +
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
                method.args = list(family = Gamma(link = "identity"))) +
    scat_theme +
    ggtitle("growth")
  
  pr_plot <- ggplot(temp_dat, aes(x = size, y = repro)) +
    geom_jitter(width = 0, height = 0.05, size = 1.5) +
    scat_theme +
    ggtitle("Pr(flowering)")
  
  fn_plot <- ggplot(temp_dat, aes(x = size, y = flower_n)) +
    geom_point(size = 1.5) +
    scat_theme +
    ggtitle("Flower #")
  
  rc_plot <- ggplot(temp_recr, aes(x = size_next)) +
    geom_histogram() +
    geom_vline(aes(xintercept = mean(size_next)),
               color = "red", 
               linetype = "dashed",
               size = 2) +
    scat_theme +
    ggtitle("Newly Observed Plants")
  
  dummy <- ggplot(check_new, aes(x = size_next)) +
    geom_histogram() +
    geom_vline(aes(xintercept = mean(size_next)),
               color = "red", 
               linetype = "dashed",
               size = 2) +
    scat_theme +
    ggtitle("Big New Plants")
  
  pdf(glue("ipms/Figures/prelim_plots_{i}.pdf"))
  
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

dev.off()