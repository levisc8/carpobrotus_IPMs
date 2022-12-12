
all_ramets <- readRDS("ipms/Data/all_ramets_di.rds") %>%
  mutate(
    native = case_when(
      site %in% c("Melkboss",      "Vogelgat", 
                  "St_Francis",    "Struisbaai",
                  "Springfontein", "Rooisand")   ~ 1,
      TRUE ~ 0
    )
  )

surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig.rds")

use_surv_mod <- surv_mods[[1]][["times_sw3_seas"]]
use_grow_mod <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
use_repr_mod <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")
use_flow_mod <- flow_mods[[1]][["times_sw2_ann"]]

surv_pred <- fitted(use_surv_mod)
grow_pred <- fitted(use_grow_mod)
repr_pred <- fitted(use_repr_mod)
flow_pred <- fitted(use_flow_mod)

surv_data <- cbind(use_surv_mod$data, 
                   surv_pred = surv_pred[ ,1])
grow_data <- cbind(use_grow_mod$data,
                   grow_pred = grow_pred[ ,1])
repr_data <- cbind(use_repr_mod$data, 
                   repr_pred = repr_pred[ ,1])
flow_data <- cbind(use_flow_mod$data,
                   flow_pred = flow_pred[ ,1])

# loop over populations so that each one gets its own pdf

for(i in seq_along(unique(all_ramets$site))) {
  
  pop <- unique(all_ramets$site)[i]

  surv_temp <- filter(surv_data, site == pop)
  grow_temp <- filter(grow_data, site == pop)
  repr_temp <- filter(repr_data, site == pop)
  flow_temp <- filter(flow_data, site == pop)
  
  
  surv_plt <- ggplot(surv_temp) +
    geom_line(aes(x = log_size, y = surv_pred), 
              color = "red", size = 1.5) +
    geom_jitter(aes(x = log_size, y = alive), alpha = 0.5,
                height = 0.05, width = 0) +
    theme_bw() +
    ylab("Survival (t + 1)") + 
    xlab("Size (t)") +
    ggtitle(pop)
  
  grow_plt <- ggplot(grow_temp) +
    geom_line(aes(x = log_size, y = grow_pred),
              color = "red", size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "grey50") +
    geom_point(aes(x = log_size, y = log_size_next), alpha = 0.5) +
    theme_bw() +
    ylab("Size (t + 1)") + 
    xlab("Size (t)")
  
  repr_plt <- ggplot(repr_temp) +
    geom_line(aes(x = log_size, y = repr_pred),
              color = "red", size = 1.5) +
    geom_jitter(aes(x = log_size, y = repro),
              alpha = 0.5, height = 0.05, width = 0) +
    theme_bw() +
    ylab("Pr(Reproductive) (t)") + 
    xlab("Size (t)") +
    ggtitle(pop)
  
  flow_plt <- ggplot(flow_temp) +
    geom_line(aes(x = log_size, y = flow_pred),
              color = "red", size = 1.5) +
    geom_point(aes(x = log_size, y = flower_n),
               alpha = 0.5) +
    theme_bw() +
    ylab("# of Flowers (t)") + 
    xlab("Size (t)")
  
  
  pdf(glue('ipms/Figures/vr_models/gam/{pop}_p_kern_preliminary_plots.pdf'),
      width = 8,
      height = 8)
  
    grid.arrange(surv_plt, grow_plt, nrow = 2, ncol = 1)
    
  dev.off()
  
  pdf(glue('ipms/Figures/vr_models/gam/{pop}_f_kern_preliminary_plots.pdf'),
      width = 8,
      height = 8)
  
    grid.arrange(repr_plt, flow_plt, nrow = 2, ncol = 1)
  
  dev.off()
  
  
}


use_sites <- c("Havatselet", "Foxton", "Rooisand", "Colares", "Rarangi")

use_surv <- filter(surv_data, site %in% use_sites)
use_grow <- filter(grow_data, site %in% use_sites)
use_repr <- filter(repr_data, site %in% use_sites)
use_flow <- filter(flow_data, site %in% use_sites)

surv_plt <- ggplot(use_surv) +
  geom_line(aes(x = log_size, y = surv_pred, color = site),
            size = 1.5, alpha = 0.7) +
  theme_bw() +
  ylab("Survival (t + 1") +
  xlab("Size (t)")

grow_plt <- ggplot(use_grow) +
  geom_line(aes(x = log_size, y = grow_pred, color = site),
            size = 1.5, alpha = 0.7) +
  theme_bw() +
  ylab("Size (t + 1") +
  xlab("Size (t)")

repr_plt <- ggplot(use_repr) +
  geom_line(aes(x = log_size, y = repr_pred, color = site),
            size = 1.5, alpha = 0.7) +
  theme_bw() +
  ylab("Pr(Reproductive) (t + 1") +
  xlab("Size (t)")

flow_plt <- ggplot(use_flow) +
  geom_line(aes(x = log_size, y = flow_pred, color = site),
            size = 1.5, alpha = 0.7) +
  theme_bw() +
  ylab("Pr(Reproductive) (t + 1") +
  xlab("Size (t)")

pdf(glue('ipms/Figures/vr_models/gam/all_pops_f_kern_preliminary_plots.pdf'),
    width = 8,
    height = 8)

  grid.arrange(repr_plt, flow_plt, nrow = 2, ncol = 1)

dev.off()

pdf(glue('ipms/Figures/vr_models/gam/all_pops_p_kern_preliminary_plots.pdf'),
    width = 8,
    height = 8)

  grid.arrange(surv_plt, grow_plt, nrow = 2, ncol = 1)

dev.off()


#  FInally, create pp_check plots for everything for appendix of publication.

# Survival

pdf("ipms/Figures/vr_models/best_mods/survival.pdf")

pp_check(use_surv_mod,
         type   = 'bars_grouped',
         group  = 'site', 
         freq   = FALSE,
         ndraws = 100L) 

dev.off()

# 2_1_test is best for growth is best overall. Create pp_check for it!

pdf("ipms/Figures/vr_models/best_mods/growth.pdf")

pp_check(use_grow_mod,
         type   = 'scatter_avg_grouped',
         group  = 'site',
         ndraws = 100L) + 
  geom_abline(slope = 1, intercept = 0)

plot(conditional_smooths(use_grow_mod,
                         surface = TRUE), 
     stype = 'raster',
     theme = theme_bw(),
     ask = FALSE)

dev.off()

# Reproduction

pdf("ipms/Figures/vr_models/best_mods/reproduction.pdf")

pp_check(use_repr_mod,
         type   = 'bars_grouped',
         group  = 'site',
         freq   = FALSE,
         ndraws = 100L) 

plot(conditional_smooths(use_repr_mod,
                         surface = TRUE), 
     stype = 'raster',
     theme = theme_bw(),
     ask = FALSE)

dev.off()

# Flower Number

pdf("ipms/Figures/vr_models/best_mods/flower_n.pdf")

pp_check(use_flow_mod,
         type   = 'scatter_avg_grouped',
         group  = 'site',
         ndraws = 100L) + 
  geom_abline(slope = 1, intercept = 0)

dev.off()
