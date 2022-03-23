
all_ramets <- readRDS("ipms/Data/all_ramets_di.rds") %>%
  mutate(
    native = case_when(
      site %in% c("Melkboss",      "Vogelgat", 
                  "St_Francis",    "Struisbaai",
                  "Springfontein", "Rooisand")   ~ 1,
      TRUE ~ 0
    )
  )

clim_data <- all_ramets %>% 
  select(site:ann_pc_2_t, native) %>%
  group_by(site) %>%
  summarise_all(.funs = unique) %>%
  ungroup()

pred_data <- all_ramets %>%
  group_by(site) %>%
  summarise(
    med_size = median(log_size, na.rm = TRUE),
    avg_size = mean(log_size, na.rm = TRUE),
    l1_size  = quantile(log_size, 0.1, na.rm = TRUE),
    u9_size  = quantile(log_size, 0.9, na.rm = TRUE)
  ) %>%
  left_join(clim_data)

# Create plot for median, mean, .9, and .1 quantile sized plants for each
# climate and population setup
# 
surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
grow_mods <- readRDS("ipms/Model_Fits/ramet_growth_list_krig_gam.rds")
repr_mods <- readRDS("ipms/Model_Fits/ramet_repro_list_krig_gam.rds")
flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig_gam.rds")
# 


use_surv_mod <- surv_mods[[1]][["times_sw2_seas"]]
use_grow_mod <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
use_repr_mod <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")
use_flow_mod <- flow_mods[[1]][["times_sw2_seas"]]

clim_names <- names(clim_data) %>%
  .[!. %in% c("site", "native")]

surv_terms <- mod_clim_terms(use_surv_mod, clim_names)
grow_terms <- mod_clim_terms(use_grow_mod, clim_names)
repr_terms <- mod_clim_terms(use_repr_mod, clim_names)
flow_terms <- mod_clim_terms(use_flow_mod, clim_names)


# loop over populations so that each one gets its own pdf

for(i in seq_along(unique(all_ramets$site))) {
  
  pop <- unique(all_ramets$site)[i]
  
  
  ramet_temp <- filter(pred_data, site == pop)

  sizes     <- select(ramet_temp, site:u9_size, native) %>%
    pivot_longer(-c(site, native), names_to = "size_qtle", values_to = "log_size")
  
  surv_clim_seqs <- make_clim_seqs(surv_terms, ramet_temp) %>%
    cbind(site = pop)
  
  grow_clim_seqs <- make_clim_seqs(grow_terms, ramet_temp) %>%
    cbind(site = pop)
  
  repr_clim_seqs <- make_clim_seqs(repr_terms, ramet_temp) %>%
    cbind(site = pop)
  
  flow_clim_seqs <- make_clim_seqs(flow_terms, ramet_temp) %>%
    cbind(site = pop)
  
  pdf(glue('ipms/Figures/vr_models/gam/{pop}_p_kern_preliminary_plots.pdf'),
      width = 8,
      height = 8)
  
    clim_plots(sizes, surv_clim_seqs, use_surv_mod, surv_terms)
    clim_plots(sizes, grow_clim_seqs, use_grow_mod, grow_terms)
    
  dev.off()
  
  pdf(glue('ipms/Figures/vr_models/gam/{pop}_f_kern_preliminary_plots.pdf'),
      width = 8,
      height = 8)
  
    clim_plots(sizes, repr_clim_seqs, use_repr_mod, repr_terms)
    clim_plots(sizes, flow_clim_seqs, use_flow_mod, flow_terms)
  
  dev.off()
  
  
}

sizes <- c(all_ramets$log_size, all_ramets$log_size_next)

use_sites <- c("Havatselet", "Foxton", "Rooisand", "Colares", "Rarangi")
all_sizes <- expand.grid(site = use_sites,
                         log_size = seq(min(sizes, na.rm = TRUE),
                                        max(sizes, na.rm = TRUE),
                                        length.out = 50)) %>%
  mutate(native = case_when(
    site %in% c("Rooisand") ~ 1,
    TRUE ~ 0
  ))


pdf(glue('ipms/Figures/vr_models/gam/all_pops_f_kern_preliminary_plots.pdf'),
    width = 8,
    height = 8)

  size_clim_plots(all_sizes, clim_data, use_repr_mod, repr_terms, use_sites)
  size_clim_plots(all_sizes, clim_data, use_flow_mod, flow_terms, use_sites)

dev.off()

pdf(glue('ipms/Figures/vr_models/gam/all_pops_p_kern_preliminary_plots.pdf'),
    width = 8,
    height = 8)

  size_clim_plots(all_sizes, clim_data, use_surv_mod, surv_terms, use_sites)
  size_clim_plots(all_sizes, clim_data, use_grow_mod, grow_terms, use_sites)

dev.off()
