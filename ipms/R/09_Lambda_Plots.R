# Plot IPM Outputs

all_lams <- read.csv("ipms/Model_Fits/ipms/lin_mod_site_lambdas.csv",
                     stringsAsFactors = FALSE) %>%
  pivot_longer(cols = -c(obs),
               names_to = "site",
               values_to = "lambda") %>%
  mutate(site = gsub("lambda_", "", site))

clim_vars <- expand.grid(c("temp", "prec", "sw1", "sw3"), c("dry", "wet"), c("t", "t_1")) %>%
  apply(1, function(x) paste(x, collapse = "_")) %>% 
  c("site") %>%
  syms()


all_clim <- read.csv("ipms/Model_Fits/ipms/site_climate_values.csv",
                     stringsAsFactors = FALSE) %>%
  select(!!! clim_vars) %>%
  pivot_longer(-site,
               names_to = "clim_var", 
               values_to = "clim_val")

all_data <- all_lams %>%
  group_by(site) %>%
  summarise(obs_val  = first(lambda[obs == "yes"]),
            up_ci    = quantile(lambda[obs == "no"], 0.975),
            lo_ci    = quantile(lambda[obs == "no"], 0.025)) %>%
  left_join(all_clim) %>%
  mutate(native = case_when(
    site %in% c("Struisbaai", "Rooisand",
                "Vogelgat", "Springfontein",
                "Melkboss", "St_Francis") ~ "Native",
    TRUE ~ "Invasive"
  ))


plts <- list()

for(i in seq_along(unique(all_data$clim_var))) {
  
  clim <- unique(all_data$clim_var)[i]
  
  temp_dat <- filter(all_data, clim_var == clim) 
  
  plts[[i]] <- ggplot(temp_dat, aes(x = clim_val, y = obs_val)) +
    geom_point(aes(shape = native, color = site), size = 5) +
    geom_linerange(aes(ymin = lo_ci, ymax = up_ci, color = site)) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, k = 10),
                color = "black",
                linetype = "dashed") +
    #geom_smooth(method = "lm", color = "red", linetype = "dotted") + 
    theme_bw() +
    ggtitle(clim)
  
}

pdf(file = "ipms/Figures/ipm_output/lambdas_by_clim.pdf",
    height = 9, width = 9)
for(i in 1:4) {
  
  ind <- seq(i * 4 - 3, i * 4, 1)
  
  temp_plts <- plts[ind]
  
  print(((temp_plts[[1]] + temp_plts[[2]]) / 
      (temp_plts[[3]] + temp_plts[[4]])) + 
        plot_layout(guides = "collect"))
  
}

dev.off()

plts <- list()

for(i in seq_along(unique(all_data$clim_var))) {
  
  clim <- unique(all_data$clim_var)[i]
  
  temp_dat <- filter(all_data, clim_var == clim & site != "Havatselet") 
  
  plts[[i]] <- ggplot(temp_dat, aes(x = clim_val, y = obs_val)) +
    geom_point(aes(shape = native, color = site), size = 5) +
    geom_linerange(aes(ymin = lo_ci, ymax = up_ci, color = site)) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, k = 10),
                color = "black",
                linetype = "dashed") +
    #geom_smooth(method = "lm", color = "red", linetype = "dotted") + 
    theme_bw() +
    ggtitle(clim)
  
}


pdf(file = "ipms/Figures/ipm_output/lambdas_by_clim_sans_israel.pdf",
    height = 9, width = 9)
for(i in 1:4) {
  
  ind <- seq(i * 4 - 3, i * 4, 1)
  
  temp_plts <- plts[ind]
  
  print(((temp_plts[[1]] + temp_plts[[2]]) / 
           (temp_plts[[3]] + temp_plts[[4]])) + 
          plot_layout(guides = "collect"))
  
}

dev.off()

# GAM IPM Outputs -----------

# Plot IPM Outputs

all_lams <- read.csv("ipms/Model_Fits/ipms/gam_mod_site_lambdas.csv",
                     stringsAsFactors = FALSE) %>%
  pivot_longer(cols = -c(obs),
               names_to = "site",
               values_to = "lambda") %>%
  mutate(site = gsub("lambda_", "", site))

clim_vars <- expand.grid(c("temp", "prec", "sw1", "sw3"),
                         c("dry", "wet"), 
                         c("t", "t_1")) %>%
  apply(1, function(x) paste(x, collapse = "_")) %>% 
  c("site") %>%
  syms()


all_clim <- read.csv("ipms/Model_Fits/ipms/site_climate_values.csv",
                     stringsAsFactors = FALSE) %>%
  select(!!! clim_vars) %>%
  pivot_longer(-site,
               names_to = "clim_var", 
               values_to = "clim_val")

all_data <- all_lams %>%
  group_by(site) %>%
  summarise(obs_val  = first(lambda[obs == "yes"]),
            up_ci    = quantile(lambda[obs == "no"], 0.975),
            lo_ci    = quantile(lambda[obs == "no"], 0.025)) %>%
  left_join(all_clim) %>%
  mutate(native = case_when(
    site %in% c("Struisbaai", "Rooisand",
                "Vogelgat", "Springfontein",
                "Melkboss", "St_Francis") ~ "Native",
    TRUE ~ "Invasive"
  ))


plts <- list()

for(i in seq_along(unique(all_data$clim_var))) {
  
  clim <- unique(all_data$clim_var)[i]
  
  temp_dat <- filter(all_data, clim_var == clim) 
  
  plts[[i]] <- ggplot(temp_dat, aes(x = clim_val, y = obs_val)) +
    geom_point(aes(shape = native, color = site), size = 5) +
    geom_linerange(aes(ymin = lo_ci, ymax = up_ci, color = site)) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, k = 10),
                color = "black",
                linetype = "dashed") +
    #geom_smooth(method = "lm", color = "red", linetype = "dotted") + 
    theme_bw() +
    ggtitle(clim)
  
}

pdf(file = "ipms/Figures/ipm_output/lambdas_by_clim_gam.pdf",
    height = 9, width = 9)
for(i in 1:4) {
  
  ind <- seq(i * 4 - 3, i * 4, 1)
  
  temp_plts <- plts[ind]
  
  print(((temp_plts[[1]] + temp_plts[[2]]) / 
           (temp_plts[[3]] + temp_plts[[4]])) + 
          plot_layout(guides = "collect"))
  
}

dev.off()

# plts <- list()
# 
# for(i in seq_along(unique(all_data$clim_var))) {
#   
#   clim <- unique(all_data$clim_var)[i]
#   
#   temp_dat <- filter(all_data, clim_var == clim & site != "Havatselet") 
#   
#   plts[[i]] <- ggplot(temp_dat, aes(x = clim_val, y = obs_val)) +
#     geom_point(aes(shape = native, color = site), size = 5) +
#     geom_linerange(aes(ymin = lo_ci, ymax = up_ci, color = site)) +
#     geom_smooth(method = "gam",
#                 formula = y ~ s(x, k = 10),
#                 color = "black",
#                 linetype = "dashed") +
#     #geom_smooth(method = "lm", color = "red", linetype = "dotted") + 
#     theme_bw() +
#     ggtitle(clim)
#   
# }
# 
# 
# pdf(file = "ipms/Figures/ipm_output/lambdas_by_clim_sans_israel.pdf",
#     height = 9, width = 9)
# for(i in 1:4) {
#   
#   ind <- seq(i * 4 - 3, i * 4, 1)
#   
#   temp_plts <- plts[ind]
#   
#   print(((temp_plts[[1]] + temp_plts[[2]]) / 
#            (temp_plts[[3]] + temp_plts[[4]])) + 
#           plot_layout(guides = "collect"))
#   
# }
# 
# dev.off()



all_lams <- read.csv("ipms/Model_Fits/ipms/gam_mod_general_ipm_site_lambdas.csv",
                     stringsAsFactors = FALSE) %>%
  pivot_longer(cols = -c(obs),
               names_to = "site",
               values_to = "lambda") %>%
  mutate(site = gsub("lambda_", "", site))

clim_vars <- expand.grid(c("temp", "prec", "sw1", "sw3"),
                         c("dry", "wet"), 
                         c("t", "t_1")) %>%
  apply(1, function(x) paste(x, collapse = "_")) %>% 
  c("site") %>%
  syms()


all_clim <- read.csv("ipms/Model_Fits/ipms/site_climate_values.csv",
                     stringsAsFactors = FALSE) %>%
  select(!!! clim_vars) %>%
  pivot_longer(-site,
               names_to = "clim_var", 
               values_to = "clim_val")

all_data <- all_lams %>%
  group_by(site) %>%
  summarise(obs_val  = first(lambda[obs == "yes"]),
            up_ci    = quantile(lambda[obs == "no"], 0.975),
            lo_ci    = quantile(lambda[obs == "no"], 0.025)) %>%
  left_join(all_clim) %>%
  mutate(native = case_when(
    site %in% c("Struisbaai", "Rooisand",
                "Vogelgat", "Springfontein",
                "Melkboss", "St_Francis") ~ "Native",
    TRUE ~ "Invasive"
  ))


plts <- list()

for(i in seq_along(unique(all_data$clim_var))) {
  
  clim <- unique(all_data$clim_var)[i]
  
  temp_dat <- filter(all_data, clim_var == clim) 
  
  plts[[i]] <- ggplot(temp_dat, aes(x = clim_val, y = obs_val)) +
    geom_point(aes(shape = native, color = site), size = 5) +
    geom_linerange(aes(ymin = lo_ci, ymax = up_ci, color = site)) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, k = 10),
                color = "black",
                linetype = "dashed") +
    #geom_smooth(method = "lm", color = "red", linetype = "dotted") + 
    theme_bw() +
    ggtitle(clim)
  
}

pdf(file = "ipms/Figures/ipm_output/lambdas_by_clim_gam_general.pdf",
    height = 9, width = 9)

for(i in 1:4) {
  
  ind <- seq(i * 4 - 3, i * 4, 1)
  
  temp_plts <- plts[ind]
  
  print(((temp_plts[[1]] + temp_plts[[2]]) / 
           (temp_plts[[3]] + temp_plts[[4]])) + 
          plot_layout(guides = "collect"))
  
}

dev.off()
