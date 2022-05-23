# LTRE Calculations: Random Forest with complete posterior + lambdas
all_surv_mods <- readRDS("ipms/Model_Fits/ramet_survival_list_krig.rds")
all_flow_mods <- readRDS("ipms/Model_Fits/ramet_flower_n_list_krig.rds")

surv_mod <- all_surv_mods[[1]][["times_sw2_seas"]]
flow_mod <- all_flow_mods[[1]][["times_sw2_ann"]]

grow_mod <- readRDS("ipms/Model_Fits/grow_2_1_lin_gam_mix.rds")
repr_mod <- readRDS("ipms/Model_Fits/repro_lin_gam_mix.rds")
recr_mod <- readRDS("ipms/Model_Fits/recr_size_brm.rds")

surv_draws <- as.data.frame(surv_mod) %>%
  select(c(1:16, 19:31)) %>%
  rename_draws("surv") 
names(surv_draws) <- gsub(pattern = "^surv_z$", 
                          "surv_log_size", 
                          x = names(surv_draws))
surv_obs <- summarise_all(surv_draws, mean)

flow_draws <- as.data.frame(flow_mod) %>%
  select(c(1:16, 20:32)) %>%
  rename_draws("flow") 
names(flow_draws) <- gsub(pattern = "^flow_z$", 
                          "flow_log_size", 
                          x = names(flow_draws))
flow_obs <- summarise_all(flow_draws, mean)

repr_draws <- as.data.frame(repr_mod) %>%
  select(1:10, 15:33) %>%
  rename_draws("repr")

repr_obs <- summarise_all(repr_draws, mean)

grow_draws <- as.data.frame(grow_mod) %>%
  select(1:10, 20:99) %>%
  rename_draws("grow")

# Replace the 11-13's on the trailing end of the final smooth names for growth
names(grow_draws)[85:90] <- gsub("t_1", "t_", names(grow_draws)[85:90]) 

grow_obs <- summarise_all(grow_draws, mean)

recr_draws <- as.data.frame(recr_mod) %>%
  select(1:2) %>%
  setNames(c("recr_mean", "recr_sd"))

recr_obs <- summarise_all(recr_draws, mean)

lambdas <- read.csv("ipms/Model_Fits/ipms/gam_mod_general_ipm_site_lambdas.csv",
                    stringsAsFactors = FALSE) %>%
  filter(obs == "no") %>%
  select(-obs)

rf_list <- vector("list", ncol(lambdas))

start <- Sys.time()

for(i in seq_along(lambdas)) {
  
  pop_nm <- names(lambdas)[i] %>%
    gsub(pattern = "^lambda_", replacement = "", x = .)
  
  pop_ran <- paste0(c("grow", "grow_sigma", "surv", "flow", "repr"), "_intercept_", pop_nm)
  
  temp_surv <- surv_draws[, c(1:16, grep(pop_ran[3], names(surv_draws)))]
  temp_flow <- flow_draws[, c(1:16, grep(pop_ran[4], names(flow_draws)))] 
  temp_repr <- repr_draws[, c(1:10, grep(pop_ran[5], names(repr_draws)), 24:29)]
  
  grow_ind  <- c(1:10, 
                 grep(pop_ran[1], names(grow_draws)),
                 grep(pop_ran[2], names(grow_draws)),
                 37:90)
  
  temp_grow <- grow_draws[ , grow_ind]
    
  pop_draws <- cbind(temp_surv, temp_grow, temp_repr, temp_flow, recr_draws)

  lams   <- lambdas[ , i]
  
  tuned <- tuneRF(x          = pop_draws, 
                  y          = lams,
                  plot       = FALSE,
                  ntreeTry   = 500,
                  stepFactor = 2, 
                  improve    = 1e-6,
                  mtryStart  = 100)
  
  optimal_mtry <- which.min(tuned[, 2])
  
  temp_rf <- randomForest(x          = pop_draws,
                          y          = lams,
                          ntree      = 500, 
                          importance = TRUE,
                          mtry       = tuned[optimal_mtry, 1])
  
  sink(paste0("ipms/Model_Summaries/ltre/", pop_nm, "_ltre.txt"))
  
    cat("Tuning Results -----------\n\n")
    print(tuned)
    cat("\n\nRandom Forest results-----------\n\n")
    print(temp_rf)
    cat("\n\n END -----------\n")
  sink()
  
  imps <- temp_rf$importance %>%
    as.data.frame() %>%
    mutate(var = rownames(.)) %>%
    select(c(1, 3)) %>%
    setNames(c("inc_mse", "var")) %>%
    arrange(desc(inc_mse))
  
  # Save top 15 most important variables for plotting
  
  if(i == 1) {
    
    out <- data.frame(site = pop_nm,
                      var  = imps[1:15, 2],
                      imp  = imps[1:15, 1])
    out <- cbind(out, out_sd = temp_rf$importanceSD[out$var])
    
    var_exp <- data.frame(site = pop_nm,
                          val  = temp_rf$rsq[length(temp_rf$rsq)] * 100)
    
  } else {
    
    temp_out <- data.frame(site = pop_nm,
                           var  = imps[1:15, 2],
                           imp  = imps[1:15, 1])
    temp_out <- cbind(temp_out, out_sd = temp_rf$importanceSD[temp_out$var])
    
    out <- rbind(out,
                 temp_out)
    
    var_exp <- rbind(var_exp,
                     data.frame(site = pop_nm,
                                val  = temp_rf$rsq[length(temp_rf$rsq)] * 100))
  }
  
  rf_list[[i]] <- temp_rf
  names(rf_list)[i] <- pop_nm
  
  
}

dif <- Sys.time() - start

dif

var_exp$native <- ifelse(var_exp$site %in% c("Struisbaai", "Rooisand",
                                             "Vogelgat", "Springfontein",
                                             "Melkboss", "St_Francis"),
                         "Native",
                         "Invaded")

saveRDS(rf_list, file = "ipms/Model_Fits/ltre/all_random_forests.rds")

saveRDS(out, 
        file = "ipms/Model_Fits/ltre/random_forest_importance_table.rds")
saveRDS(var_exp,
        file = "ipms/Model_Fits/ltre/variance_explained.rds")

pdf("ipms/Figures/ipm_output/LTRE_RF.pdf")
  for(i in unique(out$site)) {
    
    use_dat <- data.frame(x = 1:500, y = rf_list[[i]]$rsq * 100)
    
    r2_plt <- ggplot(use_dat, aes(x = x, y = y)) + 
      geom_line(color = "red", size = 1.5) +
      theme_bw() +
      ylab("% Variance Explained") +
      xlab("# of Trees")
    
    temp <- filter(out, site == i)
    vi_plt <- ggplot(temp, aes(x = var, y = imp)) + 
      geom_point(size = 2) + 
      geom_linerange(aes(ymin = imp - 2 * out_sd, ymax = imp + 2 * out_sd),
                     size = 1.2) +
      coord_flip() +
      theme_bw() +
      ggtitle(i) +
      ylab("Importance") +
      xlab("Vital Rate Regression Parameter")
    
    
    print((vi_plt / r2_plt))
  }

  print(
    ggplot(var_exp, aes(x = site, y = val, color = native, fill = native)) + 
      geom_bar(stat = "identity") + 
      theme_bw() +
      ylab("% Variance Explained by Random Forest") +
      xlab("Site") +
      scale_discrete_manual(breaks = c("Native", "Invaded"),
                            values = viridis::inferno(2, begin = 0.2, end = 0.5),
                            aesthetics = c("fill", "color")) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16))
    
  )
dev.off()
