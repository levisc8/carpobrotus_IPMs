# Bootstrapping code

# Split the data into plants that existed at T and new recruits at T+1

exists_t_1  <- filter(all_data, !is.na(log_size))
recruit_t_2 <- filter(all_data, is.na(log_size) & id > 8000)

# Number of existing plants to resample

n_resamp_existing <- dim(exists_t_1)[1]
n_resamp_recruits <- dim(recruit_t_2)[1]

n_resamp <- 1000


# Set up place to hold boot strapping outputs

temp_exp_var_out <- splice(grow_exp_var_coef_list,
                           surv_coef_list,
                           fec_coef_list,
                           lambda = lambda_exp_var,
                           p_elas = sum(P_elas_ev_mat) * h ^ 2,
                           f_elas = sum(F_elas_ev_mat) * h ^ 2)

exp_var_out <- lapply(temp_exp_var_out,
                      function(x, n_resamp) {
                        c(x, rep(NA_real_, n_resamp))
                      },
                      n_resamp = n_resamp)

for(i in seq_len(n_resamp)) {
  
  # Resample our data with replacement and create a single data frame
  
  ex_resamp_ind <- sample(1:n_resamp_existing,
                          n_resamp_existing, 
                          replace = TRUE)
  re_resamp_ind <- sample(1:n_resamp_recruits,
                          n_resamp_recruits, 
                          replace = TRUE)
  
  boot_existing <- exists_t_1[ex_resamp_ind, ]
  boot_recruits <- recruit_t_2[re_resamp_ind, ]
  
  boot_data <- rbind(boot_existing, boot_recruits)
  
  # Bootstrap growth model
  
  boot_grow_mod_exp <- gls(log_size_next ~ log_size, data = boot_data,
                           weights = varExp(),
                           na.action = na.omit,
                           method = 'ML')
  
  # Bootstrap survival model
  
  boot_surv_mod <- glm(survival ~ log_size, data = boot_data, family = binomial())
  
  # Bootstrap fecundity parameters
  
  n_flow <- sum(boot_data$flower_n, na.rm = TRUE)
  
  # Not computing a number of new recruits because we've constrained that in our
  # procedure. Thus, only the number of flowers can vary from sample to sample
  
  boot_f_r <- n_new / n_flow
  
  boot_p_r_mod <- glm(repro ~ log_size, data = boot_data, family = binomial())
  
  boot_f_d_mu <- mean(boot_recruits$log_size_next, na.rm = TRUE)
  boot_f_d_sd <- sd(boot_recruits$log_size_next, na.rm = TRUE)
  
  boot_f_s_mod <- glm(flower_n ~ log_size, 
                      data = boot_data, 
                      family = quasipoisson())
  
  # We already have the implementation details from when we generated the point
  # estimate of lambda, so no need to redefine those.
  
   P_exp_var_boot <- outer(d1, d2,
                          FUN = p_ev_z1z,
                          s_int    = coef(boot_surv_mod)[1],
                          s_slope  = coef(boot_surv_mod)[2],
                          g_int    = coef(boot_grow_mod_exp)[1],
                          g_slope  = coef(boot_grow_mod_exp)[2],
                          var_coef = as.numeric(boot_grow_mod_exp$modelStruct$varStruct),
                          L = L,
                          U = U,
                          h = h) %>%
    t()
  
  
  Fm_boot <- outer(d1, d2,
                   FUN = f_z1z,
                   pr_int = coef(boot_p_r_mod)[1],
                   pr_slope = coef(boot_p_r_mod)[2],
                   f_s_int = coef(boot_f_s_mod)[1],
                   f_s_slope = coef(boot_f_s_mod)[2],
                   fr = boot_f_r,
                   fd_mu = boot_f_d_mu,
                   fd_sd = boot_f_d_sd,
                   s_int = coef(boot_surv_mod)[1],
                   s_slope = coef(boot_surv_mod)[2],
                   L = L,
                   h = h) %>%
    t()
  
  K_exp_var_boot   <- (P_exp_var_boot   + Fm_boot)
  
  
  # Store lambdas
  exp_var_out$lambda[(i + 1)]   <- l_ev_boot <- Re(eigen(K_exp_var_boot)$values[1])
  
  k_sens_ev_boot <- sensitivity(K_exp_var_boot, h, level = 'kernel')
  
  P_fun_ev_boot <- P_exp_var_boot / h
  F_fun_boot    <- Fm_boot / h
  
  p_elas_ev_mat_boot <- P_fun_ev_boot * k_sens_ev_boot / l_ev_boot
  f_elas_ev_mat_boot <- F_fun_boot * k_sens_ev_boot / l_ev_boot
  
  exp_var_out$p_elas[(i + 1)] <- sum(p_elas_ev_mat_boot) * h^2
  exp_var_out$f_elas[(i + 1)] <- sum(f_elas_ev_mat_boot) * h^2
  
  # Store exponential variance growth parameters
  exp_var_out$g_int[(i + 1)]       <- coef(boot_grow_mod_exp)[1]
  exp_var_out$g_slope[(i + 1)]     <- coef(boot_grow_mod_exp)[2]  
  exp_var_out$g_sigma_par[(i + 1)] <- as.numeric(boot_grow_mod_exp$modelStruct$varStruct)
  
  # Survival models are the same across both 
  exp_var_out$s_int[(i + 1)]   <- coef(boot_surv_mod)[1]
  exp_var_out$s_slope[(i + 1)] <- coef(boot_surv_mod)[2] 
  
  # Same with pr(repro) and the other fecundity parameters
  exp_var_out$p_r_int[(i + 1)]   <- coef(boot_p_r_mod)[1]
  exp_var_out$p_r_slope[(i + 1)] <- coef(boot_p_r_mod)[2] 
  
  exp_var_out$f_s_int[(i + 1)]   <- coef(boot_f_s_mod)[1]
  exp_var_out$f_s_slope[(i + 1)] <- coef(boot_f_s_mod)[2]
  
  exp_var_out$f_r[(i + 1)]    <- boot_f_r
  exp_var_out$f_d_mu[(i + 1)] <- boot_f_d_mu
  exp_var_out$f_d_sd[(i + 1)] <- boot_f_d_sd
  
}


