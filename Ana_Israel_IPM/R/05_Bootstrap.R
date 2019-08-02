# Bootstrapping code

# Split the data into plants that existed at T and new recruits at T+1

exists_t_1  <- filter(all_data, !is.na(log_size))
recruit_t_2 <- filter(all_data, is.na(log_size) & id > 8000)

# Number of existing plants to resample

n_resamp_existing <- dim(exists_t_1)[1]
n_resamp_recruits <- dim(recruit_t_2)[1]

n_resamp <- 1000


# Set up place to hold boot strapping outputs

temp_const_var_out <- splice(grow_const_var_coef_list,
                             surv_coef_list,
                             fec_coef_list,
                             lambda = lambda_const_var)

const_var_out <- lapply(temp_const_var_out,
                        function(x, n_resamp) {
                          c(x, rep(NA_real_, n_resamp))
                        },
                        n_resamp = n_resamp)
temp_exp_var_out <- splice(grow_exp_var_coef_list,
                           surv_coef_list,
                           fec_coef_list,
                           lambda = lambda_exp_var)

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
  
  boot_grow_mod_lin <- lm(log_size_next ~ log_size, data = boot_data)
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
  
  boot_f_s_mod <- glm(flower_n ~ log_size, data = boot_data, family = poisson())
  
  # We already have the implementation details from when we generated the point
  # estimate of lambda, so no need to redefine those.
  
  P_const_var <- p_cv_z1z(d1, d2,
                          coef(boot_surv_mod)[1],
                          coef(boot_surv_mod)[2],
                          coef(boot_grow_mod_lin)[1],
                          coef(boot_grow_mod_lin)[2],
                          sd(resid(boot_grow_mod_lin)),
                          h) 
  
  P_exp_var <- p_ev_z1z(d1, d2,
                        coef(boot_surv_mod)[1],
                        coef(boot_surv_mod)[2],
                        coef(boot_grow_mod_exp)[1],
                        coef(boot_grow_mod_exp)[2],
                        as.numeric(boot_grow_mod_exp$modelStruct$varStruct),
                        h) 
  
  
  Fm <- outer(d1, d2,
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
              h = h)
  
  K_const_var <- P_const_var + Fm
  K_exp_var   <- P_exp_var   + Fm
  
  # Store lambdas
  const_var_out$lambda[(i + 1)] <- Re(eigen(K_const_var)$values[1])
  exp_var_out$lambda[(i + 1)]   <- Re(eigen(K_exp_var)$values[1])
  
  # Store constant variace growth parameters
  const_var_out$g_int[(i + 1)]    <- coef(boot_grow_mod_lin)[1]
  const_var_out$g_slope[(i + 1)]  <- coef(boot_grow_mod_lin)[2] 
  const_var_out$sd_g[(i + 1)]     <- sd(resid(boot_grow_mod_lin))
  
  # Store exponential variance growth parameters
  exp_var_out$g_int[(i + 1)]       <- coef(boot_grow_mod_exp)[1]
  exp_var_out$g_slope[(i + 1)]     <- coef(boot_grow_mod_exp)[2]  
  exp_var_out$g_sigma_par[(i + 1)] <- as.numeric(boot_grow_mod_exp$modelStruct$varStruct)
  
  # Survival models are the same across both 
  const_var_out$s_int[(i + 1)]   <- exp_var_out$s_int[(i + 1)]   <- coef(boot_surv_mod)[1]
  const_var_out$s_slope[(i + 1)] <- exp_var_out$s_slope[(i + 1)] <- coef(boot_surv_mod)[2] 
  
  # Same with pr(repro) and the other fecundity parameters
  const_var_out$p_r_int[(i + 1)]   <- exp_var_out$p_r_int[(i + 1)]   <- coef(boot_p_r_mod)[1]
  const_var_out$p_r_slope[(i + 1)] <- exp_var_out$p_r_slope[(i + 1)] <- coef(boot_p_r_mod)[2] 
  
  const_var_out$f_s_int[(i + 1)]   <- exp_var_out$f_s_int[(i + 1)]   <- coef(boot_f_s_mod)[1]
  const_var_out$f_s_slope[(i + 1)] <- exp_var_out$f_s_slope[(i + 1)] <- coef(boot_f_s_mod)[2]
  
  const_var_out$f_r[(i + 1)]    <- exp_var_out$f_r[(i + 1)]    <- boot_f_r
  const_var_out$f_d_mu[(i + 1)] <- exp_var_out$f_d_mu[(i + 1)] <- boot_f_d_mu
  const_var_out$f_d_sd[(i + 1)] <- exp_var_out$f_d_sd[(i + 1)] <- boot_f_d_sd
  
}


