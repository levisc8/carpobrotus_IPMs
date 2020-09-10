# Sensitivity and elasticity

# First, do kernel sensitivities. Function definitions are 
# in 01_Utils_and_Dependencies

K_exp_var <- carp_ipmr$iterators$K
P_exp_var <- carp_ipmr$sub_kernels$P
F_exp_var <- carp_ipmr$sub_kernels$F

bs <- seq(L, U, length.out = n_mesh_p + 1)
mesh <- d1 <- (bs[2:(n_mesh_p + 1)] + bs[1:n_mesh_p]) * 0.5

h <- bs[2] - bs[1]

k_sens_exp_mat   <- sensitivity(K_exp_var, h, level = 'kernel') 

# Elasticities

k_elas_exp_mat   <- elasticity(K_exp_var, h, level = 'kernel') 


# Sub-kernel elasticities

P_fun_ev <- P_exp_var / h
F_fun    <- F_exp_var / h

# Individual kernel elasticities
P_elas_ev_mat <- P_fun_ev * k_sens_exp_mat / lambda_exp_var

F_elas_ev_mat <- F_fun * k_sens_exp_mat / lambda_exp_var


# Transform for ggplot'ing

mesh_df <- expand.grid(cols = mesh,
                       rows = mesh)

k_sens_exp <- k_sens_exp_mat %>%
  mat_to_df() %>%
  cbind(mesh_df) %>%
  mutate(Model = "Exponential Variance") %>%
  setNames(
    c(
      'value', 'x', 'y', 'Model'
    )
  )


k_elas_exp <- k_elas_exp_mat %>%
  mat_to_df() %>%
  cbind(mesh_df) %>%
  mutate(Model = "Exponential Variance") %>%
  setNames(
    c(
      'value', 'x', 'y', 'Model'
    )
  )


P_elas_ev <- P_elas_ev_mat %>%
  mat_to_df() %>%
  cbind(mesh_df) %>%
  mutate(Model = "Exponential Variance") %>%
  setNames(
    c(
      'value', 'x', 'y', 'Model'
    )
  )

F_elas_ev <- F_elas_ev_mat %>%
  mat_to_df() %>%
  cbind(mesh_df) %>%
  mutate(Model = "Exponential Variance") %>%
  setNames(
    c(
      'value', 'x', 'y', 'Model'
    )
  )


# Sensitivity of lambda to flower_n


flow_unc_lambda <- rep(NA_real_, 1000)


unc_param_list <- all_param_list
unc_param_list$f_r <- NULL
unc_param_list$n_new  <- n_new

# Function to generate a named list of new parameter values at every
# iteration of the stochastic model

par_sampler <- function(int_mu, int_sd, slope_mu, slope_sd) {
  
  out <- list(
    f_s_int_stoch   = rnorm(1, int_mu, int_sd),
    f_s_slope_stoch = rnorm(1, slope_mu, slope_sd)
  )
  
  return(out)
  
}

for(i in seq_len(1000)) {
  
  unc_temp  <- par_sampler(unc_param_list$f_s_int,
                           unc_param_list$f_s_int_sd,
                           unc_param_list$f_s_slope,
                           unc_param_list$f_s_slope_sd)
  
  unc_param_list <- c(unc_param_list, unc_temp)
  
  carp_ipmr <- init_ipm('simple_di_det') %>%
    define_kernel(
      name      = "P",
      formula   = S * G,
      family    = "CC",
      G         = dnorm(sa_2, mu_g, sd_g),
      sd_g      = sqrt(exp(2 * g_sigma_par * sa_1)),
      mu_g      = g_int + g_slope * sa_1,
      S         = inv_logit(s_int, s_slope, sa_1),
      data_list = unc_param_list,
      states    = list(c('sa')),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>%
    define_kernel(
      name      = "F",
      formula   = f_r_stoch * f_s * f_d * p_r,
      family    = "CC",
      
      # Divide by 100 because there will be 100 duplicated values for every level of sa_1
      
      f_r_stoch = n_new / (sum(f_s) / 100), 
      f_s       = exp(f_s_int_stoch + f_s_slope_stoch * sa_1),
      f_d       = dnorm(sa_2, f_d_mu, f_d_sd),
      p_r       = inv_logit(p_r_int, p_r_slope, sa_1),
      data_list = unc_param_list,
      states    = list(c("sa")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "f_d")
    ) %>%
    define_k(
      name      = "K",
      family    = "IPM",
      K         = P + F,
      n_sa_t_1  = K %*% n_sa_t,
      data_list = unc_param_list,
      states    = list(c('sa')),
      evict_cor = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(  
        kernel_names = c("K", "P", "F"),
        int_rule     = rep('midpoint', 3),
        dom_start    = rep('sa', 3),
        dom_end      = rep('sa', 3) 
      ) 
    ) %>%
    define_domains(
      sa = c(L, U, n_mesh_p)
    )  %>%
    define_pop_state(
      n_sa = rep(1/100, 100)
    ) %>% 
    make_ipm(usr_funs   = list(inv_logit = s_z),
             iterate    = TRUE,
             iterations = 100L)
  
  # Exponentiate for comparison to deterministic lambda
  flow_unc_lambda[i] <- lambda(carp_ipmr)
  
  unc_param_list <- unc_param_list[!grepl("stoch", names(unc_param_list))]
  
  if(i %% 100 == 0) message("\n", i, "/1000 iterations done.\n")
}

flow_unc_lambda <- sort(flow_unc_lambda)

unc_plot <- data.frame(obs = lambda_ipmr,
                       up_ci = flow_unc_lambda[975],
                       lo_ci = flow_unc_lambda[25])

