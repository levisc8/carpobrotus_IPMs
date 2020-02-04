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


