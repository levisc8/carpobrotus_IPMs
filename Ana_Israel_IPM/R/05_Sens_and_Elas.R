# Sensitivity and elasticity

# First, do kernel sensitivities. Function definitions are 
# in 01_Utils_and_Dependencies

k_sens_exp_mat   <- sensitivity(K_exp_var, h, level = 'kernel') 

# Elasticities

k_elas_exp_mat   <- elasticity(K_exp_var, h, level = 'kernel') 

# Sub-kernel elasticities

P_fun_ev <- P_exp_var / h
F_fun    <- Fm / h

# Individual kernel elasticities
P_elas_ev_mat <- P_fun_ev * k_sens_exp_mat / lambda_exp_var

F_elas_ev_mat <- F_fun * k_sens_exp_mat / lambda_exp_var


# Transform for ggplot'ing

k_sens_exp <- k_sens_exp_mat %>%
  mat_to_df(meshp = d1) %>%
  mutate(Model = "Exponential Variance")

k_elas_exp <- k_elas_exp_mat %>%
  mat_to_df(meshp = d1) %>%
  mutate(Model = "Exponential Variance")


P_elas_ev <- P_elas_ev_mat %>%
  mat_to_df(meshp = d1) %>%
  mutate(Model = "Exponential Variance")

F_elas_ev <- F_elas_ev_mat %>%
  mat_to_df(meshp = d1) %>%
  mutate(Model = "Exponential Variance")

