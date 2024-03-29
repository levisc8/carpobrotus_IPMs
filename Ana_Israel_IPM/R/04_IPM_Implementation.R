
s_z <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

L <- min(c(all_data$log_size,
           all_data$log_size_next),
         na.rm = TRUE) * 1.2

U <- max(c(all_data$log_size,
           all_data$log_size_next),
         na.rm = TRUE) * 1.2

n_mesh_p <- 100

carp_ipmr <- init_ipm("simple", "di", "det") %>%
  define_kernel(
    name = "P",
    formula = S * G,
    family = "CC",
    G = dnorm(sa_2, mu_g, sd_g),
    sd_g = sqrt(exp(2 * g_sigma_par * sa_1)),
    mu_g = g_int + g_slope * sa_1,
    S = inv_logit(s_int, s_slope, sa_1),
    data_list = all_param_list,
    states = list(c('sa')),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "G")
  ) %>%
  define_kernel(
    name = "F",
    formula = f_r * f_s * f_d * p_r,
    family = "CC",
    f_s = exp(f_s_int + f_s_slope * sa_1),
    f_d = dnorm(sa_2, f_d_mu, f_d_sd),
    p_r = inv_logit(p_r_int, p_r_slope, sa_1),
    f_r = n_new / sum(f_s),
    data_list = all_param_list,
    states = list(c("sa")),
    evict_cor = TRUE,
    evict_fun = truncated_distributions("norm", "f_d")
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"),
      int_rule = rep('midpoint', 2),
      state_start = rep('sa', 2),
      state_end = rep('sa', 2) 
    ) 
  ) %>%
  define_domains(
    sa = c(L, U, n_mesh_p)
  )  %>%
  define_pop_state(
    n_sa = rep(1/100, n_mesh_p)
  ) %>%
  make_ipm(usr_funs = list(inv_logit = s_z),
           iterate = TRUE,
           iterations = 100,
           return_all_envs = TRUE)

lambda_ipmr <- lambda_exp_var <- lambda(carp_ipmr)

# Iterate over meshpoints to make sure matrix size doesn't matter

mesh_p <- seq(100, 500, by = 50)

lambdas <- data.frame(mesh_p = mesh_p,
                      lambdas = c(lambda_ipmr, rep(NA_real_, 8)))

it <- 1

for(i in mesh_p) {
  

  carp_ipmr_test <- init_ipm("simple", "di", "det") %>%
    define_kernel(
      name = "P",
      formula = S * G,
      family = "CC",
      G = dnorm(sa_2, mu_g, sd_g),
      sd_g = sqrt(exp(2 * g_sigma_par * sa_1)),
      mu_g = g_int + g_slope * sa_1,
      S = inv_logit(s_int, s_slope, sa_1),
      data_list = all_param_list,
      states = list(c('sa')),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "G")
    ) %>%
    define_kernel(
      name = "F",
      formula = f_r * f_s * f_d * p_r,
      family = "CC",
      f_s = exp(f_s_int + f_s_slope * sa_1),
      f_d = dnorm(sa_2, f_d_mu, f_d_sd),
      p_r = inv_logit(p_r_int, p_r_slope, sa_1),
      f_r = n_new / sum(f_s),
      data_list = all_param_list,
      states = list(c("sa")),
      evict_cor = TRUE,
      evict_fun = truncated_distributions("norm", "f_d")
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "F"),
        int_rule = rep('midpoint', 2),
        state_start = rep('sa', 2),
        state_end = rep('sa', 2) 
      ) 
    ) %>%
    define_domains(
      sa = c(L, U, i) # iterate over different meshpoint numbers
    ) %>%
    define_pop_state(
      n_sa = rep(1/i, i)
    ) %>%
    make_ipm(usr_funs = list(inv_logit = s_z))
  
  lambdas$lambdas[it] <- lambda(carp_ipmr_test)
  
  it <- it + 1
  
}

plot(lambdas ~ mesh_p, data = lambdas)

diff(lambdas$lambdas)

# increasing from 100 to 150 does not really make a difference for lambda, so
# i'm going to keep it at 100 meshpoints

