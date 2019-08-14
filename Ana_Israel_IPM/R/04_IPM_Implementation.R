
s_z <- function(int, slope, sv1) {
  1/(1 + exp(-(int + slope * sv1)))
}

g_ev_z1z <- function(sv1, sv2, int, slope, var_coef, L, U, h) {
  
  mu   <- int + slope * sv1
  sd_g <- sqrt(exp(2 * var_coef * sv1))
  
  ev <- pnorm(U, mu, sd_g) - pnorm(L, mu, sd_g)
  
  out  <- (h * dnorm(sv2, mu, sd_g)) / ev
  
  return(out)
}

p_ev_z1z <- function(sv1, sv2, s_int, s_slope, g_int, g_slope, var_coef, L, U, h) {
  
  g <- g_ev_z1z(sv1, sv2, g_int, g_slope, var_coef, L, U, h)
  
  return(s_z(s_int, s_slope, sv1) * g)
}

f_z1z <- function(sv1, sv2, pr_int, pr_slope, f_s_int, f_s_slope,
                  fr, fd_mu, fd_sd, s_int, s_slope, L, h) {
  
  pr_est <- 1/(1 + exp(-(pr_int + pr_slope * sv1)))
  s_est <- s_z(s_int, s_slope, sv1)
  f_est <- exp(f_s_int + f_s_slope * sv1)
  d_est <- (h * dnorm(sv2, fd_mu, fd_sd)) / (1 - pnorm(L, fd_mu, fd_sd))
  out <- pr_est * s_est * fr * f_est * d_est
  return(out)
  
}

L <- min(all_data$log_size_next, na.rm = TRUE) * 0.8
U <- max(all_data$log_size_next, na.rm = TRUE) * 1.2

n_mesh_p <- 100

b <- seq(L, U, length.out = n_mesh_p + 1)
d1 <- d2 <- (b[2:(n_mesh_p + 1)] + b[1:n_mesh_p]) * 0.5
h <- d1[2] - d1[1]

P_exp_var <- outer(d1, d2,
                   FUN      = p_ev_z1z,
                   s_int    = surv_coef_list$s_int,
                   s_slope  = surv_coef_list$s_slope,
                   g_int    = grow_exp_var_coef_list$g_int,
                   g_slope  = grow_exp_var_coef_list$g_slope,
                   var_coef = grow_exp_var_coef_list$g_sigma_par,
                   L = L,
                   U = U,
                   h = h) %>%
  t()


Fm <- outer(d1, d2,
            FUN       = f_z1z,
            pr_int    = coef(p_r_mod_lin)[1],
            pr_slope  = coef(p_r_mod_lin)[2],
            f_s_int   = coef(f_s_mod)[1],
            f_s_slope = coef(f_s_mod)[2],
            fr        = f_r,
            fd_mu     = f_d_mu,
            fd_sd     = f_d_sd,
            s_int     = coef(surv_mod_lin)[1],
            s_slope   = coef(surv_mod_lin)[2],
            L = L,
            h = h) %>% 
  t()

K_exp_var   <- (P_exp_var   + Fm)

lambda_exp_var   <- Re(eigen(K_exp_var)$values[1])

message('\n\nLambda for exponential variance growth: ', round(lambda_exp_var, 3))

