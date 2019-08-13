

# Now we are ready to explore the data a bit! first, fit some
# growth models, then plot the data to see how it looks. The data are not
# completely digitized, so the survival model is really a place holder until
# that process is completed. 

grow_mod_int  <- lm(log_size_next ~ 1,
                    data = all_data)

grow_mod_lin  <- lm(log_size_next ~ log_size, 
                    data = all_data)

grow_mod_quad <- lm(log_size_next ~ poly(log_size, 2), 
                    data = all_data[!is.na(all_data$log_size), ])

# Now, fit a gam with a few knots to see if the linear is doing well
grow_gam <- gam(log_size_next ~ s(log_size, k = 8), data = all_data)

xx <- seq(-5, 5, 0.1)
plot(log_size_next ~ log_size, data = all_data)

abline(grow_mod_lin, col = 'red')
abline(grow_mod_int, col = 'blue')
lines(xx, predict(grow_gam, 
                  newdata = data.frame(log_size = xx),
                  type = 'response'),
      lty = 2,
      col = 'red')
lines(xx, predict(grow_mod_quad,
                  data.frame(log_size = xx),
                  type = 'response'),
      col = 'blue',
      lty = 2)

abline(a = 0, b = 1)

legend('topleft',
       legend = c('Linear Model',
                  'Intercept only',
                  'GAM', 
                  'Quadratic fit',
                  '1:1 Line'),
       col = c('red', 'blue', 'red', 'blue', 'black'),
       lty = c(1, 1, 2, 2, 1))

print(summary(grow_mod_int)) 
print(summary(grow_mod_lin))
print(summary(grow_mod_quad))
print(summary(grow_gam))

par(mfrow = c(2,2))
plot(grow_mod_lin, ask = FALSE)

grow_aic <- AIC(grow_mod_int, grow_mod_lin, grow_mod_quad, grow_gam)

# Based on the residual plots, it looks like we have a slightly decreasing
# variance with increasing size. Next, I'll try fitting a gls with non-constant
# variance and see how that improves the fit.

grow_exp_var <- gls(log_size_next ~ log_size,
                    data = all_data,
                    weights = varExp(),
                    na.action = na.omit,
                    method = 'ML')


grow_aic <- AIC(grow_mod_int, grow_mod_lin,
                grow_mod_quad, grow_gam,
                grow_exp_var)

par(mfrow = c(1, 1))

plot(log_size_next ~ log_size, data = all_data)

abline(grow_mod_lin, col = 'red')

lines(xx, 
      predict(grow_exp_var,
              data.frame(log_size = xx),
              type = 'response'),
      col = 'blue',
      lty = 2)

abline(a = 0, b = 1)

legend('topleft',
       legend = c('Linear Model - constant variance',
                  'Linear Model - non constant variance',
                  '1:1 Line'),
       col = c('red', 'blue', 'black'),
       lty = c(1, 2, 1))

# Survival models

surv_mod_int  <- glm(survival ~ 1,
                     data   = all_data,
                     family = binomial())
surv_mod_lin  <- glm(survival ~ log_size, 
                     data   = all_data, 
                     family = binomial())
surv_mod_quad <- glm(survival ~ log_size + I(log_size^2), 
                     data   = all_data, 
                     family = binomial())

surv_gam      <- gam(survival ~ s(log_size, k = 6), 
                     data   = all_data[!is.na(all_data$log_size), ],
                     family = binomial())

par(mfrow = c(1,1))
plot(survival ~ log_size, data = all_data)

lines(xx, 
      rep(1/(1 + exp(-(coef(surv_mod_int)[1]))), length(xx)), 
      col = 'blue')

lines(xx, predict(surv_mod_lin, 
                  data.frame(log_size = xx),
                  type = 'response'),
      col = 'red')

lines(xx,
      predict(surv_mod_quad, 
              data.frame(log_size = xx), 
              type = 'response'),
      col = 'blue', 
      lty = 2)
lines(xx,
      predict(surv_gam,
              data.frame(log_size = xx),
              type = 'response'),
      col = 'red',
      lty = 2)

legend('bottomright',
       legend = c('Linear Model',
                  'Intercept only',
                  'GAM', 
                  'Quadratic fit'),
       col = c('red', 'blue', 'red', 'blue'),
       lty = c(1, 1, 2, 2))

# par(mfrow = c(2,2))
# plot(surv_mod, ask = FALSE)

print(summary(surv_mod_int)) 
print(summary(surv_mod_lin))
print(summary(surv_mod_quad))
print(summary(surv_gam))

surv_aic <- AIC(surv_mod_int,
                surv_mod_lin,
                surv_mod_quad,
                surv_gam)
message('Survival AIC scores\n\n')
print(surv_aic)
message('\n\nGrowth AIC scores\n\n')
print(grow_aic)

# Extract coefficients so we can generate a "data_list" for ipmr
grow_const_var_coef_list <- c(coef(grow_mod_lin),
                              sd(resid(grow_mod_lin))) %>%
   as.list() %>% 
   setNames(c('g_int', 'g_slope', 'sd_g'))

grow_exp_var_coef_list <- c(coef(grow_exp_var),
                            as.numeric(grow_exp_var$modelStruct$varStruct)) %>%
   as.list() %>%
   setNames(c('g_int', 'g_slope', 'g_sigma_par'))


surv_coef_list <- coef(surv_mod_lin) %>%
   as.list() %>%
   setNames(c('s_int', 's_slope'))

# Fecundity parameters ----------

# Recruitment is defined as the rate of new plants(T+1) per flower at T.
# Because of the time lagged recruitment due to seed maturation time (~1.25 years
# between fertilization to seed maturation), the core assumption here is that 
# flower production is time invariant (e.g. because recruits we see actually
# come from seeds that must be at least as old as 2017).

n_flow     <- sum(all_data$flower_n, na.rm = TRUE)
new_plants <- all_data %>%
   filter(id > 8000 & is.na(log_size))

n_new <- dim(new_plants)[1]

f_r     <- n_new / n_flow

# Probability of reproducing. f_r gets multiplied by s_z and p_r
# to generate a fecundity kernel
p_r_mod_int <- glm(repro ~ 1,        data = all_data, family = binomial())
p_r_mod_lin <- glm(repro ~ log_size, data = all_data, family = binomial())

p_r_aic <- AIC(p_r_mod_int, p_r_mod_lin)

print(summary(p_r_mod_int))
print(summary(p_r_mod_lin))

message('\n\nPr(reproduction) AIC table\n\n')
print(p_r_aic)

par(mfrow = c(1, 1))
plot(repro ~ log_size, data = all_data)
lines(xx, 
      predict(p_r_mod_lin, 
              data.frame(log_size = xx),
              type = 'response'),
      col = 'red')

lines(xx, 
      rep(1/(1 + exp(-(coef(p_r_mod_int)[1]))), length(xx)), 
      col = 'blue')

# New plant size distribution --------

f_d_mu <- mean(new_plants$log_size_next, na.rm = TRUE)
f_d_sd <- sd(new_plants$log_size_next, na.rm = TRUE)

# Flower production ~ size 

f_s_mod <- glm(flower_n ~ log_size, data = all_data, family = quasipoisson())

print(summary(f_s_mod))

fec_coef_list <- c(coef(p_r_mod_lin),
                   coef(f_s_mod),
                   f_r,
                   f_d_mu,
                   f_d_sd) %>%
   as.list() %>%
   setNames(c('p_r_int', 'p_r_slope',
              'f_s_int', 'f_s_slope',
              'f_r', 'f_d_mu', 'f_d_sd'))

# End parameter estimation!