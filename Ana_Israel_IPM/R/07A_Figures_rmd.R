# Final plots

# First, get the outputs into a reasonable format

exp_var_temp <- lapply(exp_var_out,
                       function(x) {
                         temp <- sort(x[2:1001])
                         return(c(x[1], temp[25], temp[975]))
                       }) %>%
  as_tibble() %>%
  mutate(boot_obs = c('Observed', 'Lower_CI', 'Upper_CI'))

