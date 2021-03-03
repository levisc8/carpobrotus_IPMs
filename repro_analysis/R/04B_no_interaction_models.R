# Same as other models, but with centered variables

# Fit the vital rate models using brms

all_ramets_temp <- readRDS("repro_analysis/Data/demography/all_ramets_clim.rds")

# Get all occurrence data so we can properly scale/center our sites.

all_occ_clim <- read.csv("repro_analysis/Data/all_gbif_field_sites.csv",
                         stringsAsFactors = FALSE)

temp <- select(all_occ_clim,
               site, 
               t_co_qu_rec, mat_rec, map_rec) %>%
  left_join(all_ramets_temp, by = c("site" = "population"))

center_vars <- c("t_co_qu_rec.y", "mat_rec.y", "map_rec.y")

for(i in seq_along(center_vars)) {
  
  temp[ , center_vars[i]] <- scale(temp[ , center_vars[i]], 
                                   center = TRUE,
                                   scale = TRUE)
  
}

all_ramets <- temp %>%
  filter(!is.na(id)) %>%
  select(-c(t_co_qu_rec.x:map_rec.x)) %>%
  setNames(
    names(all_ramets_temp)
  )

# Next, we want to fit the climate ~ native interaction models. Keeping the 
# climate vars separate for now, because I'm not sure how ugly that interaction
# would get. Start w/ temp_coldest_quarter, mean annual temp, and mean annual
# precip. DO NOT combine MAT and T_CO_QU! T_seas and P_seas may be needed later,
# depending on results.


# Temp Coldest Quarter ------------- 

repro_mod_t_co <- brm(repro ~ log_size + t_co_qu_rec,
                      data = all_ramets,
                      family = bernoulli(),
                      chains = 4,
                      iter = 4000,
                      warmup = 2000,
                      cores = getOption("mc.cores", 4L),
                      save_model = 'repro_analysis/Stan/no_nat_clim_pr_flower_mod_t_co.stan',
                      control = list(adapt_delta = 0.99,
                                     max_treedepth = 15))

## Mean Annual Temperature ------------------

repro_mod_mat <- brm(repro ~ log_size + mat_rec,
                     data = all_ramets,
                     family = bernoulli(),
                     chains = 4,
                     iter = 4000,
                     warmup = 2000,
                     cores = getOption("mc.cores", 4L),
                     save_model = 'repro_analysis/Stan/no_nat_clim_pr_flower_mod_mat.stan',
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 15))


# Mean annual precip ----------------

repro_mod_map <- brm(repro ~ log_size + map_rec,
                     data = all_ramets,
                     family = bernoulli(),
                     chains = 4,
                     iter = 4000,
                     warmup = 2000,
                     cores = getOption("mc.cores", 4L),
                     save_model = 'repro_analysis/Stan/no_nat_clim_pr_flower_mod_map.stan',
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 15))




# Flower models! This uses a 0-truncated negative binomial, so had to adapt
# some stan code, as stan crashes when using the Y|trunc(lb = 0) ... syntax
# in brms. 

fl_dat <- filter(all_ramets, !is.na(flower_n))

t_co_dat <- list(
  Y = fl_dat$flower_n,
  N = length(fl_dat$flower_n),
  X = model.matrix(flower_n ~ log_size + t_co_qu_rec,
                   data = fl_dat)
)

t_co_dat$K <- dim(t_co_dat$X)[2]

mat_dat <- list(
  Y = fl_dat$flower_n,
  N = length(fl_dat$flower_n),
  X = model.matrix(flower_n ~ log_size + mat_rec,
                   data = fl_dat)
)

mat_dat$K <- dim(mat_dat$X)[2]

map_dat <- list(
  Y = fl_dat$flower_n,
  N = length(fl_dat$flower_n),
  X = model.matrix(flower_n ~ log_size + map_rec,
                   data = fl_dat)
)

map_dat$K <- dim(map_dat$X)[2]

fl_mod <- stan_model("repro_analysis/Stan/flower_n_0_trunc_mat.stan")

fl_mod_t_co <- sampling(
  fl_mod, 
  data = t_co_dat,
  chains = 4,
  cores = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98, max_treedepth = 15)
)


fl_mod_mat <- sampling(
  fl_mod, 
  data = mat_dat,
  chains = 4,
  cores = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98, max_treedepth = 15)
)


fl_mod_map <- sampling(
  fl_mod, 
  data = map_dat,
  chains = 4,
  cores = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98, max_treedepth = 15)
)



# Save everything -------------

pdf("repro_analysis/Manuscript/Figures/repro_no_nat_mod_diagnostics_centered.pdf")

plot(repro_mod_t_co,
     ask = FALSE)


p <- pp_check(repro_mod_t_co,
              type     = 'bars',
              nsamples = 100L,
              freq     = FALSE)

print(p)

p <- pp_check(repro_mod_t_co,
              type     = 'dens_overlay',
              nsamples = 100L)

print(p)

plot(repro_mod_map,
     ask = FALSE)

p <- pp_check(repro_mod_map,
              type     = 'bars',
              nsamples = 100L,
              freq     = FALSE)

print(p)

p <- pp_check(repro_mod_map,
              type     = 'dens_overlay',
              nsamples = 100L)

print(p)


plot(repro_mod_mat,
     ask = FALSE)

p <- pp_check(repro_mod_mat,
              type     = 'bars',
              nsamples = 100L,
              freq     = FALSE)

print(p)

p <- pp_check(repro_mod_mat,
              type     = 'dens_overlay',
              nsamples = 100L)

print(p)

dev.off()

sink(file   = "repro_analysis/model_summaries/repro_no_nat_mod_diagnostics_centered.txt",
     append = FALSE)

cat("\n\n****************** Temp Coldest Quarter ************\n\n")
print(summary(repro_mod_t_co))

cat("\n\n************** End Temp Coldest Quarter ************\n\n")

cat("\n\n****************** Mean Precip   *******************\n\n")
print(summary(repro_mod_map))

cat("\n\n************** End Mean Precip *********************\n\n")


cat("\n\n****************** Mean Temp ***********************\n\n")
print(summary(repro_mod_mat))
cat("\n\n************** End Mean Temp ***********************\n\n")


sink()

par_nms <- c("b_Intercept",
             paste("b[", 1:2, "]", sep = ""),
             "shape")


make_y_reps <- function(model, n_reps) {
  
  pars <- rstan::extract(model, pars = c("mu", "shape"))
  
  n_obs <- dim(pars$mu)[2]
  
  y_reps <- matrix(0, nrow = n_reps, ncol = n_obs)
  
  for(i in seq_len(n_reps)) {
    
    mus <- exp(pars$mu[i, ])
    
    y_reps[i, ] <- rztnbinom(n_obs, mus, pars$shape[i])
    
  }
  
  return(y_reps)
  
}

y_rep_t_co <- make_y_reps(fl_mod_t_co, n_reps = 500)
y_rep_mat  <- make_y_reps(fl_mod_mat, n_reps = 500)
y_rep_map  <- make_y_reps(fl_mod_map, n_reps = 500)

pdf("repro_analysis/Manuscript/Figures/flower_no_nat_mod_diagnostics_centered.pdf")

mcmc_dens(as.array(fl_mod_t_co),
          pars = par_nms)

mcmc_trace(as.array(fl_mod_t_co),
           pars = par_nms)


p1 <- ppc_dens_overlay(t_co_dat$Y,
                       y_rep_t_co) + xlim(c(0, 200))

grid.arrange(
  p1, p1 + xlim(c(0, 10))
)

ppc_bars(t_co_dat$Y,
         y_rep_t_co) + xlim(c(0, 30))


mcmc_dens(as.array(fl_mod_map),
          pars = par_nms)

mcmc_trace(as.array(fl_mod_map),
           pars = par_nms)


p1 <- ppc_dens_overlay(map_dat$Y,
                       y_rep_map) + xlim(c(0, 200))

grid.arrange(
  p1, p1 + xlim(c(0, 10))
)

ppc_bars(map_dat$Y,
         y_rep_t_co) + xlim(c(0, 30))

mcmc_dens(as.array(fl_mod_mat),
          pars = par_nms)

mcmc_trace(as.array(fl_mod_mat),
           pars = par_nms)


p1 <- ppc_dens_overlay(mat_dat$Y,
                       y_rep_mat) + xlim(c(0, 200))

grid.arrange(
  p1, p1 + xlim(c(0, 10))
)

ppc_bars(mat_dat$Y,
         y_rep_t_co) + xlim(c(0, 30))

dev.off()

sink(file   = "repro_analysis/model_summaries/flower_no_nat_mod_diagnostics_centered.txt",
     append = FALSE)

cat("\n\n****************** Temp Coldest Quarter ************\n\n")
print(fl_mod_t_co, pars = c("b_Intercept", "b", "shape"))

nms <- paste(paste("\nb[",1:2, "]: ", sep = ""),
             dimnames(t_co_dat$X)[[2]][-1],
             sep = "")
cat(nms)

cat("\n\n************** End Temp Coldest Quarter ************\n\n")

cat("\n\n****************** Mean Precip   *******************\n\n")
print(fl_mod_map, pars = c("b_Intercept", "b", "shape"))

nms <- paste(paste("\nb[",1:2, "]: ", sep = ""),
             dimnames(map_dat$X)[[2]][-1], 
             sep = "")
cat(nms)

cat("\n\n************** End Mean Precip *********************\n\n")


cat("\n\n****************** Mean Temp ***********************\n\n")
print(fl_mod_mat, pars = c("b_Intercept", "b", "shape"))
nms <- paste(paste("\nb[",1:2, "]: ", sep = ""), 
             dimnames(mat_dat$X)[[2]][-1],
             sep = "")
cat(nms)
cat("\n\n************** End Mean Temp ***********************\n\n")

sink()

vr_mod_list <- list(repro_mod_t_co  = repro_mod_t_co,
                    repro_mod_map   = repro_mod_map,
                    repro_mod_mat   = repro_mod_mat,
                    flower_mod_t_co = fl_mod_t_co,
                    flower_mod_map  = fl_mod_map,
                    flower_mod_mat  = fl_mod_mat)

saveRDS(vr_mod_list, file = "repro_analysis/model_fits/vr_no_nat_mod_list_centered.rds")
