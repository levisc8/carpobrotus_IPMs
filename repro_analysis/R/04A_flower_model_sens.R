# Fit the vital rate models using brms

all_ramets_temp <- readRDS("repro_analysis/Data/demography/all_ramets_clim.rds")

# Get all occurrence data so we can properly scale/center our sites. This centers
# and scales based on the complete range, rather than just the range of our field
# sites. 

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

# The purpose of this script is to see if the interaction between nativity * clim
# is still significant when invaded climate is restricted to roughly the same range
# as native climate

center_vars <- gsub("\\.y", "", center_vars)

natives <- filter(all_ramets, native == 1)

ranges <- vapply(natives[ , center_vars], range, numeric(2L))

ranges

# the natives have the warmest t_co_qu, middle of the road mat, and low precip.
# Taking +/- 0.2 on those values for subsetting data

min_t_co <- ranges[1, 1] - 0.2
max_t_co <- ranges[2, 1] + 0.2
min_mat  <- ranges[1, 2] - 0.2
max_mat  <- ranges[2, 2] + 0.2
min_map  <- ranges[1, 3] - 0.2
max_map  <- ranges[2, 3] + 0.2

repro_t_co <- filter(all_ramets,
                        native == 1 |
                          (t_co_qu_rec >= min_t_co & t_co_qu_rec <= max_t_co))

use_t_co_rams <- repro_t_co%>%
  select(flower_n, log_size, t_co_qu_rec, native) %>%
  filter(!is.na(flower_n)) %>%
  lapply(function(x) {
    attributes(x) <- NULL
    return(x)
  }) %>%
  setNames(
    c("Y", "size_l", "clim_val", "native")
  )

t_co_dat <- list(Y = use_t_co_rams$Y)

t_co_dat$N <- length(use_t_co_rams$Y)
t_co_dat$X <- model.matrix(Y ~ size_l + clim_val * native, data = use_t_co_rams)
t_co_dat$K <- dim(t_co_dat$X)[2]


repro_map <- filter(all_ramets,
                       native == 1 |
                         (map_rec >= min_map & map_rec <= max_map))


use_map_rams <- repro_map %>%
  select(flower_n, log_size, map_rec, native) %>%
  filter(!is.na(flower_n)) %>%
  lapply(function(x) {
    attributes(x) <- NULL
    return(x)
  }) %>%
  setNames(
    c("Y", "size_l", "clim_val", "native")
  )

map_dat <- list(Y = use_map_rams$Y)

map_dat$N <- length(use_map_rams$Y)
map_dat$X <- model.matrix(Y ~ size_l + clim_val * native, data = use_map_rams)
map_dat$K <- dim(map_dat$X)[2]


repro_mat <- filter(all_ramets,
                       native == 1 |
                         (mat_rec >= min_mat & mat_rec <= max_mat))

use_mat_rams <- repro_mat %>%
  select(flower_n, log_size, mat_rec, native) %>%
  filter(!is.na(flower_n)) %>%
  lapply(function(x) {
    attributes(x) <- NULL
    return(x)
  }) %>%
  setNames(
    c("Y", "size_l", "clim_val", "native")
  )

mat_dat <- list(Y = use_mat_rams$Y)

mat_dat$N <- length(use_mat_rams$Y)
mat_dat$X <- model.matrix(Y ~ size_l + clim_val * native, data = use_mat_rams)
mat_dat$K <- dim(mat_dat$X)[2]


trunc_fl_mod <- stan_model("repro_analysis/Stan/flower_n_0_trunc_mat.stan",
                           model_name = "nb_0_trunc_mat")

flower_n_t_co <- sampling(
  trunc_fl_mod,
  data = t_co_dat,
  chains = 4,
  cores  = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98,
  #                max_treedepth = 15)
)

flower_n_mat <- sampling(
  trunc_fl_mod,
  data = mat_dat,
  chains = 4,
  cores  = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98,
  #                max_treedepth = 15)
)

flower_n_map <- sampling(
  trunc_fl_mod,
  data = map_dat,
  chains = 4,
  cores  = 4L,
  iter = 4000#,
  # control = list(adapt_delta = 0.98,
  #                max_treedepth = 15)
)

print(flower_n_t_co, pars = c("b", 'shape'))


# 
# # Save everything -------------
# 
# pdf("repro_analysis/Manuscript/Figures/repro_mod_diagnostics.pdf")
# 
# plot(repro_mod_t_co,
#      ask = FALSE)
# 
# 
# pp_check(repro_mod_t_co,
#          type     = 'bars_grouped',
#          group    = 'native',
#          nsamples = 100L,
#          freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_t_co,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# plot(repro_mod_map,
#      ask = FALSE)
# 
# p <- pp_check(repro_mod_map,
#               type     = 'bars_grouped',
#               group    = 'native',
#               nsamples = 100L,
#               freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_map,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# 
# plot(repro_mod_mat,
#      ask = FALSE)
# 
# p <- pp_check(repro_mod_mat,
#               type     = 'bars_grouped',
#               group    = 'native',
#               nsamples = 100L,
#               freq     = FALSE)
# 
# print(p)
# 
# p <- pp_check(repro_mod_mat,
#               type     = 'dens_overlay',
#               nsamples = 100L)
# 
# print(p)
# 
# dev.off()
# 
# sink(file   = "repro_analysis/Manuscript/Figures/repro_mod_diagnostics.txt",
#      append = FALSE)
# 
# cat("\n\n****************** Temp Coldest Quarter ************\n\n")
# print(summary(repro_mod_t_co))
# cat("\n\n************** End Temp Coldest Quarter ************\n\n")
# 
# cat("\n\n****************** Mean Precip   *******************\n\n")
# print(summary(repro_mod_map))
# cat("\n\n************** End Mean Precip *********************\n\n")
# 
# 
# cat("\n\n****************** Mean Temp ***********************\n\n")
# print(summary(repro_mod_mat))
# cat("\n\n************** End Mean Temp ***********************\n\n")
# 
# 
# sink()
# 
# pdf("repro_analysis/Manuscript/Figures/flower_mod_diagnostics.pdf")
# 
# plot(flower_mod_t_co,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_t_co,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# plot(flower_mod_map,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_map,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# 
# plot(flower_mod_mat,
#      ask = FALSE)
# 
# p <- pp_check(flower_mod_mat,
#               type     = 'scatter_avg_grouped',
#               group    = 'native',
#               nsamples = 100L)
# 
# print(p + geom_abline(slope = 1, intercept = 0))
# 
# dev.off()
# 
# sink(file   = "repro_analysis/Manuscript/Figures/flower_mod_diagnostics.txt",
#      append = FALSE)
# 
# cat("\n\n****************** Temp Coldest Quarter ************\n\n")
# print(summary(flower_mod_t_co))
# cat("\n\n************** End Temp Coldest Quarter ************\n\n")
# 
# cat("\n\n****************** Mean Precip   *******************\n\n")
# print(summary(flower_mod_map))
# cat("\n\n************** End Mean Precip *********************\n\n")
# 
# 
# cat("\n\n****************** Mean Temp ***********************\n\n")
# print(summary(flower_mod_mat))
# cat("\n\n************** End Mean Temp ***********************\n\n")
# 
# sink()
# 
# vr_mod_list <- list(repro_mod_t_co  = repro_mod_t_co,
#                     repro_mod_map   = repro_mod_map,
#                     repro_mod_mat   = repro_mod_mat,
#                     flower_mod_t_co = flower_mod_t_co,
#                     flower_mod_map  = flower_mod_map,
#                     flower_mod_mat  = flower_mod_mat)
# 
# saveRDS(vr_mod_list, file = "repro_analysis/model_fits/vr_mod_list.rds")