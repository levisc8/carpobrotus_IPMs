# Shrinkage models test script


suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

# Print attached packages
sessionInfo()

# Package versions for all loaded namespaces (including attached packages and
# packages loaded by attached namespaces, but not actually on the search path)
for (package_name in sort(loadedNamespaces())) {
  print(paste(package_name, packageVersion(package_name)))
}

start <- Sys.time()

print(getwd())

#  ----------------------------------------------------------------------------------------------------------------------------
# parsing arguments from submit script in EVE
#  ----------------------------------------------------------------------------------------------------------------------------

Parsoptions <- list (
  
  make_option(
    opt_str = c("-c", "--climate"),
    dest    = "climate_used",
    help    = "Specify the climate variable that will be used",
    metavar = "t_max|t_min|prec"),
  make_option(
    opt_str = c("-f", "--func"),
    dest    = "pred_fun",
    help    = "Specify the functional form of the climate predictor (linear or quadratic)",
    metavar = "lin|quad"),
  make_option(
    opt_str = c("-i", "--demo"),
    dest    = "demo_fp",
    help    = "Specify the file path for the incoming demographic data",
    metavar = "",
    default = "/data/levin/iceplant_shrinkage/all_ramets_clim.rds"),
  make_option(
    opt_str = c("-z", "--occ"),
    dest    = "occ_fp",
    help    = "Specify the file path for the incoming occurrence data",
    metavar = "",
    default = "/data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv")
  
)

parser <- OptionParser(
  usage       = "Rscript %prog [options] output",
  option_list = Parsoptions,
  description = "\nRun Bayesian shrinkage analysis",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 1)

print(cli)

clim    <- tolower(cli$options$climate_used)
func    <- tolower(cli$options$pred_fun)
demo_fp <- tolower(cli$options$demo_fp)
occ_fp  <- tolower(cli$options$occ_fp)
out_dir <- cli$args[1]

all_ramets_temp <- readRDS(demo_fp)

# Get all occurrence data so we can properly scale/center our sites.

all_occ_clim <- read.csv(occ_fp,
                         stringsAsFactors = FALSE)

temp <- select(all_occ_clim,
               site, 
               t_min_Jan_2016:prec_Dec_2018) %>%
  left_join(all_ramets_temp, by = c("site" = "population")) %>%
  select(-c(id_pop, mat_rec:native))


center_ind <- grep(clim, names(temp))


for(i in center_ind) {
  
  temp[ , i] <- scale(temp[ , i],
                      center = TRUE,
                      scale = TRUE)
  
}

# Set up data sets so we can use resp ~ . notation. 

all_data <- temp %>%
  filter(!is.na(id)) %>%
  select(flower_n, repro, log_size, starts_with(clim)) 

f_n_data <- select(all_data, -repro)
repr_data <- select(all_data, -flower_n)

# if we have a quadratic model, then we want to set up the 
# clim_var + I(clim_var)^2 for each one on the RHS of the formula

if(func == "quad") {
  
  mon_base <- expand.grid(list(month.abb,
                               2016:2018))
  
  base_form <- paste(clim, mon_base[ , 1], mon_base[ , 2], sep = "_")
  
  sq_rhs    <- paste(base_form, paste("I(", base_form, ")^2"), sep = " + ") %>%
    paste(collapse = " + ")
  
  repr_form <- paste("repro ~ log_size + ", sq_rhs, sep = "")
  f_n_form  <- paste("flower_n | trunc(lb = 1) ~ log_size + ", sq_rhs, sep = "")
  
  repr_bf <- bf(as.formula(repr_form),
                family = bernoulli())
  
  f_n_bf  <- bf(as.formula(f_n_form),
                family = poisson())
  
} else {
  
  f_n_bf <- bf(flower_n | trunc(lb = 1) ~ .,
               family = poisson())
  
  repr_bf <- bf(repro ~ ., family = bernoulli())
  
}

# Pr(repro) model ------------- 

n_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
options(mc.cores = n_cores)


repro_mod<- brm(repr_bf,
                data = repr_data,
                prior = prior("horseshoe(df = 2)"),
                chains = 8,
                iter = 4000,
                warmup = 2000,
                cores = n_cores,
                save_model = paste(out_dir,
                                   '/stan_code/repro_code.stan',
                                   sep = ""),
                control = list(adapt_delta = 0.99,
                               max_treedepth = 20))


saveRDS(repro_mod, file = paste(out_dir, 
                                "/models/repro_model_",
                                clim, 
                                "_", 
                                func,
                                ".rds", sep = ""))

## Flower model ------------------

flower_n_mod <- brm(f_n_bf,
                    data = f_n_data,
                    chains = 8,
                    prior = prior("horseshoe(df = 2)"),
                    iter = 4000,
                    warmup = 2000,
                    cores = n_cores,
                    save_model = paste(out_dir,
                                       '/stan_code/f_n_code.stan',
                                       sep = ""),
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 20))

mod_list <- list(
  repro_mod    = repro_mod,
  flower_n_mod = flower_n_mod
)

out_path <- paste(out_dir, 
                  "/models/model_objects_",
                  clim, 
                  "_", func,
                  ".rds", sep = "")

saveRDS(mod_list, file = out_path)

cat("***************** REPRO MODEL*************")

summary(repro_mod)

cat("\n\n\n***************** FLOWER MODEL*************")

summary(flower_n_mod)

