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

#  -----------------------------------------------------------------------------
# parsing arguments from submit script in EVE
#  -----------------------------------------------------------------------------

Parsoptions <- list (
  
  make_option(
    opt_str = c("-c", "--climate"),
    dest    = "climate_used",
    help    = "Specify the climate variable that will be used",
    metavar = "tmax|tmin|prec"),
  make_option(
    opt_str = c("-f", "--func"),
    dest    = "pred_fun",
    help    = "Specify the functional form of the climate predictor (linear or quadratic)",
    metavar = "lin|quad"),
  make_option(
    opt_str = c("-v", "--vital"),
    dest    = "vr",
    help    = "Specify the vital rate to model (pr(flowering) or # flowers)",
    metavar = "repro|flower_n"),
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
vr      <- tolower(cli$options$vr)
out_dir <- cli$args[1]

all_ramets_temp <- readRDS(demo_fp)

# Get all occurrence data so we can properly scale/center our sites.

all_occ_clim <- read.csv(occ_fp,
                         stringsAsFactors = FALSE) 

all_occ_clim <- all_occ_clim[complete.cases(all_occ_clim), ]

temp <- select(all_occ_clim,
               site, 
               prec_t_minus_1:tmin_t_minus_40) %>%
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

if(vr == "flower_n"){
  
  use_data <- select(all_data, -repro)
  
  if(func == "quad") {
    
    
    base_form <- paste(clim, "_t_minus_", 1:40, sep = "")
    
    sq_rhs    <- paste(base_form,
                       paste("I(", base_form, "^2)", sep = ""),
                       sep = " + ") %>%
      paste(collapse = " + ")
    
    f_n_form  <- paste("flower_n | trunc(lb = 1) ~ log_size + ", sq_rhs, sep = "")
  
    use_bf  <- bf(as.formula(f_n_form),
                  family = poisson())
    
  } else {
    
    use_bf <- bf(flower_n | trunc(lb = 1) ~ .,
                 family = poisson())
  }

} else {
  
  use_data <- select(all_data, -flower_n)
  
  if(func == "quad") {
    
    
    base_form <- paste(clim, "_t_minus_", 1:40, sep = "")
    
    sq_rhs    <- paste(base_form,
                       paste("I(", base_form, "^2)", sep = ""),
                       sep = " + ") %>%
      paste(collapse = " + ")
    
    repr_form <- paste("repro ~ log_size + ", sq_rhs, sep = "")

    use_bf <- bf(as.formula(repr_form),
                  family = bernoulli())
    
  } else {
    
    
    use_bf <- bf(repro ~ ., family = bernoulli())
    
  }
  
  
  
}
# if we have a quadratic model, then we want to set up the 
# clim_var + I(clim_var^2) for each one on the RHS of the formula


# Pr(repro) model ------------- 

n_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
options(mc.cores = n_cores)


vr_mod <- brm(use_bf,
              data = use_data,
              prior = prior("horseshoe(df = 2)"),
              chains = 8,
              iter = 4000,
              warmup = 2000,
              cores = n_cores,
              save_model = paste(out_dir,
                                 '/stan_code/',
                                 vr, '
                                   _code.stan',
                                 sep = ""),
              control = list(adapt_delta = 0.99,
                             max_treedepth = 20))


saveRDS(vr_mod, file = paste(out_dir, 
                             "/models/vr_model_",
                             clim,
                             "_",
                             vr,
                             "_", 
                             func,
                             ".rds", sep = ""))

cat("***************** Results summary*************")

summary(vr_mod)


cat("\n\ntotal time: ", Sys.time() - start)