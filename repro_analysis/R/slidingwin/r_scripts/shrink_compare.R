# Shrinkage model comparison

suppressPackageStartupMessages(library(brms))
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

#  --------------------------------------------------------------------------
# parsing arguments from submit script in EVE
#  --------------------------------------------------------------------------

Parsoptions <- list(
  
  make_option(
    opt_str = c("-v", "--vital"),
    dest    = "vital_rate",
    help    = "Specify the vital rate (pr(flowering) or flower #",
    metavar = "repro|flower_n"
  )
)

parser <- OptionParser(
  usage       = "Rscript %prog [options] output",
  option_list = Parsoptions,
  description = "\nCompare all models using WAIC and LOO CV",
  epilogue    = ""
)

cli <- parse_args(parser, positional_arguments = 1)

print(cli)

vr      <- tolower(cli$options$vital_rate)

out_dir <- cli$args[1]

mod_list_fp <- paste("/data/levin/iceplant_shrinkage/model_lists/all_models_",
                     vital_rate,
                     ".rds", 
                     sep = "")

mod_list <- readRDS(mod_list_fp)

n_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
options(mc.cores = n_cores)

waics <- waic(
  mod_list[[1]], # NULL size only model
  mod_list[[2]], # Best model using MAT
  mod_list[[3]], # Best model using MAP
  mod_list[[4]], # Best model using T_Co
  mod_list[[5]], # Shrinkage model lin
  mod_list[[6]], # Shrinkage model quad
  compare = TRUE
)

out_waic <- paste(out_dir, "/waic_", vr, sep = "")

saveRDS(waics, file = paste(out_waic, ".rds", sep = ""))

loos <- loo(
  mod_list[[1]], # NULL size only model
  mod_list[[2]], # Best model using MAT
  mod_list[[3]], # Best model using MAP
  mod_list[[4]], # Best model using T_Co
  mod_list[[5]], # Shrinkage model lin
  mod_list[[6]], # Shrinkage model quad
  compare    = TRUE,
  reloo      = TRUE,
  reloo_args = list(chains = 8,
                    iter   = 4000,
                    warmup = 2000,
                    cores  = n_cores)
)


cat("************WAIC RESULTS***************\n")
print(waics)

cat("\n************LOO RESULTS***************\n")
print(loos)

out_fp <- paste(out_dir, "/cv_results/", clim, "_", vr, sep = "")

saveRDS(list(waic = waics,
             loos = loos),
        file = paste(out_fp, ".rds", sep = ""))

cat("\n\nTotal time: ", Sys.time() - start)