# Shrinkage model comparison

suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(future))
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
    metavar = "repro|flower"
  ),
  make_option(
    opt_str = c("-c", "--climate"),
    dest    = "climate",
    help    = "Specify the climate to use (temperature/precipitation)",
    metavar = "prec|temp"
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

vital_rate <- tolower(cli$options$vital_rate)
climate    <- tolower(cli$options$climate)

out_dir <- cli$args[1]

mod_list_fp <- paste("/data/levin/iceplant_shrinkage/",
                     vital_rate,
                     "_",
                     climate,
                     "_mod_list.rds", 
                     sep = "")

mod_list <- readRDS(mod_list_fp)

n_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
options(mc.cores = n_cores)

waics <- waic(
  mod_list[[1]], # Shrinkage model lin  
  mod_list[[2]], # Shrinkage model quad
  compare = TRUE
)


cat("************WAIC RESULTS***************\n")
print(waics)
cat("\n***************************************")

out_waic <- paste(out_dir, "/cv_results/waic_", vital_rate, "_", climate, sep = "")

saveRDS(waics, file = paste(out_waic, ".rds", sep = ""))

plan("multicore")

if(vital_rate == "repro"){
  
  
  loos <- loo(
    mod_list[[1]], # Shrinkage model lin  
    mod_list[[2]], # Shrinkage model quad
    compare    = TRUE,
    reloo      = TRUE,
    reloo_args = list(chains = 8,
                      iter   = 4000,
                      warmup = 2000,
                      cores  = n_cores)
  )
} else {

  loos <- kfold(
    mod_list[[1]], # Shrinkage model lin
    mod_list[[2]], # Shrinkage model quad
    chains  = 8, 
    iter    = 4000,
    warmup  = 2000,
    cores   = n_cores,
    K       = 10,
    compare = TRUE
  )
  
}

plan("sequential")


cat("\n************LOO RESULTS***************\n")
print(loos)

out_fp <- paste(out_dir, "/cv_results/",  vital_rate, "_", climate, sep = "")

saveRDS(list(waic = waics,
             loos = loos),
        file = paste(out_fp, ".rds", sep = ""))

cat("\n\nTotal time: ", Sys.time() - start)