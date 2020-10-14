# Main script. Sources smaller components and controls 
# the order of dependency loading

source("repro_analysis/R/01_depends.R")
source("repro_analysis/R/02_gbif_cleaning.R")
source("repro_analysis/R/03_clim_envelopes.R")
source("repro_analysis/R/04_flower_models.R")
source("repro_analysis/R/05_flower_model_visualization.R")
source("repro_analysis/R/05A_model_clim_range_sensitivity.R")
