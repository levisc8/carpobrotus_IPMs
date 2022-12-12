# Apparently these packages do not play nice w/ KrigR, so 
# shift attaching until after 01_Climate_Data_Download.R

library(stringr)   # string manipulation
library(tidyr)     # Climate data munging
library(dplyr)     # data manipulation
library(rlang)
library(fs)        # file manipulation
library(ggplot2)   # plotting
library(gridExtra) # plotting
library(patchwork)
library(viridis)
library(glue)      # string manipulation
library(purrr)
library(sf)        # Handle spatial data
library(brms)  
library(bayesplot)
library(mgcv)
library(rmarkdown)
library(ipmr)
library(data.table)
library(randomForest)
library(forcats)