# Complete range map for "Carpobrotus.

# Data citation (downloaded from GBIF manually):
# GBIF.org (18 September 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.q8kmsy
# Unless GBIF discovers citations of this download, the data file is eligible
# for deletion after March 18, 2021.
#
# Direct download URL: https://www.gbif.org/occurrence/download/0064105-200613084148143
# 
# full_archive <- dwca_read("repro_analysis/Data/gbif/gbif_occ_all_raw.zip",
#                           read = TRUE)
# 
# spp_nm_regex <- "edulis|acinaciformis|chilensis"
# 
# all_occ      <- full_archive$data$occurrence.txt %>%
#   .[grepl(pattern = spp_nm_regex, x = .$specificEpithet), ]
# 
# write.csv(all_occ, file = "repro_analysis/Data/gbif/gbif_occ_data.csv",
#           row.names = FALSE)
# 
# all_occ <- read.csv("repro_analysis/Data/gbif/gbif_occ_data.csv",
#                     stringsAsFactors = FALSE)
# 
# # columns with all NAs aren't going to be useful for us I don't think.
# test_nas <- vapply(all_occ, function(x) all(is.na(x)), logical(1L))
# 
# temp_occ <- all_occ[ , !test_nas]
# 
# # Quality filters
# # First criteria - remove ones that have coordinate uncertainty greater than
# # 5 km. This will likely throw off our chelsa estimates (grids are ~1km square)
# 
# temp_occ <- temp_occ %>%
#   select(species, 
#          decimalLongitude, 
#          decimalLatitude, 
#          countryCode,
#          individualCount,
#          gbifID,
#          family, taxonRank,
#          coordinateUncertaintyInMeters, year,
#          basisOfRecord, institutionCode, datasetName) %>%
#   filter(
#     coordinateUncertaintyInMeters <= 5000
#   )
# 
# test_occ <- clean_coordinates(
#   temp_occ,
#   lon = "decimalLongitude",
#   lat = "decimalLatitude",
#   tests = c(
#     "capitals",
#     "centroids",
#     "equal", 
#     "gbif", 
#     "institutions",
#     "outliers",
#     "seas",
#     "zeros"
#   )
# )
# 
# # Show results
# plot(test_occ, 
#      lat = "decimalLatitude",
#      lon = "decimalLongitude")
# 
# # create vector of weird ones we need to get rid of. Capitals and centroids are easy.
# 
# rm_ind <- unique(
#   c(
#     which(!test_occ$.cap),
#     which(!test_occ$.cen),
#     which(!test_occ$.inst)
#   )
# )
# 
# # Next, lets plot the ones flagged as ocean. A lot of those look legit (e.g.
# # Azores, Canaries, or just close to the coastline which is where it grows anyway).
# 
# test_sea <- test_occ[!test_occ$.sea, ]
# 
# plot(test_sea, 
#      lat = "decimalLatitude",
#      lon = "decimalLongitude")
# 
# # I'm going to keep these, I don't see any immediate problems. The last flag is
# # outliers.
# 
# test_out <- test_occ[!test_occ$.otl, ]
# 
# plot(test_out, 
#      lat = "decimalLatitude",
#      lon = "decimalLongitude")
# 
# # These don't look problematic at all. Not really sure what this is testing for 
# # anyway. 
# 
# clean_occ <- temp_occ[-c(rm_ind), ]
# 
# write.csv(clean_occ, 
#           file = "repro_analysis/Data/gbif/gbif_occ_clean.csv",
#           row.names = FALSE)
