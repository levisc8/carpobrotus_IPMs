# Complete range map for "Carpobrotus.

# Data citation (downloaded from GBIF manually):
# GBIF.org (18 September 2020) GBIF Occurrence Download https://doi.org/10.15468/dl.q8kmsy
# Unless GBIF discovers citations of this download, the data file is eligible
# for deletion after March 18, 2021.
# 
# Direct download URL: https://www.gbif.org/occurrence/download/0064105-200613084148143

full_archive <- dwca_read("repro_analysis/Data/gbif_occ_all_raw.zip",
                          read = TRUE)

spp_nm_regex <- "edulis|acinaciformis|chilensis"

all_occ      <- full_archive$data$occurrence.txt %>%
  .[grepl(pattern = spp_nm_regex, x = .$specificEpithet), ]


# Need to work out quality filters