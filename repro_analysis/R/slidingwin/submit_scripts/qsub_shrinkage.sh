qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "t_max" "quad" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv

qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "t_min" "quad" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv

qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "prec" "quad" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv

qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "t_max" "lin" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv

qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "t_min" "lin" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv

qsub /home/$USER/iceplant_shrinkage/submit_shrinkage.sh "prec" "lin" /data/levin/iceplant_shrinkage/all_ramets_clim.rds /data/levin/iceplant_shrinkage/all_gbif_field_sites_plus_monthly_weather.csv
