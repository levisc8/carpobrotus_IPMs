# Copy shape files of carpobrotus to appropriate folder

# Update shape file directory. Not sure why there
# is no overwrite argument in dir_copy

if(dir_exists(glue('{getwd()}/Polygons'))){
  
  dir_delete(glue('{getwd()}/Polygons'))

}

dir_copy("I:/sie/102_data_SL/PhD_Processed_Data/Polygons",
         new_path = getwd())


# Download latest version of field site overview - only if current data set is not to date

local_mod_time <- file_info('Data/All_Field_Sites.csv')$modification_time
i_mod_time <- file_info('I:/sie/102_data_SL/PhD_Processed_Data/Csv_Data/All/All_Field_Sites.csv')$modification_time

if(local_mod_time < i_mod_time) {
  file_copy('I:/sie/102_data_SL/PhD_Processed_Data/Csv_Data/All/All_Field_Sites.csv',
            new_path = 'Data/',
            overwrite = TRUE)
}

