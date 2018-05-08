# Copy shape files of carpobrotus to appropriate folder

library(fs)

# Update shape file directory. Not sure why there
# is no overwrite argument in dir_copy
if(dir_exists(paste0(getwd(), '/Polygons'))){
  
  dir_delete(paste0(getwd(), '/Polygons'))

}

dir_copy("I:/sie/102_data_SL/PhD_Processed_Data/Polygons",
         new_path = getwd())


