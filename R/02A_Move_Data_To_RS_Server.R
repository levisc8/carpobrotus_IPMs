

i_cd <- 'I:/sie/102_data_SL/carpobrotus_ipms'

if(dir_exists(glue('{i_cd}/Data'))) {
  dir_delete(glue('{i_cd}/Data'))
  
  dir_copy('Data', glue('{i_cd}/Data'))
}
