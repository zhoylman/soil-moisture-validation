library(terra)
library(tidyverse)
library(doParallel)
library(magrittr)
library(sf)

tmp_folder = '/mnt/data1/soil-moisture-models/SPoRT/temp'
processed_folder = '/mnt/data1/soil-moisture-models/SPoRT/processed_1-100cm_mean'
raw_files = list.files('/mnt/data1/soil-moisture-models/SPoRT/raw', full.names = T)

cl = makeCluster(30)
registerDoParallel(cl)

tictoc::tic()

####
out = foreach(i = 1:length(raw_files), .packages = c('terra', 'tidyverse')) %dopar% {
  tmp_folder_i = paste0(tmp_folder, '/', i)
  
  temp_file = raw_files[i]
  
  untar(temp_file,exdir=tmp_folder_i)
  
  internal_files_temp = list.files(tmp_folder_i) %>%
    paste0(tmp_folder_i, '/', .) %>%
    tibble(full_path = list.files(., full.names = T)) %>%
    mutate(short_name = list.files(tmp_folder_i) %>%
             paste0(tmp_folder_i, '/', .) %>%
             list.files(., full.names = F)) %>%
    dplyr::select(-.) %>%
    mutate(time = substr(short_name, 10, 17))
  
  #for file
  for(f in 1:length(internal_files_temp$time)){
    temp_rast = rast(internal_files_temp$full_path[f])
    
    temp_sm_1 = temp_rast$`0-10[cm] DBLY (layer between 2 depths below land surface); Soil moisture content [kg/m^2]`
    temp_sm_2 = temp_rast$`10-40[cm] DBLY (layer between 2 depths below land surface); Soil moisture content [kg/m^2]`
    temp_sm_3 = temp_rast$`40-100[cm] DBLY (layer between 2 depths below land surface); Soil moisture content [kg/m^2]`
    #temp_sm_4 = temp_rast$`100-200[cm] DBLY (layer between 2 depths below land surface); Soil moisture content [kg/m^2]`
    
    temp_ave = list(temp_sm_1, temp_sm_2, temp_sm_3) %>%
      terra::rast() %>%
      mean()
    
    writeRaster(temp_ave, paste0(processed_folder, '/SPoRT_mean_sm_0-100cm_', internal_files_temp$time[f], '.tif'), overwrite=TRUE) 
  }
  
  unlink(tmp_folder_i, recursive = TRUE)
  gc()
}

stopCluster(cl)

tictoc::toc()


#######################################################
#######################################################
#        clip, merge and export as nc by year
#######################################################
#######################################################

`%notin%` = Negate(`%in%`) 

conus = read_sf('/home/zhoylman/mco-drought-indicators/processing/base-data/raw/states.shp') %>%
  filter(STATE_ABBR %notin% c('AK', 'HI', 'VI'))

import_project_clip = function(img){
  temp = rast(img) %>% 
    project(., crs('EPSG:4326')) %>% 
    terra::mask(., conus)
  return(temp)
}

names = tibble(name = list.files('/mnt/data1/soil-moisture-models/SPoRT/processed_1-100cm_mean')) %>%
  mutate(time_char = substr(name,23, 30),
         time = as.Date(time_char, format = "%Y%m%d"),
         year = lubridate::year(time),
         groups = rep(1:100, each=100)[1:length(year)])

groups = unique(names$groups)

for( i in 1:length(groups) ){
  temp_time = names[which(names$groups == groups[i]),]
  
  # convert to nc 
  prcoessed_files = list.files('/mnt/data1/soil-moisture-models/SPoRT/processed_1-100cm_mean', full.names = T)[which(names$groups == groups[i])] %>%
    purrr::map(., import_project_clip) %>%
    rast() 
  
  names(prcoessed_files) = temp_time$time_char
  
  # prcoessed_files = prcoessed_files %T>% 
  #   terra::set.values()
  
  #names = list.files('/mnt/data1/soil-moisture-models/SPoRT/processed')
  #names(prcoessed_files) = names
  #writeRaster(prcoessed_files, paste0('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc/SPoRT_mean_soil_moisture_',years[i],'.tif'), overwrite = T)
  writeCDF(prcoessed_files, paste0('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm/SPoRT_mean_soil_moisture_0-100cm_',groups[i],'.nc'), 
           overwrite = T, varname = 'SPoRT_mean_soil_moisture_0-100cm')
  write_csv(temp_time, paste0('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm/SPoRT_mean_soil_moisture_time_0-100cm_',groups[i],'.csv'))
  rm(prcoessed_files)
  gc()
}

#finally, merge datasets with cdo and write time

nc_files = tibble(nc_path = list.files('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm', pattern = '.nc', full.names = T),
                  nc_name = list.files('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm', pattern = '.nc', full.names = F),
                  time_path = list.files('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm', pattern = '.csv', full.names = T),
                  time_name = list.files('/mnt/data1/soil-moisture-models/SPoRT/temp_chunk_nc_0-100cm', pattern = '.csv', full.names = T)) %>%
  mutate(id = substr(nc_name, 34, 35) %>% as.numeric()) %>%
  arrange(id)

# use CDO to merge Climate Data Operator (CDO) command line tool
# https://www.isimip.org/protocol/preparing-simulation-files/cdo-help/
command = paste0('cdo mergetime', paste0(' ', nc_files$nc_path, collapse = ''),
                 paste0(' /mnt/data1/soil-moisture-models/nc/SPoRT__mean_soil_moisture_0-100cm.nc'),
                 collapse = '')

system(command)

time = purrr::map(nc_files$time_path, read_csv) %>%
  bind_rows()

write_csv(time, '/mnt/data1/soil-moisture-models/nc/SPoRT__mean_soil_moisture_0-100cm.csv')
