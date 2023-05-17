library(terra)
library(tidyverse)
library(magrittr)
library(httr)
library(raster)

password = read_table('/home/zhoylman/private/gesdisc_creds.txt')
base_url = "https://hydro1.gesdisc.eosdis.nasa.gov/data/NLDAS/NLDAS_MOS0125_H.2.0/"
out_dir = "/mnt/data1/soil-moisture-models/nldas2/mosaic"

#or download using file zilla

dates_of_interest = seq(as.Date('1991-01-01'), as.Date('2022-12-31'), by = 'day')

for(i in 1 : length(dates_of_interest)){
  print(dates_of_interest[i])
  full_url = paste0(base_url, 
                    dates_of_interest[i] %>% lubridate::year(),
                    '/',
                    dates_of_interest[i] %>% lubridate::yday() %>% str_pad(., 3, pad = "0"),
                    '/',
                    'NLDAS_MOS0125_H.A',
                    dates_of_interest[i] %>% format(., format = '%Y%m%d'),
                    '.0000.020.nc')
  
  command_start = 'wget -q --user zhoylman --password '  
  
  command_mid = ' --content-disposition -P '
  
  system(paste0(command_start, password$pswrd, command_mid, out_dir, ' ', full_url))
}

#parallel process and write out individual tifs
files = list.files('/mnt/data1/soil-moisture-models/nldas2/mosaic/', full.names = T)
files_short = list.files('/mnt/data1/soil-moisture-models/nldas2/mosaic/', full.names = F)

for(i in 1:length(files)){
  print(files_short[i])
  temp_terra = rast(files[i])
  temp_root = temp_terra$SoilM_0_100cm
  writeRaster(temp_root, paste0('/mnt/data1/soil-moisture-models/nldas2/mosaic_0-100cm_sm/', files_short[i], '_0-100cm_sm.tif'))
}

#write out nc in chunks terra cant write out nc's in +1000 layers it doesnt seam (wierd errors)
num_files = list.files('/mnt/data1/soil-moisture-models/nldas2/mosaic_0-100cm_sm') %>% length()

nc_breaks = c(seq(1,num_files, by = 300), num_files)

for(i in 2:length(nc_breaks)){
  print((i/length(nc_breaks))*100)
  #write out CDF by groups
  group_index = seq(nc_breaks[i-1], nc_breaks[i]-1, by = 1)
  if(i == length(nc_breaks)){
    group_index = seq(nc_breaks[i-1], nc_breaks[i], by = 1)
  }
  #import and turn into raster stack
  nldas2_rasterStack = list.files('/mnt/data1/soil-moisture-models/nldas2/mosaic_0-100cm_sm', full.names = T)[group_index] %>%
    purrr::map(., raster) %>%
    stack() %>%
    rast()
  
  time = nldas2_rasterStack %>%
    sources() %>%
    substr(.,75,82) %>%
    as.Date(., format = '%Y%m%d') %>%
    as_tibble() %>%
    mutate(date = value) %>%
    dplyr::select(date)
  
  write_csv(time, paste0('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_time/mosaic_time_group_', i-1 , '.csv'))
  writeCDF(nldas2_rasterStack, paste0('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_nc/nldas2_mosaic_soil_moisture_0-100cm_group_', i-1 ,'.nc'), overwrite = T)
}

#finally, merge datasets with cdo and write time
nc_files = tibble(nc_path = list.files('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_nc', pattern = '.nc', full.names = T),
                  nc_name = list.files('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_nc', pattern = '.nc', full.names = F),
                  time_path = list.files('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_time', pattern = '.csv', full.names = T),
                  time_name = list.files('/mnt/data1/soil-moisture-models/nldas2/temp/mosaic_time', pattern = '.csv', full.names = T)) %>%
  mutate(id = substr(nc_name, 42, 999) %>% parse_number()) %>%
  arrange(id)

# use CDO to merge Climate Data Operator (CDO) command line tool
# https://www.isimip.org/protocol/preparing-simulation-files/cdo-help/
command = paste0('cdo mergetime', paste0(' ', nc_files$nc_path, collapse = ''),
                 paste0(' /mnt/data1/soil-moisture-models/nc/NLDAS2_MOSAIC_soil_moisture_0-100cm.nc'),
                 collapse = '')

system(command)

time = purrr::map(nc_files$time_path, read_csv) %>%
  bind_rows()

write_csv(time, '/mnt/data1/soil-moisture-models/nc/NLDAS2_MOSAIC_soil_moisture_0-100cm_time.csv')

