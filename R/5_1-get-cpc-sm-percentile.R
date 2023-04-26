library(terra)
library(tidyverse)
library(magrittr)
library(httr)

#or download using file zilla

url = "ftp://ftp.cpc.ncep.noaa.gov/GIS/USDM_Products/soil/percentile/daily/"

files = read_table(url) %>%
  dplyr::select(`w.rank.20080205.zip`) %>%
  rename(file = `w.rank.20080205.zip`) %>%
  filter(str_detect(file, '.tif')) %>%
  mutate(full_url = paste0(url, file))

#file_location_vec = paste0(url, files$full_url)

for(i in 1:length(files$file)){
  GET(files$full_url[i], write_disk(paste0('/mnt/data1/soil-moisture-models/cpc/', files$file[i])))
}

#pull in files to convert to ncdf
cpc_data = list.files('/mnt/data1/soil-moisture-models/cpc/', full.names = T) %>%
  purrr::map(., rast) %>%
  rast()

cpc_time = cpc_data %>% 
  names() %>%
  gsub("\\D", "", .) %>%
  as_tibble() %>%
  mutate(date = as.Date(value, format = '%Y%m%d'))

write_csv(cpc_time, file = '/mnt/data1/soil-moisture-models/nc/cpc_soil_moisture_time_percentile.csv', )
terra::writeCDF(cpc_data, '/mnt/data1/soil-moisture-models/nc/cpc_soil_moisture_percentile.nc')
