library(tidyverse)
library(terra)
library(httr)
library(xml2)
library(rvest)

files = read_html('https://nasagrace.unl.edu/data/') 

dates = html_attr(html_nodes(files, "a"), "href") %>%
  gsub("\\D", "", .) %>%
  as.numeric() %>%
  as_tibble() %>%
  drop_na() %>%
  mutate(url = paste0('https://nasagrace.unl.edu/data/', value, 
                      '/rtzsm_perc_0125deg_US_', value, '.tif'))

for(i in 1:length(dates$url)){
  GET(dates$url[i], write_disk(paste0('/mnt/data1/soil-moisture-models/grace/rtzn/grace_rtzn_soil_moisture_', dates$value[i], 
                                      '.tif'),
                               overwrite = T))
}

data = list.files('/mnt/data1/soil-moisture-models/grace/rtzn/', full.names = T) %>%
  purrr::map(., rast) %>%
  rast

grace_time = list.files('/mnt/data1/soil-moisture-models/grace/rtzn/', full.names = F) %>%
  gsub("\\D", "", .) %>%
  as_tibble() %>%
  mutate(date = as.Date(value, format = '%Y%m%d'))

write_csv(grace_time, file = '/mnt/data1/soil-moisture-models/nc/grace_rtzn_soil_moisture_time.csv', )
terra::writeCDF(data, '/mnt/data1/soil-moisture-models/nc/grace_rtzn_soil_moisture.nc')
