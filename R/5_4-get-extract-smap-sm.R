#load required libraries
library(reticulate) # allows for python interfacing
library(rgee) # R wrapper for the python GEE library
library(sf) # simple feature library - used for vectors
library(tidyverse) # package for tidy syntax etc
library(geojsonio) # package to send ROI SF objects to GEE

#set up the gee environment
use_condaenv("gee", conda = "auto",required = TRUE)
ee = import("ee")
ee_Initialize(drive = TRUE)

#import stations of interest
stations = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station_meta.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs('EPSG:4326'))

#define years of interest
years = 2015:2022
months = 1:12

index = tibble(date = seq(as.Date('2015-04-01'), as.Date('2022-12-01'), by = 'month')) %>%
  mutate(month = lubridate::month(date),
         year = lubridate::year(date))

#define list to store results
out_data = list()

#loop through months (GEE memory constraint)
for(i in 1:length(index$date)){
  print(index$date[i])
  
  #define GEE images
  smap_rootzone = ee$ImageCollection("NASA/SMAP/SPL4SMGP/007")$
    select('sm_rootzone')$
    filter(ee$Filter$calendarRange(index$year[i], field = "year"))$
    filter(ee$Filter$calendarRange(index$month[i], field = "month"))$
    filter(ee$Filter$calendarRange(1, field = "hour"))$
    toBands()
  
  smap_rootzone_extract = ee_extract(x = smap_rootzone,
                                y = stations,
                                fun = ee$Reducer$mean(),
                                scale = 11000,
                                sf = T)
  
  #clean up results and compute time
  final = smap_rootzone_extract %>%
    st_drop_geometry() %>%
    pivot_longer(., cols = -c('network', 'site_id', 'elevation_ft'), values_to = 'rootzone_sm') %>%
    mutate(date = substr(name, 2, 9) %>% as.Date(format = '%Y%m%d')) %>%
    dplyr::select(-name) %>%
    dplyr::select(network, site_id, elevation_ft, date, rootzone_sm)
  
  #save!
  out_data[[i]] = final
}

final_bind = out_data %>%
  bind_rows()

#write out final csv
write_csv(final_bind, '/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/smap-soil-moisture.csv')
