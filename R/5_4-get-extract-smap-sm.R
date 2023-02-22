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

#define list to store results
out_years = list()

#loop through years (GEE memory constraint)
for(i in 1:length(years)){
  print(years[i])
  
  #define GEE images
  smap_susm = ee$ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture")$
    select('susm')$
    filter(ee$Filter$calendarRange(years[i], field = "year"))$
    toBands()
  
  smap_ssm = ee$ImageCollection("NASA_USDA/HSL/SMAP10KM_soil_moisture")$
    select('ssm')$
    filter(ee$Filter$calendarRange(years[i], field = "year"))$
    toBands()
  
  #stract at point
  smap_susm_extract = ee_extract(x = smap_susm,
                                 y = stations,
                                 fun = ee$Reducer$mean(),
                                 scale = 10000,
                                 sf = T)
  
  smap_ssm_extract = ee_extract(x = smap_ssm,
                                y = stations,
                                fun = ee$Reducer$mean(),
                                scale = 10000,
                                sf = T)
  
  #clean up results and compute time
  smap_susm_extract_clean = smap_susm_extract %>%
    st_drop_geometry() %>%
    pivot_longer(., cols = -c('network', 'site_id', 'elevation_ft'), values_to = 'susm') %>%
    mutate(time_start = substr(name, 18, 25) %>% as.Date(format = '%Y%m%d'),
           time_end = substr(name, 27, 35) %>% as.Date(format = '%Y%m%d')) %>%
    dplyr::select(-name)
  
  smap_ssm_extract_clean = smap_ssm_extract %>%
    st_drop_geometry() %>%
    pivot_longer(., cols = -c('network', 'site_id', 'elevation_ft'), values_to = 'ssm') %>%
    mutate(time_start = substr(name, 18, 25) %>% as.Date(format = '%Y%m%d'),
           time_end = substr(name, 27, 35) %>% as.Date(format = '%Y%m%d')) %>%
    dplyr::select(-name)
  
  #left join
  final = left_join(smap_ssm_extract_clean, smap_susm_extract_clean, 
                    by = c("network","site_id","elevation_ft","time_start","time_end")) %>%
    dplyr::select(network, site_id, elevation_ft, time_start, time_end, ssm, susm)
  
  #save!
  out_years[[i]] = final
  
}

#write out final csv
write_csv(final_bind, '/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/smap-soil-moisture.csv')
