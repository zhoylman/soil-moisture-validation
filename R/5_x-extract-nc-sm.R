library(tidyverse)
library(ncdf4)
library(terra)
library(sf)

stations = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station-meta-6-years-min-CDF-w-mean.csv')

extract_nc = function(nc_id, nc_path, nc_lon_name, nc_lat_name, nc_var_name, nc_time, 
                      stations, station_lat_name, station_lon_name) {
  
  #define ancellary functions
  wherenearest = function(val,matrix){
    dist  = abs(matrix-val)
    index = which.min(dist)
    return( index )
  }
  
  nc_data = nc_open(nc_path)
  lon_data = ncvar_get(nc_data,varid=nc_lon_name)
  lat_data = ncvar_get(nc_data,varid=nc_lat_name)
  
  data_extract = apply(stations,MARGIN=1,FUN=function(x){  
    return(ncvar_get(nc_data,varid=nc_var_name,start=c(wherenearest(x[station_lon_name] %>% as.numeric(),lon_data),wherenearest(x[station_lat_name] %>% as.numeric(),lat_data),1),count=c(1,1,-1))) 
  }) %>%
    as_tibble() %>%
    `colnames<-`(c(stations$site_id)) %>%
    mutate(time = nc_time,
           nc_id = nc_id) %>%
    dplyr::select(nc_id, time, everything())
  
  nc_close(nc_data)
  return(data_extract)
}

#extract SPoRT
SPoRT_extract = extract_nc(nc_id = 'SPoRT_mean_soil_moisture_0-100cm',
                           nc_path = '/mnt/data1/soil-moisture-models/nc/SPoRT__mean_soil_moisture_0-100cm.nc', 
                           nc_lon_name = 'longitude', 
                           nc_lat_name = 'latitude', 
                           nc_var_name = 'SPoRT_mean_soil_moisture_0-100cm', 
                           nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/SPoRT__mean_soil_moisture_0-100cm.csv')$time, 
                           stations = stations, 
                           station_lat_name = 'latitude', 
                           station_lon_name = 'longitude')

write_csv(SPoRT_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/SPoRT-soil-moisture.csv')

#extract cpc
cpc_extract = extract_nc(nc_id = 'cpc_soil_moisture_percentile',
                           nc_path = '/mnt/data1/soil-moisture-models/nc/cpc_soil_moisture_percentile.nc', 
                           nc_lon_name = 'longitude', 
                           nc_lat_name = 'latitude', 
                           nc_var_name = 'cpc_soil_moisture_percentile', 
                           nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/cpc_soil_moisture_time_percentile.csv')$date, 
                           stations = stations, 
                           station_lat_name = 'latitude', 
                           station_lon_name = 'longitude')

#cpc has percentiles greater than 100 and less than 0!!! Not good!
cpc_extract = cpc_extract %>% mutate_at(vars(-c(nc_id, time)), ~replace(., . > 100, 100))
cpc_extract = cpc_extract %>% mutate_at(vars(-c(nc_id, time)), ~replace(., . < 0, 0))

write_csv(cpc_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/cpc-soil-moisture-percentile.csv')

#extract grace
grace_extract = extract_nc(nc_id = 'grace_rtzn_soil_moisture',
                         nc_path = '/mnt/data1/soil-moisture-models/nc/grace_rtzn_soil_moisture.nc', 
                         nc_lon_name = 'longitude', 
                         nc_lat_name = 'latitude', 
                         nc_var_name = 'grace_rtzn_soil_moisture', 
                         nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/grace_rtzn_soil_moisture_time.csv')$date, 
                         stations = stations, 
                         station_lat_name = 'latitude', 
                         station_lon_name = 'longitude')

write_csv(grace_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/grace-soil-moisture-percentile.csv')

# extract the NLDAS2 models
#VIC
nldas2_vic_extract = extract_nc(nc_id = 'NLDAS2_VIC_soil_moisture_0-100cm',
                           nc_path = '/mnt/data1/soil-moisture-models/nc/NLDAS2_VIC_soil_moisture_0-100cm.nc', 
                           nc_lon_name = 'longitude', 
                           nc_lat_name = 'latitude', 
                           nc_var_name = 'nldas2_vic_soil_moisture_0-100cm_group_1', 
                           nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/NLDAS2_VIC_soil_moisture_0-100cm_time.csv')$date, 
                           stations = stations, 
                           station_lat_name = 'latitude', 
                           station_lon_name = 'longitude')

write_csv(nldas2_vic_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_VIC_soil_moisture_0-100cm.csv')

#NOAH
nldas2_noah_extract = extract_nc(nc_id = 'NLDAS2_NOAH_soil_moisture_0-100cm',
                                nc_path = '/mnt/data1/soil-moisture-models/nc/NLDAS2_NOAH_soil_moisture_0-100cm.nc', 
                                nc_lon_name = 'longitude', 
                                nc_lat_name = 'latitude', 
                                nc_var_name = 'nldas2_noah_soil_moisture_0-100cm_group_1', 
                                nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/NLDAS2_NOAH_soil_moisture_0-100cm_time.csv')$date, 
                                stations = stations, 
                                station_lat_name = 'latitude', 
                                station_lon_name = 'longitude')

write_csv(nldas2_noah_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_NOAH_soil_moisture_0-100cm.csv')

#MOSAIC
nldas2_mosaic_extract = extract_nc(nc_id = 'NLDAS2_MOSAIC_soil_moisture_0-100cm',
                                nc_path = '/mnt/data1/soil-moisture-models/nc/NLDAS2_MOSAIC_soil_moisture_0-100cm.nc', 
                                nc_lon_name = 'longitude', 
                                nc_lat_name = 'latitude', 
                                nc_var_name = 'nldas2_mosaic_soil_moisture_0-100cm_group_1', 
                                nc_time = read_csv('/mnt/data1/soil-moisture-models/nc/NLDAS2_MOSAIC_soil_moisture_0-100cm_time.csv')$date, 
                                stations = stations, 
                                station_lat_name = 'latitude', 
                                station_lon_name = 'longitude')

write_csv(nldas2_mosaic_extract, file = '~/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_MOSAIC_soil_moisture_0-100cm.csv')
