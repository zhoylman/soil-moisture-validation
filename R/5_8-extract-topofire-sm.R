library(tidyverse)
library(ncdf4)
library(terra)
library(magrittr)
library(doParallel)
library(foreach)
library(parallel)
library(doSNOW)

files = list.files('/mnt/beegfs/data/public_data/topofire_weather/daily_grids', full.names = T, pattern = 'soil') %>%
	as_tibble() %>%
	mutate(year = parse_number(value) %>% substr(., 1,4) %>% as.numeric()) %>%
	filter(year > 1990) 

files = files$value

stations = read_csv('/home/zach.hoylman/soil-moisture-validation/standardized-station_meta.csv')

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


cl = makeSOCKcluster(20)
registerDoSNOW(cl)
pb = txtProgressBar(min=1, max=length(files), style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
out = list()

topofire_extract = list()

topofire_extract = foreach(i = 1:length(files), .packages = c('tidyverse', 'ncdf4', 'terra'), .options.snow=opts) %dopar% {
	tryCatch({
	temp_time = rast(files[i]) %>% 
			names() %>%
			parse_number() %>%
			as.Date(., origin="1979-01-01")

	#extract topofire
	topofire_extract_temp = extract_nc(nc_id = 'topofire_soil_moisture',
		                   nc_path = files[i], 
		                   nc_lon_name = 'longitude', 
		                   nc_lat_name = 'latitude', 
		                   nc_var_name = 'soil', 
		                   nc_time = temp_time, 
		                   stations = stations, 
		                   station_lat_name = 'latitude', 
		                   station_lon_name = 'longitude')

	write_csv(topofire_extract_temp, paste0('/home/zach.hoylman/soil-moisture-validation/data/topofire_',i,'.csv'))

	topofire_extract_temp
	},
        error=function(cond) {

	topofire_extract_temp = NA
	topofire_extract_temp

	})
}


stopCluster(cl)


binded = topofire_extract %>% bind_rows()
write_csv(binded, file = '/home/zach.hoylman/soil-moisture-validation/topofire-soil-moisture.csv')

