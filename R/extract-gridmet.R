library(tidyverse)
library(ncdf4)

#define ancellary functions
wherenearest = function(val,matrix){
  dist  = abs(matrix-val)
  index = which.min(dist)
  return( index )
}

stations = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-meta.csv')

nc_precip  = nc_open('/home/zhoylman/soil-moisture-validation-data/gridmet/precip_gridmet.nc')
lon.precip = ncvar_get(nc_precip,varid='lon')
lat.precip = ncvar_get(nc_precip,varid='lat')

precip_time = read_csv('/home/zhoylman/soil-moisture-validation-data/gridmet/precip_gridmet_time.csv')

#exctract topo vals from disk
data_precip = apply(stations,MARGIN=1,FUN=function(x){  
  return(ncvar_get(nc_precip,varid='precip',start=c(wherenearest(x['longitude'] %>% as.numeric(),lon.precip),wherenearest(x['latitude'] %>% as.numeric(),lat.precip),1),count=c(1,1,-1))) 
}) %>%
  as_tibble() %>%
  `colnames<-`(c(stations$site_id)) %>%
  mutate(time = precip_time$datetime) %>%
  pivot_longer(cols = -c(time))

nc_close(nc_precip)

