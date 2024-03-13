#this script employs array slicing of netcdf files to extract gridmet data
#at the points of interest and then computes the Standardized Precipitation Index

library(tidyverse)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(foreach)
library(doParallel)
library(magrittr)
library(doSNOW)

#import locations
stations = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station_meta.csv')
station_data = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide.csv') %>%
  pivot_longer(., cols = -c(site_id, date))

#read in soil moisture data to see what the oldest record is
#this tells us how deep of an SPI time series we need to compute
min_date = min(station_data$date)

#compute the unique ids
ids = unique(station_data$site_id)
#compute the dates of interest for each site
start_end_dates = station_data %>%
  group_by(site_id) %>%
  summarise(min_date = min(date),
            max_date = max(date))

#import drought functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

#define an extra ancellary function (this is to do the array slice for the gridmet pixel centroid closest to station)
wherenearest = function(val,matrix){
  dist  = abs(matrix-val)
  index = which.min(dist)
  return( index )
}

#open the data from the nc file
nc_precip  = nc_open('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pr.nc')
#get vectors from nc file for lat and lon
lon.precip = ncvar_get(nc_precip,varid='lon')
lat.precip = ncvar_get(nc_precip,varid='lat')

#get the time ids of the z-dimention
precip_time = read_csv('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pr_time.csv')

#exctract precip vals from disk using array slicing
data_precip = apply(stations,MARGIN=1,FUN=function(x){  
  return(ncvar_get(nc_precip,varid='precipitation_amount',start=c(wherenearest(x['longitude'] %>% as.numeric(),lon.precip),wherenearest(x['latitude'] %>% as.numeric(),lat.precip),1),count=c(1,1,-1))) 
}) %>%
  as_tibble() %>%
  `colnames<-`(c(stations$site_id)) %>%
  mutate(time = precip_time$datetime) %>%
  pivot_longer(cols = -c(time))
#close the nc file
nc_close(nc_precip)

#define function to compute the SPI
compute_spi = function(index, data, timescale){
  data = data %>%
    mutate(day = day(time),
           month = month(time)) %>%
    #slice out data older than index
    slice(1:index)
  
  #compute indexes for time breaks
  first_date_breaks = which(data$day == data$day[index] & data$month == data$month[index])
  second_date_breaks = first_date_breaks-(timescale-1)
  
  #if there are negative indexes remove last year (incomplete data range)
  #change this to remove all indexes from both vectors that are negative
  if(!all(second_date_breaks < 0)){
    pos_index = which(second_date_breaks > 0)
    first_date_breaks = first_date_breaks[c(pos_index)]
    second_date_breaks = second_date_breaks[c(pos_index)]
  }
  
  #create slice vectors and group by vectors which define which data should be used where for the SPI calculation
  for(j in 1:length(first_date_breaks)){
    if(j == 1){
      slice_vec = seq(second_date_breaks[j],first_date_breaks[j], by = 1)
      group_by_vec = rep(j,(first_date_breaks[j] - second_date_breaks[j]+1))
    }
    else{
      slice_vec = append(slice_vec, seq(second_date_breaks[j],first_date_breaks[j], by = 1))
      group_by_vec = append(group_by_vec, rep(j,(first_date_breaks[j] - second_date_breaks[j]+1)))
    }
  }
  
  #pre-process the precipaitaion data before sending it to the SPI function
  processed = data %>%
    slice(slice_vec) %>%
    tibble::add_column(group_by_vec = group_by_vec)%>%
    #group by the group_by_vec
    group_by(group_by_vec)%>%
    #summation of data (and convert to mm)
    dplyr::summarise(sum = sum(value/10, na.rm = T)) %>%
    #add back time
    mutate(time = data$time[first_date_breaks]) %>%
    #select the most recent data
    tail(., 30)

  #take the preprocessed data and compute the SPI  
  final = processed %>%
    #compute spi
    mutate(spi = gamma_fit_spi(processed$sum, export_opts = 'SPI', return_latest = F, climatology_length = 30))
  
  return(final %>% mutate(site_id = data$name[1], timescale = timescale) %>% dplyr::select(site_id,time, timescale, spi))
}

#17 hrs on 30 cores
#start the timer to compute runtime
tictoc::tic()
#rec up the cluster
cl = makeSOCKcluster(30)
registerDoSNOW(cl)
#define the parameters for the progressbar
pb = txtProgressBar(min=1, max=length(ids), style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)
#define the export storage vector
out = list()
#foreach by site
out = foreach(i = 1:length(ids), .packages = c('tidyverse', 'lubridate', 'magrittr', 'Lmoments'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    #compute the end date
    temp_end_date = start_end_dates %>%
      filter(site_id == ids[i]) %$%
      min_date
    
    #select the precip data for the site of interest
    temp_data = data_precip %>%
      filter(name == ids[i])
    
    #compute the indicies of interest
    indicies_of_interest = which(
      temp_data %>%
        mutate(year = lubridate::year(time)) %$%
        year 
      == 2021
    )
    #define the export list
    export = list()
    #define the timescales of interes (10 - 730 by 10)
    timescales_of_interest = c(seq(10,730,10))
    
    #loop through timescales
    for(t in 1:length(timescales_of_interest)){
      #define the name of the timescale for dynamic name definitions
      name = paste0('t_', timescales_of_interest[t])
      #compute the SPI and dynamically name data
      export[[t]] = tibble(apply(indicies_of_interest %>% as.data.frame, MARGIN = 1, FUN = compute_spi, data = temp_data, timescale = timescales_of_interest[t])) %>%
        purrr::map(., bind_rows) %$%
        `apply(...)` %>%
        arrange(time)
    }
    
    #define the final out tibble
    out_final = export %>%
      bind_rows() %>%
      filter(time >= temp_end_date)
    
    out_final
  }, error = function(e){
    out_final = NA
    out_final
  })
  
} 
close(pb)
#stop the cluster and compute runtime
stopCluster(cl)
tictoc::toc()

#save out the data as an RDS (list)
saveRDS(out, '/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spi-data-list-10s.RDS')

#find failed indexes
failed = lapply(out, function(x) length(x)) %>%
  unlist() 

#compute the final tibble binding lists
final = out %>%
  bind_rows() 

#write out the final data
write_csv(final, '/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spi-data-long-10s.csv')