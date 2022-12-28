library(tidyverse)
library(sf)
library(raster)
library(ncdf4)
library(lubridate)
library(foreach)
library(doParallel)
library(magrittr)

#import locations
stations = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station_meta.csv')
station_data = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide.csv') %>%
  pivot_longer(., cols = -c(site_id, date))

#read in soil moisture data to see what the oldest record is
#this tells us how deep of an SPI time series we need to go
min_date = min(station_data$date)

ids = unique(station_data$site_id)
start_end_dates = station_data %>%
  group_by(site_id) %>%
  summarise(min_date = min(date),
            max_date = max(date))

#define eddi function
eddi_fun = function(x, climatology_length = 30) {
  #define coeffitients
  C0 = 2.515517
  C1 = 0.802853
  C2 = 0.010328
  d1 = 1.432788
  d2 = 0.189269
  d3 = 0.001308
  
  # following Hobbins et al., 2016
  x = as.numeric(x)
  
  #extract the "climatology length from the dataset (assumes x is ordered in time, 1991, 1992, 1993... 2020 etc)
  x = tail(x, climatology_length)
  
  if(all(is.na(x))){
    return(NA)
  } else {
    
    #Rank PET (1 = max)
    rank_1 = rank(-x)
    
    #Calcualte emperical probabilities
    prob = ((rank_1 - 0.33)/(length(rank_1) + 0.33))
    
    #compute W (whaterver that is)
    W = numeric(length(prob))
    for(i in 1: length(prob)){
      if(prob[i] <= 0.5){
        W[i] = sqrt(-2*log(prob[i]))
      } else {
        W[i] = sqrt(-2*log(1 - prob[i]))
      }
    }
    
    #Find indexes which need inverse EDDI sign
    reverse_index = which(prob > 0.5)
    
    #Compute EDDI
    EDDI = W - ((C0 + C1*W + C2*W^2)/(1 + d1*W + d2*W^2 + d3*W^3))
    
    #Reverse sign of EDDI values where prob > 0.5
    EDDI[reverse_index] = -EDDI[reverse_index]
    
    #Return full eddi vector
    return(EDDI)
  }
}


#define ancellary functions
wherenearest = function(val,matrix){
  dist  = abs(matrix-val)
  index = which.min(dist)
  return( index )
}

#extract pet
nc_pet  = nc_open('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pet.nc')
lon.pet = ncvar_get(nc_pet,varid='lon')
lat.pet = ncvar_get(nc_pet,varid='lat')

pet_time = read_csv('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pet_time.csv')

#exctract topo vals from disk
data_pet = apply(stations,MARGIN=1,FUN=function(x){  
  return(ncvar_get(nc_pet,varid='potential_evapotranspiration',start=c(wherenearest(x['longitude'] %>% as.numeric(),lon.pet),wherenearest(x['latitude'] %>% as.numeric(),lat.pet),1),count=c(1,1,-1))) 
}) %>%
  as_tibble() %>%
  `colnames<-`(c(stations$site_id)) %>%
  mutate(time = pet_time$datetime) %>%
  pivot_longer(cols = -c(time))

nc_close(nc_pet)

#testers
i = 1

data = data_pet %>%
  filter(name == ids[i])

#compute index of interest
indicies_of_interest = which(
  data %>%
    mutate(year = lubridate::year(data$time)) %$%
    year 
  == 2021
)

index = indicies_of_interest[32]

timescale = 140

#add in data length conditional (30 years)
compute_eddi = function(index, data, timescale){
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
  
  #create slice vectors and group by vectors
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
  
  final = processed %>%
    #compute eddi
    mutate(eddi = eddi_fun(processed$sum, climatology_length = 30))
  
  return(final %>% mutate(site_id = data$name[1], timescale = timescale) %>% dplyr::select(site_id,time, timescale, eddi))
}

library(doSNOW)
tictoc::tic()
cl = makeSOCKcluster(20)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=length(ids), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
#length(ids)
out = foreach(i = 1:length(ids), .packages = c('tidyverse', 'lubridate', 'magrittr', 'Lmoments', 'lmomco'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    temp_end_date = start_end_dates %>%
      filter(site_id == ids[i]) %$%
      min_date
    
    temp_data = data_pet %>%
      filter(name == ids[i])
    
    indicies_of_interest = which(
      temp_data %>%
        mutate(year = lubridate::year(time)) %$%
        year 
      == 2021
    )
    
    export = list()
    #timescales_of_interest = c(seq(5,25,5), seq(30,80,10), seq(90,360,30), 540, 730)
    timescales_of_interest = c(seq(10,730,10))
    #timescales_of_interest = c(10,20,30)
    
    
    for(t in 1:length(timescales_of_interest)){
      name = paste0('t_', timescales_of_interest[t])
      export[[t]] = tibble(apply(indicies_of_interest %>% as.data.frame, MARGIN = 1, FUN = compute_eddi, data = temp_data, timescale = timescales_of_interest[t])) %>%
        purrr::map(., bind_rows) %$%
        `apply(...)` %>%
        arrange(time)
    }
    
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
stopCluster(cl)
tictoc::toc()

saveRDS(out, '/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/eddi-data-list-10s.RDS')

final = out %>%
  bind_rows() 

write_csv(final, '/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/eddi-data-long-10s.csv')
