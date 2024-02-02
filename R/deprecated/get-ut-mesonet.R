library(tidyverse)
library(RCurl)
library(rjson)
library(magrittr)
library(data.table)

stations = getURL('https://climate.usu.edu/API/api.php/v3/key=XyJM0NexLqHUJjfJSeXkb7SaT8speO/station_search/source=UCC/state=UT') %>%
  fromJSON %$%
  payload %>%
  bind_rows() 

out = list()

for(i in 1:length(stations$station_id)){
  tryCatch({
    out[[i]] = 
      paste0('https://climate.usu.edu/API/api.php/v3/key=XyJM0NexLqHUJjfJSeXkb7SaT8speO/station_search', 
             '/network=',stations$network[i],
             '/station_id=',stations$station_id[i],
             '/get_daily/start_date=2020-9-9',
             '/end_date=2020-10-12',
             '/units=m') %>% 
      getURL() %>%
      fromJSON %$%
      payload %>%
      bind_rows() 
  }, error = function(e){
    out[[i]] = NA
  })
  print(i)
}

soil_id = function(x){
  return(c(x %like% 'soil') %>%
           any(. == TRUE))
}

cleaner = function(x){
  if(length(x) == 0){
    return(NA)
  } else {
    return(x)
  }
}

sm_stations = lapply(out, colnames) %>%
  lapply(., soil_id) 

index = which(sm_stations %>% unlist == TRUE)

test = stations[index,] %>%
  filter(network %in% c('AGWX', 'UCRN', 'UAGRIMET'))
