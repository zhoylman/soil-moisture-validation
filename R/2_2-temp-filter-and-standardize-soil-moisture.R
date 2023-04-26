library(tidyverse)
library(magrittr)
library(data.table)
library(lubridate)
library(doSNOW)

# source ancillary functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

`%notin%` = Negate(`%in%`)

vwc_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/soil-moisture-data-wide.csv')
sites = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/station_meta.csv')

#compute the cdf, this can be converted to stand anomoly with a qnorm transform
export_opts_id = 'CDF' # 'SPI'

moving_window_standardize = function(x){
  #set min data for drought anom calculations (6 years min)
  min_data_thresh = 31*6
  
  indicies = 1:length(x$date)
  out = vector()
  for(index in indicies){
    date_of_interest = x$date[index]
    #compute index specific window
    top_window = x$yday[index] + 15
    bottom_window = x$yday[index] - 15
    #correct for yday breaks ar 1 and 365
    if(top_window > 365){
      top_window = top_window - 365
      top_range = c(x$yday[index]:365,1:top_window)
    } else {
      top_range = x$yday[index]:top_window
    }
    if(bottom_window < 1){
      bottom_window = bottom_window + 365
      bottom_range = c(bottom_window:365, 1:x$yday[index])
    } else {
      bottom_range = bottom_window:x$yday[index]
    }
    #compute range of y days
    range = c(bottom_range, top_range) %>% unique()
    
    #filter data for index
    standard = x %>%
      filter(yday %in% range) %>%
      #return all data (static reference frame) and long cliamtology length
      mutate(drought_anomaly = gamma_fit_spi(.[[2]], return_latest = F, climatology_length = Inf, export_opts = export_opts_id))
    
    #find index of centroid date (date of interest)
    if(length(standard$date) >= min_data_thresh){
      out[index] = standard$drought_anomaly[which(standard$date == date_of_interest)]
    } else {
      out[index] = NA
    }
  }
  return(out)
}
temp_filter_standardize_vwc = function(site_of_interest, vwc_all){
  
  #pull in data and select site of interest
  #soil moisture data
  vwc = vwc_all %>%
    filter(site_id == site_of_interest) %>% 
    #remove columns without data
    select_if(~sum(!is.na(.)) > 0)
  
  #evalaute varialbe depths
  cols = colnames(vwc) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'date')) %>%
    mutate(depth = parse_number(value)) 
  
  #compute unique depths
  unique_depths = unique(cols$depth)
  
  temp = list()
  
  for(i in 1:length(unique_depths)){
    #extract moisture column name for specific depth
    temp_moisture_col = cols %>% 
      filter(depth == unique_depths[i],
             value %like% 'moisture') %$%
      value
    
    #extract temperature column name for specific depth
    temp_temperature_col = cols %>% 
      filter(depth == unique_depths[i],
             value %like% 'temperature') %$%
      value
    
    #filter for temperature and standardize
    storage_anom = vwc %>%
      select(date, temp_moisture_col, temp_temperature_col) %>%
      #filter data for temperature less than 34 F
      filter(!!as.name(temp_temperature_col) >= 34) %>%
      select(date, temp_moisture_col) %>%
      #drop NAs
      drop_na() %>%
      #compute standardized storage anomoly
      mutate(!!as.name(paste0('storage_anomaly_',unique_depths[i], "in")) := 
               gamma_fit_spi(!!as.name(temp_moisture_col), return_latest = F, climatology_length = Inf, export_opts = export_opts_id),
             yday = yday(date))

    #compute drought anomoly (based on a 31 day centered moving window)
    drought_anom = moving_window_standardize(storage_anom)
    
    temp[[i]] = storage_anom %>%
      mutate(!!as.name(paste0('drought_anomaly_',unique_depths[i], "in")) := 
               drought_anom) %>%
      select(-yday)
    print(i)
  }
  export = purrr::reduce(temp, left_join, by = c('date')) %>%
    mutate(site_id = site_of_interest) %>%
    relocate(site_id)
  return(export)
}

#3.3 hr on 20 cores
tictoc::tic()
cl = makeSOCKcluster(30)
registerDoSNOW(cl)
pb = txtProgressBar(min=1, max=length(sites$site_id), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
#length(ids)
out = foreach(i = 1:length(sites$site_id), .packages = c('tidyverse', 'lubridate', 'magrittr', 'data.table'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    return = temp_filter_standardize_vwc(sites$site_id[i], vwc_all)
    return
  }, error = function(e){
    return = NA
    return
  })
} 
close(pb)
stopCluster(cl)
tictoc::toc()

#function to remove lists with only NA
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

final = out %>%
  na.omit.list() %>%
  bind_rows()

sites_final = sites %>%
  filter(site_id %in% unique(final$site_id))

if(export_opts_id == 'CDF'){
  write_csv(final, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide-6-years-min-CDF.csv')
  write_csv(sites_final, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station-meta-6-years-min-CDF.csv')
}
