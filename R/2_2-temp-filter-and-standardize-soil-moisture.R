# standardize the soil moisture data using a beta distrobution

library(tidyverse)
library(magrittr)
library(data.table)
library(lubridate)
library(doSNOW)

# source ancillary functions
source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

#define special
`%notin%` = Negate(`%in%`)

#import all data (merged in script 2_1)
vwc_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/soil-moisture-data-wide.csv')
sites = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/station_meta.csv')

#compute the cdf, this can be converted to stand anomoly with a qnorm transform
export_opts_id = 'CDF' # 'SPI'

#define the function to compute the moving window standardization 
#this will be nested within the next funtion
moving_window_standardize = function(x){
  #set min data for drought anom calculations (6 years min)
  min_data_thresh = 31*6
  
  #define indicies of interest
  indicies = 1:length(x$date)
  #define the out vector
  out = vector()
  #iterate through indicies of interest
  for(index in indicies){
    #compute the day of interest for the index
    date_of_interest = x$date[index]
    #compute index specific window (31 day)
    top_window = x$yday[index] + 15
    bottom_window = x$yday[index] - 15
    #correct for yday breaks at 1 and 365 
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
      mutate(drought_anomaly = beta_fit_smi(.[[2]], return_latest = F, climatology_length = Inf, export_opts = export_opts_id))
    
    #find index of centroid date (date of interest)
    #make sure it has enough data as well (6 year minimum)
    if(length(standard$date) >= min_data_thresh){
      out[index] = standard$drought_anomaly[which(standard$date == date_of_interest)]
    } else {
      out[index] = NA
    }
  }
  #memory management
  rm(standard, x); gc(); gc()
  return(out)
}

#function that first filters out frozen soils, then computes the centered moving window 
#standardized values
temp_filter_standardize_vwc = function(site_of_interest, vwc_all){
  #pull in data and select site of interest
  #soil moisture data
  vwc = vwc_all %>%
    filter(site_id == site_of_interest) %>% 
    #remove columns without data
    select_if(~sum(!is.na(.)) > 0)
  
  #evalaute varialbe depths and compute standardized depths
  cols = colnames(vwc) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'date')) %>%
    mutate(depth = parse_number(value),
           generalized_depth = ifelse(depth <= 4, 'Shallow',
                                      ifelse(depth > 4 & depth <= 20, 'Middle',
                                             ifelse(depth > 20, 'Deep', NA))),
           depth = paste0(depth, 'in')) 
  
  
  #if there is data at each depth compute depth averaged soil moisture
  if(length(unique(cols$generalized_depth)) == 3){
    depth_averaged = vwc %>%
      pivot_longer(., cols = -c('site_id', 'date')) %>%
      mutate(depth = parse_number(name),
             generalized_depth = ifelse(depth <= 4, 'Shallow',
                                        ifelse(depth > 4 & depth <= 20, 'Middle',
                                               ifelse(depth > 20, 'Deep', NA))),
             variable = ifelse(str_detect(name, 'soil_moisture'), 'mean_soil_moisture',
                               ifelse(str_detect(name, 'soil_temperature_'), 'mean_soil_temperature',NA))) %>%
      #compute generalized depth average conditions
      group_by(site_id, date, generalized_depth, variable) %>%
      summarise(value = mean(value, na.rm = T)) %>%
      ungroup() %>%
      #now group only by date and variable to compute full depth averaged conditions
      group_by(site_id, date, variable) %>%
      summarise(value = mean(value, na.rm = F)) %>%
      ungroup() %>%
      #add a "depth" id
      drop_na(value) %>%
      pivot_wider(., names_from = 'variable', values_from = 'value')
    
    #join disperate depths to single tibble
    vwc = left_join(vwc, depth_averaged, by = c('site_id', 'date'))
  }

  #if there is mean depth data, add it to the tibble
  if(length(unique(cols$generalized_depth)) == 3){
    cols = cols %>%
      bind_rows(., tibble(value = c('mean_soil_moisture', 'mean_soil_temperature'),
                          depth = c('mean', 'mean'),
                          generalized_depth = c('Mean', 'Mean')))
  }
  
  #compute unique depths
  unique_depths = unique(cols$depth)
  
  #define temp list for storage
  temp = list()
  
  #itterate across depths
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
    
    #filter for temperature (non-frozen) and standardize
    storage_anom = vwc %>%
      dplyr::select(date, temp_moisture_col, temp_temperature_col) %>%
      #filter data for temperature less than 34 F
      filter(!!as.name(temp_temperature_col) >= 34) %>%
      dplyr::select(date, temp_moisture_col) %>%
      # if the data is scalled between 0-100, convert to 0-1
      mutate(!!as.name(temp_moisture_col) := case_when(any(!!as.name(temp_moisture_col) > 1) ~ !!as.name(temp_moisture_col)/100, 
                                                       TRUE ~!!as.name(temp_moisture_col)))%>%
      #drop NAs
      drop_na() %>%
      #compute standardized storage anomoly (no centered window) 
      #this does not account for seasonality
      mutate(!!as.name(paste0('storage_anomaly_',unique_depths[i])) := 
               #gamma_fit_spi(!!as.name(temp_moisture_col), return_latest = F, climatology_length = Inf, export_opts = export_opts_id),
               beta_fit_smi(!!as.name(temp_moisture_col), return_latest = F, climatology_length = Inf, export_opts = export_opts_id),
               yday = yday(date))

    #compute drought anomoly (based on a 31 day centered moving window)
    #this removes the seasonality
    drought_anom = moving_window_standardize(storage_anom)
    
    #store results in the temp list
    temp[[i]] = storage_anom %>%
      mutate(!!as.name(paste0('drought_anomaly_',unique_depths[i])) := 
               drought_anom) %>%
      dplyr::select(-yday)
    print(i)
    #memory management
    rm(drought_anom, storage_anom, temp_temperature_col, 
       temp_moisture_col); gc(); gc()
  }
  #Memory management
  rm(vwc, vwc_all); gc(); gc()
  #compute final export tibble
  export = purrr::reduce(temp, left_join, by = c('date')) %>%
    mutate(site_id = site_of_interest) %>%
    relocate(site_id)
  
  return(export)
}

#~ 24 hours on 20 cores - needs significant RAM (~300Gb)
#there is no l-moment method for beta distrobutions so we need to use
#aximum likelihood which itterates across the error space which increases memory
#usage significantly
#itteration is a bit goofy here, memory drift was a problem here, so we rev up the cluster
#for parallel processing, 20 at a time, then shut down the cluster to clear memory and
#then start up the cluster again. So we itterate over groups of 20, opening and closing the cluster 
#each time. 
full_start = Sys.time()
#define the groups
groups = rep(1:41,each = 20)[1:length(sites$site_id)]
#define the storage "out" list
out_list = list()
#for each group
for(g in unique(groups)){
  tictoc::tic()
  print(paste0(g, ' of ', max(groups)))
  #find the sites assosiated with the batch group
  batch_sites = sites$site_id[which(groups == g)]
  #rev up the cluster ready for 20 sites
  cl = makeSOCKcluster(20)
  registerDoSNOW(cl)
  #foreach loop across the sites in the batch group
  out_list[[g]] = foreach(i = which(groups == g), .packages = c('tidyverse', 'lubridate', 'magrittr', 'data.table')) %dopar% {
    gc()
    #filter for the sites of interest
    vwc_temp = vwc_all %>%
      filter(site_id == sites$site_id[i])
    tryCatch({
      #compute the standardization
      return = temp_filter_standardize_vwc(sites$site_id[i], vwc_temp)
      #memory management!
      rm(vwc_temp); gc(); gc()
      return
    }, error = function(e){
      return = NA
      #memory management!
      rm(vwc_temp); gc(); gc()
      return
    })
  }
  tictoc::toc()
  #shutdown the cluster to open memory
  stopCluster(cl)
}
#compute full run time
full_end = Sys.time()
#print full run time
print(full_end - full_start)

#function to remove lists with only NA
na.omit.list = function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

#compute the final tibble
final = out_list %>%
  unlist(., recursive = FALSE) %>%
  na.omit.list() %>%
  bind_rows()

#filter the final meta for sites that passed all filter checks
sites_final = sites %>%
  filter(site_id %in% unique(final$site_id))

#export!
if(export_opts_id == 'CDF'){
  write_csv(final, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/beta-standardized-soil-moisture-data-wide-6-years-min-CDF-w-mean.csv')
  write_csv(sites_final, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/beta-standardized-station-meta-6-years-min-CDF-w-mean.csv')
}

