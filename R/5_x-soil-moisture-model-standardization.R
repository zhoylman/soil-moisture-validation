library(tidyverse)
library(magrittr)
library(doParallel)
library(foreach)

source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

#define function to do the moving window drought anom
moving_window_standardize = function(x, data_gap = 1){
  #data_gap represents the gap between obs, here daily = 1, weekly = 7 and so on
  #set min data for drought anom calculations (3 years min)
  min_data_thresh = 93/data_gap

  #compute some date info
  x = x %>%
    dplyr::mutate(date = time,
                  yday = lubridate::yday(date))
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
      mutate(drought_anomaly = gamma_fit_spi(value, return_latest = F, climatology_length = Inf, export_opts = 'SPI'))

    #find index of centroid date (date of interest)
    if(length(standard$date) >= min_data_thresh){
      out[index] = standard$drought_anomaly[which(standard$date == date_of_interest)]
    } else {
      out[index] = NA
    }
  }
  return(out)
}

parallel_standardize = function(data, unique_sites, data_gap){
  out_list = foreach(sites = unique_sites, 
                     .packages = c('tidyverse', 'magrittr')) %dopar% {
                       temp_x = data %>%
                         filter(name == sites)
                       
                       out = temp_x %>%
                         mutate(drought_anomaly = moving_window_standardize(temp_x, data_gap),
                                storage_anomoly = gamma_fit_spi(value, climatology_length = Inf, export_opts = 'SPI', return_latest = F))
                       
                       out
                     }
  return(out_list)
}

## Initialize cluster
cl = makeCluster(10)
registerDoParallel(cl)
clusterExport(cl, c("moving_window_standardize", "gamma_fit_spi"))

## compute grace standardized
grace = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/grace-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

grace_sites = unique(grace$name)

standardized_grace = parallel_standardize(data = grace, unique_sites = grace_sites, data_gap = 7) %>%
  bind_rows()

write_csv(standardized_grace, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/grace-soil-moisture-standardized-w-storage.csv')

## compute cpc standardized
cpc = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/cpc-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

cpc_sites = unique(grace$name)

standardized_cpc = parallel_standardize(data = cpc, unique_sites = cpc_sites, data_gap = 1) %>%
  bind_rows()

write_csv(standardized_cpc, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/cpc-soil-moisture-standardized-w-storage.csv')

## compute SPoRT standardized
SPoRT = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/SPoRT-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

SPoRT_sites = unique(SPoRT$name)

standardized_SPoRT = parallel_standardize(data = SPoRT, unique_sites = SPoRT_sites, data_gap = 1) %>%
  bind_rows()

write_csv(standardized_SPoRT, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/SPoRT-soil-moisture-standardized-w-storage.csv')

## SMAP
smap_raw = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/smap-soil-moisture.csv') %>%
  mutate(nc_id = 'SMAP_rootzone_soil_moisture',
         time = date,
         value = rootzone_sm,
         name = site_id) %>%
  select(nc_id, time, name, value)

smap_sites = unique(smap_raw$name)

standardized_smap = parallel_standardize(data = smap_raw, unique_sites = smap_sites, data_gap = 1) %>%
  bind_rows()

write_csv(standardized_smap, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/smap-soil-moisture-standardized-w-storage.csv')

## compute TOPOFIRE standardized
topofire = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/topofire-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

topofire_sites = unique(topofire$name)

standardized_topofire= parallel_standardize(data = topofire, unique_sites = topofire_sites, data_gap = 1) %>%
  bind_rows()

write_csv(standardized_topofire, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/topofire-soil-moisture-standardized-w-storage.csv')


#stop cluster
stopCluster(cl)
