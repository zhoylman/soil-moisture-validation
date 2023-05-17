library(tidyverse)
library(magrittr)
library(doParallel)
library(foreach)

is_na <- function(x) {
  anyNA(x)
}

source('https://raw.githubusercontent.com/mt-climate-office/mco-drought-indicators/master/processing/ancillary-functions/R/drought-functions.R')

#define function to do the moving window drought anom
moving_window_standardize = function(x, data_gap = 1, export_opts = 'CDF'){
  #data_gap represents the gap between obs, here daily = 1, weekly = 7 and so on
  #set min data for drought anom calculations (6 years min)
  min_data_thresh = (31*6)/data_gap

  #compute some date info
  x = x %>%
    dplyr::mutate(date = time,
                  yday = lubridate::yday(date),
                  year = lubridate::year(date)) %>%
    arrange(date)
  
  
  #compute which years have a full dataset
  full_years = x %>%
    group_by(year) %>%
    summarise(n = length(value))
    
  first_full_year = full_years$year[which(full_years$n == 365)[1]]
  
  indicies = which(x$year == first_full_year)
  
  out = list()
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
      mutate(drought_anomaly = gamma_fit_spi(value, return_latest = F, climatology_length = Inf, export_opts = 'CDF'))

    #find index of centroid date (date of interest)
    if(length(standard$date) >= min_data_thresh){
      out[[index]] = standard[which(standard$yday == lubridate::yday(date_of_interest)),]
    } else {
      out[[index]] = NA
    }
  }
  #if its the last index
  if(index == indicies[length(indicies)]){
    #function to remove empty elements from list
    is_na <- function(x) {
      anyNA(x)
    }
    #bind the final data frame together
    out_final = out %>% 
      purrr::discard(., is_na) %>%
      bind_rows() %>%
      dplyr::arrange(time)
  }
  return(out_final)
}

parallel_standardize = function(data, unique_sites, data_gap, export_opts = 'CDF'){
  out_list = foreach(sites = unique_sites, 
                     .packages = c('tidyverse', 'magrittr', 'ggplot2')) %dopar% {
                       tryCatch({
                         export_opts_id = 'CDF'
                         
                         temp_x = data %>%
                           filter(name == sites)
                         
                         drought_anomaly = moving_window_standardize(temp_x, data_gap) %>%
                           dplyr::select(time, drought_anomaly)
                         
                         out = temp_x %>%
                           left_join(., drought_anomaly, by = 'time') %>%
                           mutate(storage_anomoly = gamma_fit_spi(value, climatology_length = Inf, export_opts = 'CDF', return_latest = F))
                         
                         # plot = ggplot(out, aes(x = storage_anomoly, y = drought_anomaly, color = lubridate::yday(time)))+
                         #   geom_point()+
                         #   ggtitle(out$name[1], out$nc_id[1])+
                         #   scale_color_gradientn(colors = rainbow(365))
                         # 
                         # ggsave(plot, file = paste0('/home/zhoylman/temp/dist_plots/', out$name[1],'_',out$nc_id[1], '.png'))
                         # 
                         out
                       }, error = function(e){
                         out = NA
                         out
                       })
                     }
  return(out_list)
}

####################### Already in Percentiles ##############################

## reorganize grace 
grace = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/grace-soil-moisture-percentile.csv') %>%
  pivot_longer(cols = -c(nc_id, time)) %>%
  mutate(drought_anomaly = value,
         storage_anomoly = NA)

write_csv(grace, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/grace-soil-moisture-standardized-percentile.csv')

## reorganize cpc 
cpc = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/cpc-soil-moisture-percentile.csv') %>%
  pivot_longer(cols = -c(nc_id, time)) %>%
  mutate(drought_anomaly = value,
         storage_anomoly = NA)

write_csv(cpc, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/cpc-soil-moisture-standardized-percentile.csv')

## reorganize SMAP
smap = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/smap-soil-moisture-percentile.csv') %>%
  mutate(nc_id = 'SMAP_rootzone_soil_moisture',
         time = date,
         value = sm_rootzone_pctl,
         name = site_id) %>%
  select(nc_id, time, name, value) %>%
  mutate(drought_anomaly = value,
         storage_anomoly = NA)

write_csv(smap, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/smap-soil-moisture-standardized-percentile.csv')

############################################################################

####################### Percentiles Unavailable ############################
#here we will compute our own according to procedure in script 2_2 
#standardization method for observed data. 

## Initialize cluster
cl = makeCluster(30)
registerDoParallel(cl)
clusterExport(cl, c("moving_window_standardize", "gamma_fit_spi"))

## compute SPoRT standardized
SPoRT = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/SPoRT-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

SPoRT_sites = unique(SPoRT$name)

standardized_SPoRT = parallel_standardize(data = SPoRT, unique_sites = SPoRT_sites, data_gap = 1) 

standardized_SPoRT_bind = Filter(function(a) any(!is.na(a)), standardized_SPoRT) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_SPoRT_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/SPoRT-soil-moisture-standardized-percentile.csv')

## compute TOPOFIRE standardized
topofire_raw = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/topofire-soil-moisture.csv') %>%
  pivot_longer(cols = -c(nc_id, time)) %>%
  arrange(time)

#filter for the most recent 30 years
topofire_years = topofire_raw %>%
  summarize(years = unique(lubridate::year(time)))

#extract the last 30 years of data to be consistant with hoylman et al., 2022
topofire = topofire_raw %>%
  filter(lubridate::year(time) %in% tail(topofire_years$years, 30))

topofire_sites = unique(topofire$name)

standardized_topofire = parallel_standardize(data = topofire, unique_sites = topofire_sites, data_gap = 1) 

standardized_topofire_bind = Filter(function(a) any(!is.na(a)), standardized_topofire) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_topofire_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/topofire-soil-moisture-standardized-percentile.csv')

######################## NLDAS2 ##################################
# starting with VIC
## compute NLDAS2 standardized
nldas2_vic = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_VIC_soil_moisture_0-100cm.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

nldas2_vic_sites = unique(nldas2_vic$name)

tictoc::tic()
standardized_nldas2_vic = parallel_standardize(data = nldas2_vic, unique_sites = nldas2_vic_sites, data_gap = 1) 
tictoc::toc()

standardized_nldas2_vic_bind = Filter(function(a) any(!is.na(a)), standardized_nldas2_vic) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_nldas2_vic_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/NLDAS2-VIC-soil-moisture-standardized-percentile.csv')

# next with NOAH
## compute NLDAS2 standardized
nldas2_noah = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_NOAH_soil_moisture_0-100cm.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

nldas2_noah_sites = unique(nldas2_noah$name)

tictoc::tic()
standardized_nldas2_noah = parallel_standardize(data = nldas2_noah, unique_sites = nldas2_noah_sites, data_gap = 1) 
tictoc::toc()

standardized_nldas2_noah_bind = Filter(function(a) any(!is.na(a)), standardized_nldas2_noah) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_nldas2_noah_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/NLDAS2-NOAH-soil-moisture-standardized-percentile.csv')

# then with MOSAIC
## compute NLDAS2 standardized
nldas2_mosaic = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/soil-moisture-model-extractions/NLDAS2_MOSAIC_soil_moisture_0-100cm.csv') %>%
  pivot_longer(cols = -c(nc_id, time))

nldas2_mosaic_sites = unique(nldas2_mosaic$name)

tictoc::tic()
standardized_nldas2_mosaic = parallel_standardize(data = nldas2_mosaic, unique_sites = nldas2_mosaic_sites, data_gap = 1) 
tictoc::toc()

standardized_nldas2_mosaic_bind = Filter(function(a) any(!is.na(a)), standardized_nldas2_mosaic) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_nldas2_mosaic_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/NLDAS2-MOSAIC-soil-moisture-standardized-percentile.csv')

#finally, based on the ensamble mean
nldas2_ensamble = nldas2_vic %>%
  select(time,name,value) %>%
  rename(vic = value) %>%
  left_join(., nldas2_noah %>%
              select(time,name,value) %>%
              rename(noah = value)) %>%
  left_join(., nldas2_mosaic %>%
              select(time,name,value) %>%
              rename(mosaic = value)) %>%
  rowwise %>% 
  mutate(value = mean(c(vic, noah, mosaic))) %>%
  select(time, name, value) %>%
  mutate(nc_id = 'NLDAS2_ensamble_soil_moisture_0-100cm') %>%
  relocate(nc_id) %>%
  drop_na(value)

nldas2_ensamble_sites = unique(nldas2_ensamble$name)

tictoc::tic()
standardized_nldas2_ensamble = parallel_standardize(data = nldas2_ensamble, unique_sites = nldas2_ensamble_sites, data_gap = 1) 
tictoc::toc()

standardized_nldas2_ensamble_bind = Filter(function(a) any(!is.na(a)), standardized_nldas2_ensamble) %>%
  bind_rows() %>%
  mutate(drought_anomaly = drought_anomaly*100,
         storage_anomoly = storage_anomoly*100)

write_csv(standardized_nldas2_ensamble_bind, '/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models/NLDAS2-ENSAMBLE-soil-moisture-standardized-percentile.csv')

#stop cluster
stopCluster(cl)
