library(tidyverse)
library(magrittr)
library(data.table)

`%notin%` = Negate(`%in%`)

#download most current data from fcfc-mcoapps.cfc.umt.edu
#system('scp -r zhoylman@fcfc-mcoapps.cfc.umt.edu:/home/zhoylman/mesonet-download-data/* /home/zhoylman/soil-moisture-validation-data/raw/mt-mesonet-flat-files')

#pull in raw mesonet data
base_files = list.files('~/soil-moisture-validation-data/raw/mt-mesonet-flat-files', full.names = T) %>%
  lapply(., read_csv) 

#bind rows to single tibble
base_data = base_files %>%
  bind_rows()

#select vars of interest
select_vars = base_data %>%
  filter(str_detect(name, 'soilt|soilw|precipit|air_temp')) %>%
  mutate(date = datetime %>% as.Date()) 

#compute daily data not in the soil  
daily_data_others = select_vars %>%
  group_by(station_key, date) %>%
  summarise(minimum_air_temperature_degF = min(value[name=="air_temp"]) *  (9/5) + 32,
            max_air_temperature_degF = max(value[name=="air_temp"]) *  (9/5) + 32,
            precipitation_increment_in = sum(value[name=="precipit"])/ 25.4) %>%
  mutate(site_id = station_key) %>%
  ungroup() %>%
  select(-c(station_key)) %>%
  pivot_longer(cols = -c(site_id, date)) %>%
  select(site_id, date, name, value)

#convert names to standardize nomencalature across networks
name_conversion_temp = tibble(names = c('station_key', 'date', unique(select_vars %>%
                                              filter(str_detect(name, 'soilt|soilw')) %$% name))) %>%
  mutate(partial_element = ifelse(names %like% 'soilt', 'soil_temperature_', 
                                  ifelse(names %like% 'soilwc', 'soil_moisture_', NA)),
         partial_depth = stringr::str_extract(names, "[[:digit:]]+"),
         new_names = ifelse(names == 'date', 'date',
                            ifelse(names == 'station_key', 'site_id', paste0(partial_element, partial_depth,'in'))))

#compute daily soil data
daily_data_soil = select_vars %>%
  filter(str_detect(name, 'soilt')) %>%
  group_by(station_key, date, name) %>%
  summarise(value = value[1] *  (9/5) + 32) %>%
  ungroup() %>%
  bind_rows(., select_vars %>%
              filter(str_detect(name, 'soilwc')) %>%
              group_by(station_key, date, name) %>%
              summarise(value = value[1]) %>%
              ungroup()) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename_at(vars(name_conversion_temp$names), function(x) name_conversion_temp$new_names) %>%
  pivot_longer(cols = -c(site_id, date))

#retrieve station data
stations = RCurl::getURL("https://mesonet.climate.umt.edu/api/v2/stations?type=csv&clean=true") %>%
  read_csv() %>%
  filter(station %in% unique(bind_rows(daily_data_others, daily_data_soil)$site_id))

# bind together final data
final_data = bind_rows(daily_data_others, daily_data_soil) %>%
  filter(site_id %in% stations$station) %>%
  arrange(site_id,date,name) %>%
  drop_na()

#write out final daily dataset
write_csv(final_data, '~/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-long-with-precip-and-airtemp.csv')

#write out with precip and temp
final_data_wo_precip_or_T = final_data %>%
  filter(name %notin% c('max_air_temperature_degF', 'minimum_air_temperature_degF', 'precipitation_increment_in'))

#write out without precip and temp
write_csv(final_data_wo_precip_or_T, '/home/zhoylman/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-long.csv')

#clean up station meta
stations_final = stations %>% 
  mutate(network = 'MT Mesonet') %>%
  mutate(elevation_ft = elevation*3.28084) %>%
  select(network, station, latitude, longitude, elevation_ft) %>%
  `colnames<-`(c('network', 'site_id', 'latitude', 'longitude', 'elevation_ft'))

#write out
write_csv(stations_final, '/home/zhoylman/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-meta.csv')
