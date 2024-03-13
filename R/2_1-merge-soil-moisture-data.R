# script to merge datasets and standardize units

# load libs
library(tidyverse)
library(spData)
library(sf)
library(magrittr)
library(data.table)

# define special
`%notin%` = Negate(`%in%`)

# there are some bad sites with soil moisture readings related to saturated conditions and wetlands
# we filter these here
sites_to_remove = c('whitshaw', 'mdaglasw') # wetlands = bad soil moisture data 
# filter out blm1arge late data for 20in

# import mesonet data
# temperature is in F
mt_mesonet = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-long.csv')%>%
  filter(site_id %notin% sites_to_remove,
         !row_number() %in% which(site_id == 'blm1arge' & date>as.Date('2021-01-01') & name == 'soil_moisture_20in'),
         name %notin% c('soil_moisture_00in', 'soil_temperature_00in')) %>%
  arrange(site_id, date, name) %>%
  mutate(name = ifelse(name == "soil_moisture_04in", "soil_moisture_4in", 
                       ifelse(name == 'soil_moisture_08in', 'soil_moisture_8in', 
                              ifelse(name == 'soil_temperature_04in', 'soil_temperature_4in', 
                                     ifelse(name == 'soil_temperature_08in', 'soil_temperature_8in', name)))))

# import NRCS data
# temperature is in F
nrcs = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-long-soil-data-only.csv') %>%
  arrange(site_id, date, name)

# import USCRN data
# temperature is in C - convert to F
uscrn = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/uscrn-soil-moisture/uscrn-soil-moisture-long.csv')%>%
  arrange(site_id, date, name) %>% 
  #convert to F
  mutate(value = ifelse(name %like% 'soil_temperature', ((value * 9/5) + 32), value))

#bind them
almost_final_data = bind_rows(mt_mesonet, nrcs, uscrn)

#compute which ones have soil moisture data
with_data = almost_final_data  %>% 
  filter(str_detect(name, 'soil_moisture')) %>%
  group_by(site_id) %>%
  summarise(n = length(value))

#define final data tibble
final_data = almost_final_data %>%
  filter(site_id %in% with_data$site_id)

#import meta for each network
mt_mesonet_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/mt-mesonet-soil-moisture/mt-mesonet-soil-moisture-meta.csv')%>%
  filter(site_id %notin% sites_to_remove) %>%
  arrange(site_id)

nrcs_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-meta.csv') %>%
  arrange(site_id)

uscrn_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/uscrn-soil-moisture/uscrn-soil-moisture-meta.csv')%>%
  arrange(site_id)

# bind the meta
final_meta = bind_rows(mt_mesonet_meta, nrcs_meta, uscrn_meta) %>%
  filter(site_id %in% unique(final_data$site_id))

#write final data
write_csv(final_meta, '/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/station_meta.csv')
write_csv(final_data, '/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/soil-moisture-data-long.csv')

final_data_wide = final_data%>%
  pivot_wider(names_from = name, values_from = value)

write_csv(final_data_wide, '/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/soil-moisture-data-wide.csv')