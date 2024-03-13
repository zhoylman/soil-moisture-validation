# download all the uscrn data

library(tidyverse)
library(RCurl)
library(magrittr)
library(sf)
library(data.table)
library(xml2)
library(rvest)

#define special
`%notin%` = Negate(`%in%`)

#define where files are
url = 'https://www.ncei.noaa.gov/pub/data/uscrn/products/soil/soilanom01/'

#compute the file names
files = read_html(url) %>%
  html_nodes(., 'a') %>%
  html_attr(., "href") %>%
  as_tibble() %>%
  filter(str_detect(value, '.csv')) %$%
  value

# download data from uscrn file server 
# only need to do this once...
# #bind for full name
# full_url = paste0(url, files)
# 
# for(i in 1:length(full_url)){
#   download.file(full_url[i], paste0('/home/zhoylman/soil-moisture-validation-data/raw/uscrn/raw/', files[i]))
# }
#download.file('https://www.ncei.noaa.gov/pub/data/uscrn/products/stations.tsv', '/home/zhoylman/soil-moisture-validation-data/raw/uscrn/stations.tsv')

#import raw data (post download)
data = list.files('/home/zhoylman/soil-moisture-validation-data/raw/uscrn/raw', full.names = T) %>%
  as_tibble() %>%
  mutate(short = list.files('/home/zhoylman/soil-moisture-validation-data/raw/uscrn/raw/')) %>%
  filter(short %notin% 'CRNSSM0101-AK_Kenai_29_ENE.csv') %$%
  value %>%
  lapply(., read_csv, show_col_types = FALSE) %>%
  lapply(., function(x){x %>% mutate(WBAN_NO = WBAN_NO %>% as.character())})

#cbind data and select start of day value
full_data = data %>%
  bind_rows() %>%
  mutate(hour = substr(DATE_TIME, 9, 10) %>% as.numeric) %>%
  filter(hour == 0)

#reformat the data
formated = full_data %>%
  mutate(date = substr(DATE_TIME, 1, 8) %>% as.Date(., format = '%Y%m%d')) %>%
  select(WBAN_NO, date, SMVWC_5_CM, SMVWC_10_CM, SMVWC_20_CM, SMVWC_50_CM, SMVWC_100_CM,
         ST_5_CM, ST_10_CM, ST_20_CM, ST_50_CM, ST_100_CM)

#compute name conversion table
name_conversion = tibble(names = c(colnames(formated))) %>%
  mutate(partial_element = ifelse(names %like% 'ST', 'soil_temperature_', 
                                  ifelse(names %like% 'SMVWC', 'soil_moisture_', NA)),
         partial_depth = stringr::str_extract(names, "[[:digit:]]+") %>% as.numeric() / 2.54 %>% round(1),
         new_names = ifelse(names == 'date', 'date',
                            ifelse(names == 'WBAN_NO', 'site_id', paste0(partial_element, partial_depth,'in'))))

#compute file data by renaming using the name conversion table
final = formated %>% 
  rename_at(vars(name_conversion$names), function(x) name_conversion$new_names) %>%
  pivot_longer(cols = -c(site_id, date))

#extract unique sites
sites = unique(final$site_id)

#compute and refomat station metadata
uscrn_meta = read_tsv('/home/zhoylman/soil-moisture-validation-data/raw/uscrn/stations.tsv') %>%
  filter(WBAN %in% sites) %>%
  mutate(site_id = WBAN,
         elevation_ft = ELEVATION %>% as.numeric(),
         latitude = LATITUDE, 
         longitude = LONGITUDE,
         network = 'USCRN') %>%
  select(network, site_id, latitude, longitude, elevation_ft)

#write out!
write_csv(final, '/home/zhoylman/soil-moisture-validation-data/processed/uscrn-soil-moisture/uscrn-soil-moisture-long.csv')
write_csv(uscrn_meta, '/home/zhoylman/soil-moisture-validation-data/processed/uscrn-soil-moisture/uscrn-soil-moisture-meta.csv')