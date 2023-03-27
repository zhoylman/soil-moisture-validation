library(tidyverse)
library(sf)

`%notin%` = Negate(`%in%`)

stations_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station_meta.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs('EPSG:4326')) %>%
  st_transform(., st_crs('EPSG:5070')) %>%
  mutate(network = ifelse(network == 'SNTLT', 'SNOTEL', network),
         network = ifelse(network == 'SNTL', 'SNOTEL', network))

states = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') %>%
  st_transform(., st_crs('EPSG:5070')) %>%
  filter(NAME %notin% c('Puerto Rico', 'Alaska', 'Hawaii'))

site_map = ggplot()+
  geom_sf(data = states, fill = 'transparent')+
  geom_sf(data = stations_meta, aes(fill = network), shape = 21)+
  theme_bw(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        plot.subtitle = element_text(hjust = 0.5),
        legend.title=element_blank())+
  scale_fill_manual(values = viridis::turbo(5), name = '')+
  ggtitle(expression(~italic(In)~italic(Situ)~' Soil Moisture Observations'), paste0('n = ', length(stations_meta$network)))+
  guides(fill = guide_legend(override.aes = list(size=5)))

ggsave(site_map, file = '/home/zhoylman/soil-moisture-validation/figs/site_map.png', width = 8, height = 8)
