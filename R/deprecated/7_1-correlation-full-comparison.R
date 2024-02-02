library(tidyverse)

spi = read_csv('~/soil-moisture-validation-data/processed/rmse-comparison/spi-6-years-min-season-rmse.csv') %>%
  mutate(drought_metric = 'SPI')
spei = read_csv('~/soil-moisture-validation-data/processed/rmse-comparison/spei-6-years-min-season-rmse.csv')%>%
  mutate(drought_metric = 'SPEI')
eddi = read_csv('~/soil-moisture-validation-data/processed/rmse-comparison/eddi-6-years-min-season-rmse.csv')%>%
  mutate(drought_metric = 'EDDI',
         median_seasonal_rigid_r = median_seasonal_rigid_r*-1,
         median_seasonal_flexible_r = median_seasonal_flexible_r * -1)

drought_metrics = bind_rows(list(spi, spei, eddi)) %>%
  #select(site_id, drought_metric, generalized_depth, median_seasonal_rigid_rmse, median_seasonal_flexible_rmse) %>%
  mutate(rmse = median_seasonal_rigid_r)%>%
  select(site_id, drought_metric, generalized_depth, rmse) 

soil_moisture_model = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/rmse-comparison/soil-moisture-model-r.csv') %>%
  rename(drought_metric = 'nc_id',
         generalized_depth = 'depth') %>%
  mutate(rmse = pearson_r) %>%
  drop_na(rmse) %>%
  mutate(generalized_depth = ifelse(generalized_depth == "Middle (8-20in)", "Middle (8in - 20in)", generalized_depth)) %>%
  select(site_id, drought_metric, generalized_depth, rmse) 
  
merged = bind_rows(drought_metrics, soil_moisture_model)

#compute density data
density_data = merged %>%
  group_by(generalized_depth, drought_metric) %>%
  do(density_x = density(.$rmse, bw = 'nrd')$x,
     density_y = density(.$rmse, bw = 'nrd')$y,
     density_bw = density(.$rmse, bw = 'nrd')$bw,
     n = length(.$rmse)) %>%
  unnest(cols = c(density_x, density_y, n, density_bw)) %>%
  mutate(generalized_depth = factor(generalized_depth, levels = c('Shallow (0-4in)', 'Middle (8in - 20in)', 'Deep (>20in)')))  %>%
  mutate(drought_metric = factor(drought_metric, 
                                 levels = c("SPI", "SPEI", "EDDI",
                                            "CPC Soil Moisture", "GRACE Rootzone Soil Moisture",
                                            "SMAP (L4) Rootzone Soil Moisture",
                                            "SPoRT 0-100cm Soil Moisture",
                                            "Topofire Soil Moisture")))
  
mode_data = density_data %>%
  group_by(generalized_depth, drought_metric) %>%
  filter(density_y == max(density_y)) %>%
  mutate(modal_r = plyr::round_any(density_x, 0.01))

density_plot = ggplot(data = density_data)+
  geom_segment(data = mode_data, aes(x = density_x, xend = density_x, y = density_y, yend = 0, color = drought_metric), alpha = 0.6)+
  geom_line(aes(x = density_x, y = density_y, color = drought_metric))+
  ggtitle(paste0('Distrobution of Station Specific Correlation (r)'), 'All Sites')+
  theme_bw(base_size = 14)+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  geom_point(data = mode_data, aes(x = density_x, y = density_y, color = drought_metric))+
  facet_grid(generalized_depth~.)+
  labs(x = 'Correlation (r)', y = 'Density')+
  scale_colour_manual(values = viridis::turbo(8), name = '')+
  xlim(-0.5,1)+
  theme(legend.position="bottom", legend.box = "horizontal",
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.background = element_rect(size = 0.2))+
  guides(color=guide_legend(nrow=4,byrow=TRUE))

density_plot

ggsave(density_plot, file = '/home/zhoylman/soil-moisture-validation/figs/pearson_comparison.png',
       width = 6.5, height = 8, dpi = 300)



### just uscrn


#compute density data
stie_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station-meta-6-years-min-CDF.csv')

uscrn = stie_meta %>% filter(network == 'USCRN') %$% 
  site_id

merged_uscrn = merged %>%
  filter(site_id %in% uscrn)

density_data_uscrn = merged_uscrn %>%
  group_by(generalized_depth, drought_metric) %>%
  do(density_x = density(.$rmse, bw = 'nrd')$x,
     density_y = density(.$rmse, bw = 'nrd')$y,
     density_bw = density(.$rmse, bw = 'nrd')$bw,
     n = length(.$rmse)) %>%
  unnest(cols = c(density_x, density_y, n, density_bw)) %>%
  mutate(generalized_depth = factor(generalized_depth, levels = c('Shallow (0-4in)', 'Middle (8in - 20in)', 'Deep (>20in)')))  %>%
  mutate(drought_metric = factor(drought_metric, 
                                 levels = c("SPI", "SPEI", "EDDI",
                                            "CPC Soil Moisture", "GRACE Rootzone Soil Moisture",
                                            "SMAP (L4) Rootzone Soil Moisture",
                                            "SPoRT 0-100cm Soil Moisture",
                                            "Topofire Soil Moisture")))

mode_data_uscrn = density_data_uscrn %>%
  group_by(generalized_depth, drought_metric) %>%
  filter(density_y == max(density_y)) %>%
  mutate(modal_r = plyr::round_any(density_x, 0.01))

density_plot = ggplot(data = density_data_uscrn)+
  geom_segment(data = mode_data_uscrn, aes(x = density_x, xend = density_x, y = density_y, yend = 0, color = drought_metric), alpha = 0.6)+
  geom_line(aes(x = density_x, y = density_y, color = drought_metric))+
  ggtitle(paste0('Distrobution of Station Specific Correlation (r)'),'USCRN')+
  theme_bw(base_size = 14)+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  geom_point(data = mode_data_uscrn, aes(x = density_x, y = density_y, color = drought_metric))+
  facet_grid(generalized_depth~.)+
  labs(x = 'Correlation (r)', y = 'Density')+
  scale_colour_manual(values = viridis::turbo(8), name = '')+
  xlim(-0.5,1)+
  theme(legend.position="bottom", legend.box = "horizontal",
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        legend.background = element_rect(size = 0.2))+
  guides(color=guide_legend(nrow=4,byrow=TRUE))

density_plot

ggsave(density_plot, file = '/home/zhoylman/soil-moisture-validation/figs/pearson_comparison_uscrn.png',
       width = 6.5, height = 8, dpi = 300)


