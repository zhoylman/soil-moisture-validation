library(tidyverse)

standardized_soil_moisture_obs = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide-6-years-min.csv')

standardized_model = list.files('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models', 
                                full.names = T, pattern = 'storage') %>%
  purrr::map(., read_csv) %>%
  bind_rows() %>%
  mutate(date = time,
         site_id = name,
         model_drought_anomaly = drought_anomaly,
         model_storage_anomaly = storage_anomoly)%>%
  select(nc_id, date, site_id, model_drought_anomaly, model_storage_anomaly)

binded_raw = left_join(standardized_soil_moisture_obs, standardized_model, by = c('site_id', 'date')) %>%
  mutate(nc_id = ifelse(nc_id == 'cpc_soil_moisture', 'CPC Soil Moisture', 
                        ifelse(nc_id == 'grace_rtzn_soil_moisture', 'GRACE Rootzone Soil Moisture', 
                               ifelse(nc_id == 'SMAP_subsurface_soil_moisture', 'SMAP Subsurface Soil Moisture', 
                                      ifelse(nc_id == 'SPoRT_mean_soil_moisture_0-100cm', 'SPoRT 0-100cm Soil Moisture', 
                                             ifelse(nc_id == 'topofire_soil_moisture', 'Topofire Soil Moisture', NA))))))

binded_raw[binded_raw == Inf] = NA
binded_raw[binded_raw == -Inf] = NA

names = c('all','uscrn', 'mt_mesonet', 'scan', 'snotel')
full_names = c('All Sites','USCRN', 'MT Mesonet', 'SCAN', 'SNOTEL / SNOTEL LITE')

uscrn = unique(standardized_soil_moisture_obs$site_id)[680:791]
mt_mesonet = unique(standardized_soil_moisture_obs$site_id)[1:78]
scan = unique(standardized_soil_moisture_obs$site_id)[79:189]
snotel = unique(standardized_soil_moisture_obs$site_id)[190:679]

data_list = list(unique(standardized_soil_moisture_obs$site_id),uscrn, mt_mesonet, scan, snotel)

`%notin%` = Negate(`%in%`)

for(i in 1:length(names)){
  print(i)  
  binded = binded_raw %>%
    filter(site_id %in% data_list[[i]]) %>%
    filter(lubridate::month(date) >= 5 & lubridate::month(date) <= 10 ) %>% 
    select(nc_id, site_id, model_drought_anomaly, drought_anomaly_4in, drought_anomaly_8in, 
           drought_anomaly_20in,drought_anomaly_36in, drought_anomaly_40in) %>%
    pivot_longer(cols = -c(nc_id, site_id, model_drought_anomaly)) %>%
    drop_na() %>%
    mutate(depth = parse_number(name),
           depth = ifelse(depth <= 4, 'Shallow (0-4in)', 
                          ifelse(depth > 4 & depth <=20, 'Middle (8-20in)',
                                 ifelse(depth > 20, 'Deep (>20in)', NA))),
           depth = factor(depth, levels = c('Shallow (0-4in)', 'Middle (8-20in)', 'Deep (>20in)'))) %>%
    mutate(obs_class = .bincode(value, breaks = c(-Inf, -2.0, -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2, Inf)),
           modelled_class = .bincode(model_drought_anomaly, breaks = c(-Inf, -2.0, -1.6, -1.3, -0.8, -0.5, 0.5, 0.8, 1.3, 1.6, 2, Inf))) 
  
  binded_clamped = binded %>%
    mutate(value = ifelse(value > 2, 2, value), value = ifelse(value < -2, -2, value),
            model_drought_anomaly = ifelse(model_drought_anomaly > 2, 2, model_drought_anomaly),
            model_drought_anomaly = ifelse(model_drought_anomaly < -2, -2, model_drought_anomaly))
  
  drought_anomoly_stats = binded %>%
    group_by(nc_id, depth) %>%
    do(r = cor(.$model_drought_anomaly, .$value),
       r_classes = cor(.$modelled_class, .$obs_class),
       #lm = lm(.$value ~ .$model_drought_anomaly),
       rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
       rmse_classes = sqrt(mean((.$obs_class - .$modelled_class)^2)),
       n = length(.$value)) %>%
    unnest(c(r, n, rmse, r_classes, rmse_classes))
  
  drought_anomoly_stats_clamped = binded_clamped %>%
    group_by(nc_id, depth) %>%
    do(r = cor(.$model_drought_anomaly, .$value),
       r_classes = cor(.$modelled_class, .$obs_class),
       #lm = lm(.$value ~ .$model_drought_anomaly),
       rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
       rmse_classes = sqrt(mean((.$obs_class - .$modelled_class)^2)),
       n = length(.$value)) %>%
    unnest(c(r, n, rmse, r_classes, rmse_classes))
  
  color_scale = khroma::color('roma')
  
  n_sites = length(unique(binded$site_id))
  n_obs = length(binded$value)
  
  # plot1 = binded %>%
  #   #sample_n(., 100000) %>%
  #   ggplot(., aes(x = value, y = model_drought_anomaly)) +
  #   geom_bin2d(bins = 100, aes(fill = after_stat(density))) +
  #   geom_smooth(method = 'lm', color = 'black', size = 0.5)+
  #   scale_fill_gradientn(colours = color_scale(100) , name = 'Density', guide = "colourbar") +
  #   geom_text(data = drought_anomoly_stats, aes(x = -6, y = 6, label = paste0("RMSE = ", round(rmse, 3))), hjust = 0, fontface = "bold")+
  #   geom_text(data = drought_anomoly_stats, aes(x = -6, y = 5, label = paste0("n =", n %>% format(., format="d", big.mark=","))), hjust = 0, fontface = "bold")+
  #   geom_text(data = drought_anomoly_stats, aes(x = -6, y = 4, label = paste0("r = ", round(r, 3))), hjust = 0, fontface = "bold")+
  #   theme_bw(base_size = 16)+
  #   geom_abline(slope=1, intercept=0, color = 'black', linetype = 'dashed')+
  #   ylim(c(-6,6))+
  #   xlim(c(-6,6))+
  #   facet_grid(depth~nc_id)+
  #   labs(x = 'Modelled Drought Anomoly', y = 'Observed Drought Anomoly')+
  #   theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
  #         legend.position = 'bottom', legend.key.width=unit(4,"cm"))+
  #   guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
  #   #ggtitle(paste0(full_names[i], ' (May - Oct)'))+
  #   ggtitle(paste0(full_names[i]), 
  #           paste0('n (sites) = ', n_sites, ', n (observations) = ', 
  #                  n_obs %>% format(.,format="d", big.mark=",")))+
  #   theme(plot.title = element_text(hjust = 0.5, size=32), plot.subtitle = element_text(hjust = 0.5))
  # #plot1
  # 
  # png(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_all_time_", names[i],".png"),
  #     width = 17, height = 13, units = 'in', res = 200)
  # print(plot1)
  # dev.off()
  # 
  
  plot1 = binded_clamped %>%
    #sample_n(., 100000) %>%
    ggplot(., aes(x = value, y = model_drought_anomaly)) +
    #geom_bin2d(bins = 25, aes(fill = after_stat(density))) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    )+
    geom_smooth(method = 'lm', color = 'black', size = 0.5)+
    scale_fill_gradientn(colours = color_scale(100) , name = 'Density', guide = "colourbar") +
    geom_text(data = drought_anomoly_stats_clamped, aes(x = -1.75, y = 1.75, label = paste0("RMSE = ", round(rmse, 3))), hjust = 0, fontface = "bold", color = 'white')+
    geom_text(data = drought_anomoly_stats_clamped, aes(x = -1.75, y = 1.5, label = paste0("n =", n %>% format(., format="d", big.mark=","))), hjust = 0, fontface = "bold", color = 'white')+
    geom_text(data = drought_anomoly_stats_clamped, aes(x = -1.75, y = 1.25, label = paste0("r = ", round(r, 3))), hjust = 0, fontface = "bold", color = 'white')+
    theme_bw(base_size = 16)+
    geom_abline(slope=1, intercept=0, color = 'black', linetype = 'dashed')+
    ylim(c(-2,2))+
    xlim(c(-2,2))+
    facet_grid(depth~nc_id)+
    labs(x = 'Modelled Drought Anomoly', y = 'Observed Drought Anomoly')+
    theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
          legend.position = 'bottom', legend.key.width=unit(4,"cm"))+
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
    ggtitle(paste0(full_names[i], ' (May - Oct)'),
    #ggtitle(paste0(full_names[i]), 
            paste0('n (sites) = ', n_sites, ', n (observations) = ', 
                   n_obs %>% format(.,format="d", big.mark=",")))+
    theme(plot.title = element_text(hjust = 0.5, size=32), plot.subtitle = element_text(hjust = 0.5))
  #plot1
  
  png(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_", names[i],".png"),
      width = 17, height = 13, units = 'in', res = 200)
  print(plot1)
  dev.off()
}


##############################################################################


binded_storage = binded_raw %>%
  filter(site_id %in% uscrn) %>%
  select(nc_id, site_id, model_storage_anomaly, storage_anomaly_4in, storage_anomaly_8in, 
         storage_anomaly_20in,storage_anomaly_36in, storage_anomaly_40in) %>%
  pivot_longer(cols = -c(nc_id, site_id, model_storage_anomaly)) %>%
  drop_na() %>%
  mutate(depth = paste0(parse_number(name), 'in'),
         depth = factor(depth, levels = c('4in', '8in', '20in', '40in')))


plot2 = binded_storage %>%
  #sample_n(., 100000) %>%
  ggplot(., aes(x = value, y = model_storage_anomaly)) +
  geom_bin2d(bins = 70) +
  #geom_point()+
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  geom_abline(slope=1, intercept=0, color = 'green')+
  #geom_smooth(method = 'lm', color = 'red')+
  ylim(c(-6,6))+
  xlim(c(-6,6))+
  facet_grid(depth~nc_id)+
  labs(x = paste0('Modeled Storage Anomoly'), y = 'Observed Storage Anomoly')+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill='transparent') ) 

plot2
