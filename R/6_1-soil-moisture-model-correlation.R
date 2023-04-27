library(tidyverse)
library(sf)
library(magrittr)
library(grid)

`%notin%` = Negate(`%in%`)

stie_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station-meta-6-years-min-CDF.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs('EPSG:4326'))

standardized_soil_moisture_obs = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide-6-years-min-CDF.csv')

standardized_model = list.files('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models', 
                                full.names = T, pattern = 'percentile') %>%
  purrr::map(., read_csv) %>%
  bind_rows() %>%
  mutate(date = time,
         site_id = name,
         model_drought_anomaly = drought_anomaly,
         model_storage_anomaly = storage_anomoly)%>%
  select(nc_id, date, site_id, model_drought_anomaly)

#read in this way becasue it multiprocesses compared to RDS
spi_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spi-data-long-10s.csv') %>%
  filter(timescale %in% c(20,60,90))  %>%
  rename(model_drought_anomaly = spi,
         date = time) %>%
  mutate(nc_id = 'Optimized SPI') %>%
  select(nc_id,date,site_id,timescale,model_drought_anomaly)

spei_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spei-data-long-10s.csv') %>%
  filter(timescale %in% c(20,60,70))  %>%
  rename(model_drought_anomaly = spei,
         date = time) %>%
  mutate(nc_id = 'Optimized SPEI') %>%
  select(nc_id,date,site_id,timescale,model_drought_anomaly)

eddi_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/eddi-data-long-10s.csv') %>%
  filter(timescale %in% c(10,50))  %>%
  rename(model_drought_anomaly = eddi,
         date = time) %>%
  mutate(nc_id = 'Optimized EDDI') %>%
  select(nc_id,date,site_id,timescale,model_drought_anomaly)

#break out data by depth (optimal timescale)
#spi
spi_shallow = spi_all %>%
  filter(timescale == 20)
spi_mid = spi_all %>%
  filter(timescale == 60)
spi_deep = spi_all %>%
  filter(timescale == 90)
#spei
spei_shallow = spei_all %>%
  filter(timescale == 20)
spei_mid = spei_all %>%
  filter(timescale == 60)
spei_deep = spei_all %>%
  filter(timescale == 70)
#eddi
eddi_shallow = eddi_all %>%
  filter(timescale == 10)
eddi_mid = eddi_all %>%
  filter(timescale == 50)
eddi_deep = eddi_all %>%
  filter(timescale == 50)

#bind it all together, starting with shallow spi as a placeholder
binded_raw = bind_rows(list(
  spi_shallow, spi_mid, spi_deep, 
  spei_shallow, spei_mid, spei_deep, 
  eddi_shallow, eddi_mid, eddi_deep,
  standardized_model)) %>%
  left_join(., standardized_soil_moisture_obs, by = c('site_id', 'date'))%>%
  mutate(nc_id = ifelse(nc_id == 'cpc_soil_moisture_percentile', 'CPC Soil Moisture', 
                        ifelse(nc_id == 'grace_rtzn_soil_moisture', 'GRACE Rootzone Soil Moisture', 
                               ifelse(nc_id == 'SMAP_rootzone_soil_moisture', 'SMAP (L4) Rootzone Soil Moisture', 
                                      ifelse(nc_id == 'SPoRT_mean_soil_moisture_0-100cm', 'SPoRT 0-100cm Soil Moisture', 
                                             ifelse(nc_id == 'topofire_soil_moisture', 'Topofire Soil Moisture', 
                                                    ifelse(nc_id == 'Optimized SPI', 'Optimized SPI', 
                                                           ifelse(nc_id == 'Optimized SPEI', 'Optimized SPEI', 
                                                                  ifelse(nc_id == 'Optimized EDDI', 'Optimized EDDI', NA)))))))))

names = c('all','uscrn', 'mt_mesonet', 'scan', 'snotel')
full_names = c('All Sites','USCRN', 'MT Mesonet', 'SCAN', 'SNOTEL')

# define sites in each network
uscrn = stie_meta %>% filter(network == 'USCRN') %$% 
  site_id
mt_mesonet = stie_meta %>% filter(network == 'MT Mesonet') %$% 
  site_id
scan = stie_meta %>% filter(network == 'SCAN') %$% 
  site_id
snotel = stie_meta %>% filter(network %in% c('SNTL', 'SNTLT')) %$% 
  site_id

data_list = list(unique(standardized_soil_moisture_obs$site_id),uscrn, mt_mesonet, scan, snotel)

`%notin%` = Negate(`%in%`)

binded_joined = binded_raw %>%
  filter(lubridate::month(date) >= 5 & lubridate::month(date) <= 10 ) %>% 
  select(nc_id, site_id, date, timescale, model_drought_anomaly, drought_anomaly_4in, drought_anomaly_8in, 
         drought_anomaly_20in,drought_anomaly_36in, drought_anomaly_40in) %>%
  pivot_longer(cols = -c(date, nc_id, site_id, timescale, model_drought_anomaly)) %>%
  #add dummy timescale id for soil moisture models
  mutate(timescale = ifelse(is.na(timescale), 'N/A', timescale)) %>%
  #drop_na() %>%
  mutate(depth = parse_number(name),
         depth = ifelse(depth <= 4, 'Shallow (0-4in)', 
                        ifelse(depth > 4 & depth <=20, 'Middle (8-20in)',
                               ifelse(depth > 20, 'Deep (>20in)', NA))),
         depth = factor(depth, levels = c('Shallow (0-4in)', 'Middle (8-20in)', 'Deep (>20in)'))) 

#select optimal timescales by depth
drought_indicies = binded_joined %>% 
  filter(nc_id %in% c("Optimized SPI", "Optimized SPEI", "Optimized EDDI")) %>%
  filter((nc_id == "Optimized SPI" & timescale == 20 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized SPI" & timescale == 60 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized SPI" & timescale == 90 & depth == 'Deep (>20in)')|
           (nc_id == "Optimized SPEI" & timescale == 20 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized SPEI" & timescale == 60 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized SPEI" & timescale == 70 & depth == 'Deep (>20in)')|
           (nc_id == "Optimized EDDI" & timescale == 10 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized EDDI" & timescale == 50 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized EDDI" & timescale == 50 & depth == 'Deep (>20in)'))

binded = binded_joined %>%
  filter(nc_id %notin% c("Optimized SPI", "Optimized SPEI", "Optimized EDDI")) %>%
  bind_rows(., drought_indicies %>% distinct()) %>%
  #convert to standard anomoly
  mutate(model_drought_anomaly = ifelse(nc_id == 'Optimized SPI'|
                                          nc_id == 'Optimized SPEI'|
                                          nc_id == 'Optimized EDDI', model_drought_anomaly, qnorm(model_drought_anomaly/100)),
         value = qnorm(value),
         #model_drought_anomaly = value,
         #-2 to 2 values only!
         value = ifelse(value < -2, -2, value),
         value = ifelse(value > 2, 2, value),
         model_drought_anomaly = ifelse(model_drought_anomaly < -2, -2, model_drought_anomaly),
         model_drought_anomaly = ifelse(model_drought_anomaly > 2, 2, model_drought_anomaly))


#define color pallette for plotting
color_scale = khroma::color("roma") # color blind safe!

for(i in 1:2){
  print(i)  
  #filter for network of interest
  binded_filtered = binded %>%
    filter(site_id %in% data_list[[i]]) %>%
    drop_na(value, model_drought_anomaly) %>%
    mutate(value = ifelse(nc_id == 'Optimized EDDI', -1*value, value)) %>%
    mutate(nc_id = factor(nc_id, levels = c("Optimized SPI", "Optimized SPEI", "Optimized EDDI",
                                            "CPC Soil Moisture", "GRACE Rootzone Soil Moisture",
                                            "SMAP (L4) Rootzone Soil Moisture", "SPoRT 0-100cm Soil Moisture",
                                            "Topofire Soil Moisture")))

    
  #compute correlation stats for correlation
  drought_anomoly_stats = binded_filtered %>%
    group_by(nc_id, depth) %>%
    do(r = cor(.$model_drought_anomaly, .$value),
       rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
       n = length(.$value)) %>%
    unnest(c(r, n, rmse))
  
  #compute some meta (n sites and obs)
  n_sites = length(unique(binded_filtered$site_id))
  n_obs = binded_filtered %>%
    select(site_id, date, value, name) %>%
    drop_na(value) %>%
    distinct() %$%
    value %>%
    length()
  
  #for plot development
  # binded_filtered_test=binded_filtered %>%
  #   sample_n(20000)
  
  #plot the results!
  plot1 = binded_filtered %>%
    ggplot(., aes(x = value, y = model_drought_anomaly)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    )+
    geom_smooth(method = 'lm', color = 'black', size = 0.5, se  = F)+
    scale_fill_gradientn(colours = color_scale(100) , name = 'Density', guide = "colourbar", limits = c(0,.2), na.value = color_scale(100)[100]) +
    geom_text(data = drought_anomoly_stats, aes(x = -1.75, y = 1.75, label = paste0("RMSE = ", round(rmse, 3))), hjust = 0, fontface = "bold", color = 'white')+
    geom_text(data = drought_anomoly_stats, aes(x = -1.75, y = 1.5, label = paste0("n =", n %>% format(., format="d", big.mark=","))), hjust = 0, fontface = "bold", color = 'white')+
    geom_text(data = drought_anomoly_stats, aes(x = -1.75, y = 1.25, label = paste0("r = ", round(r, 3))), hjust = 0, fontface = "bold", color = 'white')+
    theme_bw(base_size = 15)+
    geom_abline(slope=1, intercept=0, color = 'black', linetype = 'dashed')+
    ylim(c(-2,2))+
    xlim(c(-2,2))+
    facet_grid(depth~nc_id)+
    labs(x = bquote(Drought~Metric~or~Modelled~Soil~Moisture~Drought~Index~(SMDI[mod])), y = bquote(Observed~Soil~Moisture~Drought~Index~(SMDI[obs])))+
    theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
          legend.position = 'bottom', legend.key.width=unit(4,"cm"))+
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
    ggtitle(paste0(full_names[i], ' (May - Oct)'),
    #ggtitle(paste0(full_names[i]), 
            paste0('n (sites) = ', n_sites, ', n (unique soil moisture observations) = ', 
                   n_obs %>% format(.,format="d", big.mark=",")))+
    theme(plot.title = element_text(hjust = 0.5, size=32, margin=margin(0,0,5,0)), 
          plot.subtitle = element_text(hjust = 0.5, margin=margin(0,0,25,0)))
  #plot1
  
  # png(paste0("/home/zhoylman/temp/drought_anomoly_model_comparison_6_year_min_summer_clamped_precomputed_percentiles_", names[i],".png"),
  #     width = 23, height = 13, units = 'in', res = 200)
  # print(plot1)
  # dev.off()
  
  png(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_precomputed_percentiles_", names[i],".png"),
      width = 23, height = 13, units = 'in', res = 200)
  print(plot1)
  dev.off()
}

##############################################################################

# compute site specific RMSE

#define function to compute RMSE from lm object
RMSE = function(lm){
  #Residual sum of squares:
  RSS = c(crossprod(lm$residuals))
  #Mean squared error:
  MSE = RSS / length(lm$residuals)
  #Root MSE:
  RMSE = sqrt(MSE)
  return(RMSE)
}

states = read_sf('/home/zhoylman/soil-moisture-validation-data/raw/shp/conus_states.shp')

site_specific_results = binded %>%
  group_by(nc_id, site_id, depth) %>%
  do(rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
     linearFit = lm(.$value ~ .$model_drought_anomaly),
     pearson = cor(.$model_drought_anomaly, .$value,method="pearson")) %>%
  mutate(linearFit_RMSE = RMSE(linearFit),
         pearson_r = unlist(pearson),
         rmse = unlist(rmse))

write_csv(site_specific_results %>% select(-c('linearFit', 'pearson')), '/home/zhoylman/soil-moisture-validation-data/processed/rmse-comparison/soil-moisture-model-r.csv')

site_specific_results_spatial = left_join(stie_meta, site_specific_results, by = 'site_id')

spatial_plot = ggplot(site_specific_results_spatial %>% filter(!is.na(rmse)))+
  facet_grid(nc_id ~ depth)+
  geom_sf(data = states, fill = 'transparent')+
  geom_sf(aes(fill = rmse),color = 'black', shape = 21, size = 2)+
  scale_fill_gradientn(colours = color_scale(100) , name = 'RMSE', guide = "colourbar",
                       breaks=c(min(site_specific_results_spatial$rmse, na.rm = T),
                                mean(c(max(site_specific_results_spatial$rmse, na.rm = T),min(site_specific_results_spatial$rmse, na.rm = T))),
                                max(max(site_specific_results_spatial$rmse, na.rm = T))), labels=c(paste0("< ", min(site_specific_results_spatial$rmse, na.rm = T) %>% round(., 1)),
                                                                                           mean(c(max(site_specific_results_spatial$rmse, na.rm = T),min(site_specific_results_spatial$rmse, na.rm = T))) %>% round(., 1),
                                                                                           paste0("> ", max(site_specific_results_spatial$rmse, na.rm = T) %>% round(., 1))))+
  theme_bw(base_size = 12)+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
        legend.position = 'bottom', legend.key.width=unit(4,"cm"))+
  guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))
  
  

png(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_spatial.png"),
    width = 13, height = 13, units = 'in', res = 200)
print(spatial_plot)
dev.off()

write_sf(site_specific_results_spatial, '/home/zhoylman/temp/spatail_error.geojson')

# drought_anomoly_stats = binded %>%
#   group_by(nc_id, depth) %>%
#   do(r = cor(.$model_drought_anomaly, .$value),
#      r_classes = cor(.$modelled_class, .$obs_class),
#      #lm = lm(.$value ~ .$model_drought_anomaly),
#      rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
#      rmse_classes = sqrt(mean((.$obs_class - .$modelled_class)^2)),
#      n = length(.$value)) %>%
#   unnest(c(r, n, rmse, r_classes, rmse_classes))

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





# 
# 
# binded_storage = binded_raw %>%
#   filter(site_id %in% uscrn) %>%
#   select(nc_id, site_id, model_storage_anomaly, storage_anomaly_4in, storage_anomaly_8in, 
#          storage_anomaly_20in,storage_anomaly_36in, storage_anomaly_40in) %>%
#   pivot_longer(cols = -c(nc_id, site_id, model_storage_anomaly)) %>%
#   drop_na() %>%
#   mutate(depth = paste0(parse_number(name), 'in'),
#          depth = factor(depth, levels = c('4in', '8in', '20in', '40in')))
# 
# 
# plot2 = binded_storage %>%
#   #sample_n(., 100000) %>%
#   ggplot(., aes(x = value, y = model_storage_anomaly)) +
#   geom_bin2d(bins = 70) +
#   #geom_point()+
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()+
#   geom_abline(slope=1, intercept=0, color = 'green')+
#   #geom_smooth(method = 'lm', color = 'red')+
#   ylim(c(-6,6))+
#   xlim(c(-6,6))+
#   facet_grid(depth~nc_id)+
#   labs(x = paste0('Modeled Storage Anomoly'), y = 'Observed Storage Anomoly')+
#   theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill='transparent') ) 
# 
# plot2
