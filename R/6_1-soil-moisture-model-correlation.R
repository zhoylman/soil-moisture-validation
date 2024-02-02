library(tidyverse)
library(sf)
library(magrittr)
library(shadowtext)

plot_all_models = F

#define special 
`%notin%` = Negate(`%in%`)

#import site metadata
stie_meta = read_csv('~/soil-moisture-validation-data/processed/standardized-soil-moisture/beta-standardized-station-meta-6-years-min-CDF-w-mean.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs('EPSG:4326'))

#import observed soil moisture data
standardized_soil_moisture_obs = read_csv('~/soil-moisture-validation-data/processed/standardized-soil-moisture/beta-standardized-soil-moisture-data-wide-6-years-min-CDF-w-mean.csv')

#import soil moisture models 
standardized_model = list.files('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture-models', 
                                full.names = T, pattern = 'percentile') %>%
  purrr::map(., read_csv) %>%
  bind_rows() %>%
  mutate(date = time,
         site_id = name,
         model_drought_anomaly = drought_anomaly,
         model_storage_anomaly = storage_anomoly)%>%
  dplyr::select(nc_id, date, site_id, model_drought_anomaly)

#import all SPI data for relevant time scales from global maxima
spi_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spi-data-long-10s.csv') %>%
  filter(timescale %in% c(30,20,60,80))  %>%
  rename(model_drought_anomaly = spi,
         date = time) %>%
  mutate(nc_id = 'Optimized SPI') %>%
  dplyr::select(nc_id,date,site_id,timescale,model_drought_anomaly)

#import all SPEI data for relevant time scales from global maxima
spei_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spei-data-long-10s.csv') %>%
  filter(timescale %in% c(30,20,50,70))  %>%
  rename(model_drought_anomaly = spei,
         date = time) %>%
  mutate(nc_id = 'Optimized SPEI') %>%
  dplyr::select(nc_id,date,site_id,timescale,model_drought_anomaly)

#import all EDDI data for relevant time scales from global maxima
eddi_all = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/eddi-data-long-10s.csv') %>%
  filter(timescale %in% c(20, 10, 40, 60))  %>%
  rename(model_drought_anomaly = eddi,
         date = time) %>%
  mutate(nc_id = 'Optimized EDDI') %>%
  dplyr::select(nc_id,date,site_id,timescale,model_drought_anomaly)

#break out data by depth (optimal timescale) 
#this discreate approach will allow us to replace values for specific depths
#spi
spi_mean = spi_all %>%
  filter(timescale == 30)
spi_shallow = spi_all %>%
  filter(timescale == 20)
spi_mid = spi_all %>%
  filter(timescale == 60)
spi_deep = spi_all %>%
  filter(timescale == 80)
#spei
spei_mean = spei_all %>%
  filter(timescale == 30)
spei_shallow = spei_all %>%
  filter(timescale == 20)
spei_mid = spei_all %>%
  filter(timescale == 50)
spei_deep = spei_all %>%
  filter(timescale == 70)
#eddi
eddi_mean = eddi_all %>%
  filter(timescale == 20)
eddi_shallow = eddi_all %>%
  filter(timescale == 10)
eddi_mid = eddi_all %>%
  filter(timescale == 40)
eddi_deep = eddi_all %>%
  filter(timescale == 60)

#bind it all together, drought metrics, soil moisture models and observed values
binded_raw = bind_rows(list(
  spi_mean, spi_shallow, spi_mid, spi_deep, 
  spei_mean, spei_shallow, spei_mid, spei_deep, 
  eddi_mean, eddi_shallow, eddi_mid, eddi_deep,
  standardized_model)) %>%
  left_join(., standardized_soil_moisture_obs, by = c('site_id', 'date'))%>%
  mutate(nc_id = ifelse(nc_id == 'cpc_soil_moisture_percentile', 'CPC Soil Moisture', 
                        ifelse(nc_id == 'grace_rtzn_soil_moisture', 'GRACE Rootzone Soil Moisture', 
                               ifelse(nc_id == 'SMAP_rootzone_soil_moisture', 'SMAP (L4) Rootzone Soil Moisture', 
                                      ifelse(nc_id == 'SPoRT_mean_soil_moisture_0-100cm', 'SPoRT 0-100cm Soil Moisture', 
                                             ifelse(nc_id == 'topofire_soil_moisture', 'Topofire Soil Moisture', 
                                                    ifelse(nc_id == 'NLDAS2_ensamble_soil_moisture_0-100cm', 'NLDAS-2 Ensemble 0-100cm Soil Moisture', 
                                                         ifelse(nc_id == 'NLDAS2_MOSAIC_soil_moisture_0-100cm', 'NLDAS-2 MOSAIC 0-100cm Soil Moisture', 
                                                                ifelse(nc_id == 'NLDAS2_VIC_soil_moisture_0-100cm', 'NLDAS-2 VIC 0-100cm Soil Moisture', 
                                                                       ifelse(nc_id == 'NLDAS2_NOAH_soil_moisture_0-100cm', 'NLDAS-2 NOAH 0-100cm Soil Moisture', 
                                                                            ifelse(nc_id == 'Optimized SPI', 'Optimized SPI', 
                                                                                ifelse(nc_id == 'Optimized SPEI', 'Optimized SPEI', 
                                                                                      ifelse(nc_id == 'Optimized EDDI', 'Optimized EDDI', NA)))))))))))))

#define each network independently incase we want to look at each indivitual network
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
#generate a list of network names for filtering
data_list = list(unique(standardized_soil_moisture_obs$site_id),uscrn, mt_mesonet, scan, snotel)

#clean up final dataset and pivot to long
binded_joined = binded_raw %>%
  #filter by date (summertime)
  filter(lubridate::month(date) >= 5 & lubridate::month(date) <= 10 ) %>% 
  #dplyr::select the relevant variables and drought anomolies
  dplyr::select(nc_id, site_id, date, timescale, model_drought_anomaly, drought_anomaly_4in, drought_anomaly_8in, 
         drought_anomaly_20in,drought_anomaly_36in, drought_anomaly_40in, drought_anomaly_mean) %>%
  #pivot to longer form
  pivot_longer(cols = -c(date, nc_id, site_id, timescale, model_drought_anomaly)) %>%
  #add dummy timescale id for soil moisture models
  mutate(timescale = ifelse(is.na(timescale), 'N/A', timescale)) %>%
  #mutate, add depth, modify names
  mutate(depth = parse_number(name),
         depth = ifelse(depth <= 4, 'Shallow (0-4in)', 
                        ifelse(depth > 4 & depth <=20, 'Middle (8-20in)',
                               ifelse(depth > 20, 'Deep (>20in)', NA))),
         depth = ifelse(name == 'drought_anomaly_mean', 'Depth Averaged', depth),
         depth = factor(depth, levels = c('Depth Averaged', 'Shallow (0-4in)', 'Middle (8-20in)', 'Deep (>20in)'))) 

#now that we have the drought metric data formatted correctly we will 
#dplyr::select the drought metrics and dplyr::select optimal timescales by depth 
drought_indicies = binded_joined %>% 
  filter(nc_id %in% c("Optimized SPI", "Optimized SPEI", "Optimized EDDI")) %>%
  filter((nc_id == "Optimized SPI" & timescale == 30 & depth == 'Depth Averaged')|
           (nc_id == "Optimized SPI" & timescale == 20 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized SPI" & timescale == 60 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized SPI" & timescale == 80 & depth == 'Deep (>20in)')|
           (nc_id == "Optimized SPEI" & timescale == 30 & depth == 'Depth Averaged')|
           (nc_id == "Optimized SPEI" & timescale == 20 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized SPEI" & timescale == 50 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized SPEI" & timescale == 70 & depth == 'Deep (>20in)')|
           (nc_id == "Optimized EDDI" & timescale == 20 & depth == 'Depth Averaged')|
           (nc_id == "Optimized EDDI" & timescale == 10 & depth == 'Shallow (0-4in)')|
           (nc_id == "Optimized EDDI" & timescale == 40 & depth == 'Middle (8-20in)')|
           (nc_id == "Optimized EDDI" & timescale == 60 & depth == 'Deep (>20in)'))

#now that we have the correct drought metric data, lets rebind it all together
#and make final conversions to standard anomaly 
binded = binded_joined %>%
  #remove old drought index data
  filter(nc_id %notin% c("Optimized SPI", "Optimized SPEI", "Optimized EDDI")) %>%
  #bind rows with the new filtered drought index data (optimal timescales)
  bind_rows(., drought_indicies %>% distinct()) %>%
  #convert to soil moisture models to standard anomaly
  #drought metrics are already in standard anomaly form
  mutate(model_drought_anomaly = ifelse(nc_id == 'Optimized SPI'|
                                          nc_id == 'Optimized SPEI'|
                                          nc_id == 'Optimized EDDI', model_drought_anomaly, qnorm(model_drought_anomaly/100)),
         value = qnorm(value),
         #clamp values to -2 to 2 values only!
         value = ifelse(value < -2, -2, value),
         value = ifelse(value > 2, 2, value),
         model_drought_anomaly = ifelse(model_drought_anomaly < -2, -2, model_drought_anomaly),
         model_drought_anomaly = ifelse(model_drought_anomaly > 2, 2, model_drought_anomaly))

#define color pallette for plotting
color_scale = khroma::color("roma") # color blind safe!

#loop through all sites and uscrn only
for(i in 1:2){
  print(i)  
  #filter for network of interest
  binded_filtered = binded %>%
    #filter for all or dplyr::select networks
    filter(site_id %in% data_list[[i]]) %>%
    #drop rows where we are missing either x or y
    drop_na(value, model_drought_anomaly) %>%
    #invert eddi
    mutate(value = ifelse(nc_id == 'Optimized EDDI', -1*value, value)) %>%
    #define factor levels for plotting
    mutate(., nc_id = factor(nc_id, levels = c("Optimized SPI", "Optimized SPEI", "Optimized EDDI",
                                                                      "CPC Soil Moisture", "GRACE Rootzone Soil Moisture",
                                                                      "NLDAS-2 VIC 0-100cm Soil Moisture",
                                                                      "NLDAS-2 NOAH 0-100cm Soil Moisture",
                                                                      "NLDAS-2 MOSAIC 0-100cm Soil Moisture",
                                                                      "NLDAS-2 Ensemble 0-100cm Soil Moisture",
                                                                      "SPoRT 0-100cm Soil Moisture",
                                                                      "SMAP (L4) Rootzone Soil Moisture", 
                                                                      "Topofire Soil Moisture")))
      
  
  #compute correlation stats for correlation
  drought_anomoly_stats = binded_filtered %>%
    #make sure we have only unique observations
    distinct() %>%
    #group by model and depth
    group_by(nc_id, depth) %>%
    #compute the statistics (r, rmse, and number of obs)
    do(r = cor(.$model_drought_anomaly, .$value),
       p_value = cor.test(.$model_drought_anomaly, .$value)$p.value,
       rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
       n = length(.$value)) %>%
    unnest(c(r, p_value, n, rmse)) %>%
    mutate(p_value = ifelse(p_value < 0.01, '< 0.01', p_value))
  
  #compute some meta (n sites and obs)
  n_sites = length(unique(binded_filtered$site_id))
  #compute number of observations for each group (depth*model)
  n_obs = binded_filtered %>%
    dplyr::select(site_id, date, value, name) %>%
    drop_na(value) %>%
    distinct() %$%
    value %>%
    length()
  
  binded_filtered_test = binded_filtered %>%
    sample_n(., 10000)

  #vertical
  plot2 = binded_filtered %>%
    ggplot(., aes(x = model_drought_anomaly, y = value)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    )+
    scale_fill_gradientn(colours = color_scale(100) , name = 'Density', guide = "colourbar", limits = c(0,.2), na.value = color_scale(100)[100]) +
    geom_shadowtext(data = drought_anomoly_stats, aes(x = -1.75, y = 1.65, label = paste0("RMSE = ", round(rmse, 3))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    geom_shadowtext(data = drought_anomoly_stats, aes(x = -1.75, y = 1.15, label = paste0("n =", n %>% format(., format="d", big.mark=","))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    geom_shadowtext(data = drought_anomoly_stats, aes(x = -1.75, y = 0.65, label = paste0("r = ", round(r, 3))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    theme_bw(base_size = 17)+
    geom_abline(slope=1, intercept=0, color = 'black', linetype = 'dashed')+
    ylim(c(-2,2))+
    xlim(c(-2,2))+
    facet_grid(nc_id~depth, labeller = label_wrap_gen(width=16))+
    labs(x = bquote(Drought~Metric~or~Modelled~Soil~Moisture~Index~(SMI[mod])), y = bquote(Observed~Soil~Moisture~Index~(SMI[obs])))+
    theme(legend.key = element_blank(), strip.background =element_blank(),
          legend.position = 'bottom', legend.key.width=unit(2,"cm"))+
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
    ggtitle(paste0(full_names[i], ' (May - Oct)'),
            paste0('n (sites) = ', n_sites, ', n (unique soil moisture observations) = ', 
                   n_obs %>% format(.,format="d", big.mark=",")))+
    theme(plot.title = element_text(hjust = 0.5, size=18), 
          plot.subtitle = element_text(hjust = 0.5, size=12),
          strip.text.y.right = element_text(angle=360, vjust = 0.5, hjust = 0.5),
          strip.placement = "outside")
  
    if(plot_all_models == T){
      jpeg(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_precomputed_percentiles_verticals_", names[i],".jpg"),
          width = 12, height = 20, units = 'in', res = 600)
      print(plot2)
      dev.off()
    }

    #make table
    # Create a table using kable function
    drought_anomoly_stats_tabble = drought_anomoly_stats %>%
      mutate(r = round(r, 3),
             rmse = round(rmse, 3))%>%
      rename(Model = nc_id,
             Depth = depth,
             `Pearson's r` = r,
             `p-value` = p_value,
             RMSE = rmse) %>%
      arrange(factor(Depth, levels = c('Depth Averaged', 'Shallow (0-4in)', 'Middle (8-20in)', 'Deep (>20in)')))
    
    bold_rmse = drought_anomoly_stats_tabble %>%
      group_by(Depth) %>%
      mutate(bold = ifelse(RMSE == min(RMSE), TRUE, FALSE)) %>%
      ungroup() %>%
      filter(bold == TRUE)
    
    bold_r = drought_anomoly_stats_tabble %>%
      group_by(Depth) %>%
      mutate(bold = ifelse(`Pearson's r` == max(`Pearson's r`), TRUE, FALSE)) %>%
      ungroup() %>%
      filter(bold == TRUE)
    
    bold_index_RMSE =  which(
      paste0(drought_anomoly_stats_tabble$`RMSE`, drought_anomoly_stats_tabble$Depth) 
      %in% 
        paste0(bold_rmse$`RMSE`, bold_rmse$Depth ))
    
    bold_index_r = which(
      paste0(drought_anomoly_stats_tabble$`Pearson's r`, drought_anomoly_stats_tabble$Depth) 
      %in% 
        paste0(bold_r$`Pearson's r`, bold_r$Depth ))

    Sys.setenv("OPENSSL_CONF"="/dev/null")
    
    kableExtra::kbl(drought_anomoly_stats_tabble) %>%
      #kableExtra::kable_styling(latex_options = c("basic")) %>%
      kableExtra::kable_classic_2(full_width = F, html_font = "Cambria", bootstrap_options = c("condensed")) %>%
      kableExtra::row_spec(0, font_size=15) %>%
      kableExtra::row_spec(., bold_index_r,  bold = FALSE, color = 'blue') %>%
      kableExtra::row_spec(., bold_index_RMSE,  bold = TRUE) %>%
      kableExtra::row_spec(., c(12,24, 36),  extra_css = "border-bottom: 1px solid") %>%
      kableExtra::save_kable(file = paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/model_table_", names[i],".jpg"))
    
    #write table in csv form
    write_csv(drought_anomoly_stats_tabble, paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/model_table_", names[i],".csv"))
    
  #only mean soil moisture (plot if all models are included)
  binded_filtered_mean = binded_filtered %>%
    filter(depth == 'Depth Averaged') 
  
  drought_anomoly_stats_mean = drought_anomoly_stats %>%
    filter(depth == 'Depth Averaged')
  
  n_obs_mean = binded_filtered_mean %>%
    dplyr::select(site_id, date, value, name) %>%
    drop_na(value) %>%
    distinct() %$%
    value %>%
    length()
  
  plot3 = binded_filtered_mean %>%
    ggplot(., aes(x = model_drought_anomaly, y = value)) +
    stat_density_2d(
      geom = "raster",
      aes(fill = after_stat(density)),
      contour = FALSE
    )+
    scale_fill_gradientn(colours = color_scale(100) , name = 'Density', guide = "colourbar", limits = c(0,.2), na.value = color_scale(100)[100]) +
    geom_shadowtext(data = drought_anomoly_stats_mean, aes(x = -1.75, y = 1.65, label = paste0("RMSE = ", round(rmse, 3))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    geom_shadowtext(data = drought_anomoly_stats_mean, aes(x = -1.75, y = 1.15, label = paste0("n =", n %>% format(., format="d", big.mark=","))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    geom_shadowtext(data = drought_anomoly_stats_mean, aes(x = -1.75, y = 0.65, label = paste0("r = ", round(r, 3))), hjust = 0, fontface = "bold", color = 'white', size = 4)+
    theme_bw(base_size = 14)+
    geom_abline(slope=1, intercept=0, color = 'black', linetype = 'dashed')+
    ylim(c(-2,2))+
    xlim(c(-2,2))+
    facet_wrap(~nc_id, labeller = label_wrap_gen(width=25), nrow = 3)+
    labs(x = bquote(Drought~Metric~or~Modelled~Soil~Moisture~Index~(SMI[mod])), y = bquote(Observed~Soil~Moisture~Index~(SMI[obs])))+
    theme(legend.key = element_blank(), strip.background =element_blank(),
          legend.position = 'bottom', legend.key.width=unit(2,"cm"))+
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
    ggtitle(paste0(full_names[i], ' (May - Oct)'),
            paste0('n (sites) = ', n_sites, ', n (observations) = ', 
                   n_obs_mean %>% format(.,format="d", big.mark=",")))+
    theme(plot.title = element_text(hjust = 0.5, size=18), 
          plot.subtitle = element_text(hjust = 0.5, size=12),
          strip.text.y.right = element_text(angle=360, vjust = 0.5, hjust = 0.5),
          strip.placement = "outside")
  
    jpeg(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_precomputed_percentiles_mean_", names[i],".jpg"),
        width = 12, height = 9, units = 'in', res = 600)
    print(plot3)
    dev.off()
    
    ## ordinal classification error
    ## classify drought coniditons based on observed soil moisture 
    ## compute error based on ordinal difference in classes
    value_pairs = data.frame(
      old_value = 1:11,
      new_value = c("D4 (Exceptional Drought)", "D3 (Extreme Drought)", "D2 (Severe Drought)", "D1 (Moderate Drought)", "D0 (Abnormally Dry)",
                    "Neutral", 'W0 (Abnormally Wet)', 'W1 (Moderately Wet)', 'W2 (Severely Wet)', 'W3 (Extremely Wet)', 'W4 (Exceptionally Wet)')
    )
    
    ordinal_class_error = binded_filtered %>%
      mutate(obs_class = .bincode(value, breaks = rev(c(Inf,1.9999,1.6,1.3,0.8,0.5, -0.5, -0.8, -1.3, -1.6, -1.9999, -Inf))),
             modelled_class = .bincode(model_drought_anomaly, breaks = rev(c(Inf,1.9999,1.6,1.3,0.8,0.5, -0.5, -0.8, -1.3, -1.6, -1.9999, -Inf)))) %>%
      group_by(nc_id, obs_class, depth) %>%
      do(mse = mean((.$value - .$model_drought_anomaly)^2) %>% round(., 3),
         rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2))%>% round(., 3),
         MAE = mean(abs(.$value - .$model_drought_anomaly)) %>% round(., 3),
         n =  length(.$value)) %>%
      unnest(c('mse', 'rmse', 'MAE', 'n')) %>%
      ungroup() %>%
      group_by(obs_class, depth) %>%
      slice_min(order_by = mse) %>%
      ungroup() %>%
      left_join(value_pairs, by = c("obs_class" = "old_value")) %>%
      mutate(DM_class = coalesce(new_value, as.character(obs_class))) %>%
      dplyr::select(nc_id, depth, DM_class, mse, MAE, n) %>%
      arrange(factor(depth, levels = c('Depth Averaged', 'Shallow (0-4in)', 'Middle (8-20in)', 'Deep (>20in)'))) %>%
      rename(`Best Model` = nc_id,
             Depth = depth,
             `Theoretical Drought Class` = DM_class,
             MSE = mse) %>%
      filter(Depth == "Depth Averaged")
    
    kableExtra::kbl(ordinal_class_error) %>%
      kableExtra::kable_classic_2(full_width = F, html_font = "Cambria", bootstrap_options = c("condensed")) %>%
      kableExtra::row_spec(0, font_size=15) %>%
      kableExtra::row_spec(., c(5,6),  extra_css = "border-bottom: 1px solid") %>%
      kableExtra::save_kable(file = paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/ordinal_model_table_", names[i],".jpg"))
    
    write_csv(ordinal_class_error, paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/ordinal_model_table_", names[i],".csv"))
} 

##############################################################################

# compute site specific RMSE map

states = read_sf('/home/zhoylman/soil-moisture-validation-data/raw/shp/conus_states.shp')

site_specific_results = binded %>%
  mutate(value = ifelse(nc_id == 'Optimized EDDI', -1*value, value)) %>%
  drop_na(value, model_drought_anomaly) %>%
  group_by(nc_id, site_id, depth) %>%
  do(rmse = sqrt(mean((.$value - .$model_drought_anomaly)^2)),
     pearson = cor(.$model_drought_anomaly, .$value,method="pearson")) %>%
  mutate(pearson_r = unlist(pearson),
         rmse = unlist(rmse)) %>%
  mutate(nc_id = factor(nc_id, levels =c("Optimized SPI", "Optimized SPEI", "Optimized EDDI",
                                         "CPC Soil Moisture", "GRACE Rootzone Soil Moisture",
                                         "NLDAS-2 VIC 0-100cm Soil Moisture",
                                         "NLDAS-2 NOAH 0-100cm Soil Moisture",
                                         "NLDAS-2 MOSAIC 0-100cm Soil Moisture",
                                         "NLDAS-2 Ensemble 0-100cm Soil Moisture",
                                         "SPoRT 0-100cm Soil Moisture",
                                         "SMAP (L4) Rootzone Soil Moisture", 
                                         "Topofire Soil Moisture")))

#write_csv(site_specific_results %>% select(-c('linearFit', 'pearson')), '/home/zhoylman/soil-moisture-validation-data/processed/rmse-comparison/soil-moisture-model-r.csv')

site_specific_results_spatial = left_join(stie_meta, site_specific_results, by = 'site_id') %>%
  mutate(rmse = ifelse(rmse > quantile(rmse,0.99, na.rm = T), quantile(rmse,0.99, na.rm = T), rmse),
         rmse = ifelse(rmse < quantile(rmse,0.01, na.rm = T), quantile(rmse,0.01, na.rm = T), rmse),
         pearson_r = ifelse(pearson_r < quantile(pearson_r,0.01, na.rm = T), quantile(pearson_r,0.01, na.rm = T), pearson_r),
         pearson_r = ifelse(pearson_r > quantile(pearson_r,0.99, na.rm = T), quantile(pearson_r,0.99, na.rm = T), pearson_r))

#rmse
spatial_plot = ggplot(site_specific_results_spatial %>% filter(!is.na(rmse)))+
  facet_grid(nc_id ~ depth, labeller = label_wrap_gen(width=10))+
  geom_sf(data = states, fill = 'transparent')+
  geom_sf(aes(fill = rmse),color = 'black', shape = 21, size = 1)+
  scale_fill_gradientn(colours = color_scale(100) , name = 'RMSE', guide = "colourbar",
                       breaks=c(min(site_specific_results_spatial$rmse, na.rm = T),
                                mean(c(max(site_specific_results_spatial$rmse, na.rm = T),min(site_specific_results_spatial$rmse, na.rm = T))),
                                max(max(site_specific_results_spatial$rmse, na.rm = T))), labels=c(paste0("< ", min(site_specific_results_spatial$rmse, na.rm = T) %>% round(., 1)),
                                                                                                   mean(c(max(site_specific_results_spatial$rmse, na.rm = T),min(site_specific_results_spatial$rmse, na.rm = T))) %>% round(., 1),
                                                                                                   paste0("> ", max(site_specific_results_spatial$rmse, na.rm = T) %>% round(., 1))))+
  theme_bw(base_size = 10)+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
        legend.position = 'bottom', legend.key.width=unit(2,"cm"))+
  guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, size=18), 
        plot.subtitle = element_text(hjust = 0.5, size=12),
        strip.text.y.right = element_text(angle=360, vjust = 0.5, hjust = 0.5),
        strip.placement = "outside")



jpeg(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_spatial_rmse.jpg"),
    width = 11, height = 16, units = 'in', res = 600)
print(spatial_plot)
dev.off()

#pearson
spatial_plot = ggplot(site_specific_results_spatial %>% filter(!is.na(pearson_r)))+
  facet_grid(nc_id ~ depth, labeller = label_wrap_gen(width=10))+
  geom_sf(data = states, fill = 'transparent')+
  geom_sf(aes(fill = pearson_r),color = 'black', shape = 21, size = 1)+
  scale_fill_gradientn(colours = color_scale(100) , name = "Pearson's r", guide = "colourbar",
                       breaks=c(min(site_specific_results_spatial$pearson_r, na.rm = T),
                                mean(c(max(site_specific_results_spatial$pearson_r, na.rm = T),min(site_specific_results_spatial$pearson_r, na.rm = T))),
                                max(max(site_specific_results_spatial$pearson_r, na.rm = T))), labels=c(paste0("< ", min(site_specific_results_spatial$pearson_r, na.rm = T) %>% round(., 1)),
                                                                                                   mean(c(max(site_specific_results_spatial$pearson_r, na.rm = T),min(site_specific_results_spatial$pearson_r, na.rm = T))) %>% round(., 1),
                                                                                                   paste0("> ", max(site_specific_results_spatial$pearson_r, na.rm = T) %>% round(., 1))))+
  theme_bw(base_size = 10)+
  theme(legend.key = element_blank(), strip.background = element_rect(colour="transparent", fill="transparent"),
        legend.position = 'bottom', legend.key.width=unit(2,"cm"))+
  guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
  theme(plot.title = element_text(hjust = 0.5, size=18), 
        plot.subtitle = element_text(hjust = 0.5, size=12),
        strip.text.y.right = element_text(angle=360, vjust = 0.5, hjust = 0.5),
        strip.placement = "outside")



jpeg(paste0("/home/zhoylman/soil-moisture-validation/figs/drought_model_comaprison/drought_anomoly_model_comparison_6_year_min_summer_clamped_spatial_r.jpg"),
    width = 11, height = 16, units = 'in', res = 600)
print(spatial_plot)
dev.off()

#final site map 
sites_considered = binded %>%
  drop_na(model_drought_anomaly, value) %$%
  site_id %>%
  unique()

stations_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-station-meta-6-years-min-CDF.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(., st_crs('EPSG:4326')) %>%
  st_transform(., st_crs('EPSG:5070')) %>%
  mutate(network = ifelse(network == 'SNTLT', 'SNOTEL', network),
         network = ifelse(network == 'SNTL', 'SNOTEL', network)) %>%
  filter(site_id %in% sites_considered)

states = read_sf('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_20m.json') %>%
  st_transform(., st_crs('EPSG:5070')) %>%
  filter(NAME %notin% c('Puerto Rico', 'Alaska', 'Hawaii'))

site_map = ggplot()+
  geom_sf(data = states, fill = 'transparent')+
  geom_sf(data = stations_meta, aes(fill = network), shape = 21, size = 2)+
  theme_bw(base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom", legend.box = "horizontal",
        plot.subtitle = element_text(hjust = 0.5),
        legend.title=element_blank())+
  scale_fill_manual(values = viridis::turbo(5), name = '')+
  ggtitle(expression(~italic(In)~italic(Situ)~' Soil Moisture Observations'), paste0('n = ', length(stations_meta$network)))+
  guides(fill = guide_legend(override.aes = list(size=5)))

ggsave(site_map, file = '/home/zhoylman/soil-moisture-validation/figs/site_map.jpg', width = 8, height = 8, dpi = 600)

## compute total number of timeseries

standardized_soil_moisture_obs %>% 
  filter(site_id %in% sites_considered) %>%
  select(contains(c('date','site_id','drought'))) %>% pivot_longer(cols = -c(date, site_id)) %>%
  mutate(site_name = paste0(site_id, name)) %>% 
  drop_na() %>% 
  filter(name != 'drought_anomaly_mean') %$%
  site_name %>%
  unique() %>%
  length()
