library(data.table)
library(sf)
library(magrittr)
library(tidyverse)

# define functions

extract_density = function(x){
  x = x[[1]]
  out = data.frame(x = x$x, y = x$y) %>%
    as_tibble()
  return(out)
}

extract_mode = function(x){
  x = x[[1]]
  out = data.frame(x = x$x, y = x$y) %>%
    as_tibble()
  best = out[which(out$y == max(out$y)),]
  return(best)
}

getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

extract_density_n = function(x){
  return(x[[1]]$n)
}

extract_cors_seasonal = function(x){
  tryCatch({
    x_ = x[[1]][[1]]
    out = x_ %>%
      filter(timescale %in% paste0('t_', mode_timescale$timescale_round)) %>%
      select(contains(c('site_id' ,'timescale', 'drought'))) %>%
      pivot_longer(cols = -c(site_id, timescale)) %>%
      mutate(depth = parse_number(name)) %>%
      mutate(generalized_depth = ifelse(depth <= 4, 'Shallow',
                                        ifelse(depth > 4 & depth <= 20, 'Middle',
                                               ifelse(depth > 20, 'Deep', NA)))) %>%
      ungroup() %>%
      filter(paste0(timescale, generalized_depth) %in% paste0('t_', mode_timescale$timescale_round, mode_timescale$generalized_depth))
    return(out)
  }, error = function(e){
    out = tibble(site_id = NA,
                 timescale = NA,
                 name = NA,
                 value = NA, 
                 depth = NA,
                 generalized_depth = NA)
    return(out)
  })
}

# define variables of interest
vars = c('spi', 'spei', 'eddi')
for(i in 1:3){
  var = vars[i]
  min_stations = 0
  
  ## Seasonal best time scales
  seasonal_cor_all = readRDS(paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/',var,'-correlation-matrix-seasonal.RDS'))
  station_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/station_meta.csv') %>%
    st_as_sf(., coords = c('longitude', 'latitude')) %>%
    st_set_crs(st_crs(4326))
  
  #filter to domain of interest
  montana = st_read('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_500k.json') %>%
    filter(NAME %in% c('Montana', 'North Dakota', 'South Dakota', 'Wyoming', 'Idaho'))
  mt_sites = st_intersection(station_meta, montana)
  seasonal_cor = seasonal_cor_all %>% 
    filter(site_id %in% mt_sites$site_id)
  
  # post processes corelations using KDE to determine best timescales and compute
  # average corealtions across season and depth
  post_process_correlations = function(seasonal_cor, var){
    seasonal_expand = seasonal_cor %>%
      dplyr::select(month, cor_filtered) %>%
      unnest(cor_filtered)%>%
      unnest(cor_filtered) %>%
      dplyr::select(contains(c('month', 'site_id', 'timescale', 'drought'))) %>%
      pivot_longer(cols = -c(month, site_id, timescale)) %>%
      mutate(depth = parse_number(name)) %>%
      mutate(generalized_depth = ifelse(depth <= 4, 'Shallow',
                                        ifelse(depth > 4 & depth <= 20, 'Middle',
                                               ifelse(depth > 20, 'Deep', NA)))) 
    
    if(var == 'eddi'){
      site_specific_seasonal = seasonal_expand %>%
        group_by(site_id, generalized_depth, month) %>%
        slice(which.min(value)) %>%
        mutate(timescale_num = parse_number(timescale))
    } else {
      site_specific_seasonal = seasonal_expand %>%
        group_by(site_id, generalized_depth, month) %>%
        slice(which.max(value)) %>%
        mutate(timescale_num = parse_number(timescale))
    }
    
    kde_seasonal_optimal = site_specific_seasonal %>%
      group_by(generalized_depth, month) %>%
      do(denstiy = density(.$timescale_num)) %>%
      ungroup() %>%
      group_by(generalized_depth, month) %>%
      mutate(n = extract_density_n(denstiy)) %>%
      filter(n > min_stations) %>%
      do(mode_timescale = extract_mode(.$denstiy)) %>%
      unnest(mode_timescale) %>%
      mutate(optimal_timescale = plyr::round_any(x, 10))
    
    seasonal_general_cor = seasonal_expand %>%
      mutate(depth_timescale = paste0(generalized_depth, ':', timescale)) %>%
      filter(depth_timescale %in% paste0(kde_seasonal_optimal$generalized_depth, ':t_', kde_seasonal_optimal$optimal_timescale)) %>%
      drop_na(value) %>%
      group_by(month, generalized_depth) %>%
      summarise(median_cor = median(value))
    
    
    final = left_join(kde_seasonal_optimal, seasonal_general_cor, by = c('month', 'generalized_depth'))
    return(final)
  }
  
  #compute optinal seasonal correlations
  final = post_process_correlations(seasonal_cor, var)
  
  # make "table" plot
  table_plot = function(final, var){
    ggplot(data = final, aes(x = month, y = paste0(casefold(var, upper = TRUE),'\n', generalized_depth, ' VWC'), fill = median_cor, label = optimal_timescale))+
      geom_tile()+
      geom_text()+
      scale_x_continuous(breaks = 1:12)+
      labs(x = 'Month', y = NULL)+
      {if(var == 'eddi'){
        scale_fill_gradientn(colours = rev(c('red','orange', 'white', 'cyan', 'blue')),
                             breaks = seq(-0.65,-0.15, by = 0.1),
                             name = NULL,
                             limits = c(-0.65,-0.15))
      } else {
        scale_fill_gradientn(colours = (c('red','orange', 'white', 'cyan', 'blue')),
                             breaks = seq(0.15,0.65, by = 0.1),
                             name = NULL,
                             limits = c(0.15,0.65))
      }}+
      theme_bw(base_size = 16)+
      theme(legend.position = 'bottom',
            legend.key.width= unit(2, 'cm'))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text.y = element_text(angle=0, vjust = 0.5, hjust = 0.5))+
      ggtitle('UMRB Results', subtitle = paste('Montana', 'North Dakota', 'South Dakota', 'Wyoming', 'Idaho', sep = ', '))+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  table = table_plot(final, var)
  
  #export
  ggsave(table, file = paste0('/home/zhoylman/soil-moisture-validation/figs/',var,'_seasonal_table_umrb_final.png'), height = 5, width = 7, units = 'in')  
}


# 
# seasonal_optimal_scatter = ggplot(data = final, aes(x = month, y = optimal_timescale, color = paste0('SPI\n', generalized_depth, ' VWC'), size = median_cor))+
#   geom_point()+
#   geom_line(size = 0.5, alpha = 0.5)+
#   theme_bw(base_size = 16)+
#   theme(legend.position = 'bottom',
#         legend.key.width= unit(2, 'cm'))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.text.y = element_text(angle=0, vjust = 0.5, hjust = 0.5),
#         legend.box="vertical", legend.margin=margin())+
#   labs(x = 'Month', y = 'Optimal Timescale', color = NULL, size = 'Median Correlation (r)')+
#   scale_x_continuous(breaks = 1:12)
# 
# seasonal_optimal_scatter
# 
# ggsave(seasonal_optimal_scatter, file = '/home/zhoylman/soil-moisture-validation/figs/spi_seasonal_scatter.png', height = 5, width = 7, units = 'in')  
