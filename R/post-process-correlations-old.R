library(data.table)
library(sf)
library(magrittr)
library(tidyverse)

## Full time period

spi_cor = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/spei-correlation-matrix.csv')

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

extract_density_n = function(x){
  return(x[[1]]$n)
}

max_vals = spi_cor %>%
  pivot_longer(cols = -c(site_id, timescale)) %>%
  filter(name %like% 'drought_anomaly') %>%
  group_by(site_id, name) %>%
  slice(which.max(value)) %>%
  mutate(depth = parse_number(name)) %>%
  mutate(generalized_depth = ifelse(depth <= 4, 'Shallow',
                                    ifelse(depth > 4 & depth <= 20, 'Middle',
                                           ifelse(depth > 20, 'Deep', NA)))) %>%
  ungroup() %>%
  mutate(timescale_num = parse_number(timescale))
  
optimal_kde = max_vals %>%
  group_by(generalized_depth) %>%
  do(denstiy = density(.$timescale_num)) 

optimal_cors = max_vals %>%
  group_by(generalized_depth) %>%
  summarise(cor = median(value) %>%
              round(., 2))

mode_timescale = optimal_kde %>%
  group_by(generalized_depth) %>%
  do(mode_timescale = extract_mode(.$denstiy)) %>%
  mutate(timescale = mode_timescale  %$%
           x,
         density = mode_timescale  %$%
           y,
         timescale_round = plyr::round_any(timescale, 10))

plot = ggplot()+
  geom_line(data = extract_density(optimal_kde$denstiy[3]), aes(x = x, y = y, color = 'Shallow'))+
  geom_line(data = extract_density(optimal_kde$denstiy[2]), aes(x = x, y = y, color = 'Middle'))+
  geom_line(data = extract_density(optimal_kde$denstiy[1]), aes(x = x, y = y, color = 'Deep'))+
  #add labels for optimal
  geom_point(data = extract_mode(optimal_kde$denstiy[3]), aes(x = x, y = y, color = 'Shallow'))+
  geom_text(data = extract_mode(optimal_kde$denstiy[3]), aes(x = x+50, y = y, label = paste0(plyr::round_any(x, 10), ' Days'), color = 'Shallow'))+
  geom_point(data = extract_mode(optimal_kde$denstiy[2]), aes(x = x, y = y, color = 'Middle'))+
  geom_text(data = extract_mode(optimal_kde$denstiy[2]), aes(x = x+50, y = y, label = paste0(plyr::round_any(x, 10), ' Days'), color = 'Middle'))+
  geom_point(data = extract_mode(optimal_kde$denstiy[1]), aes(x = x, y = y, color = 'Deep'))+
  geom_text(data = extract_mode(optimal_kde$denstiy[1]), aes(x = x+50, y = y, label = paste0(plyr::round_any(x, 10), ' Days'), color = 'Deep'))+
  xlim(x = c(0,730))+
  labs(x = 'Timescale (Days)', y = 'Density')+
  theme_bw(base_size = 16)+
  scale_color_manual(breaks = c("Shallow", "Middle", "Deep"),
                     values = c("forestgreen", "blue", "purple"),
                     labels = paste0(c("Shallow: r = ", "Middle: r = ", "Deep: r = "),
                                     c(optimal_cors$cor %>% rev())),
                     name = "Probe Depth")+
  theme(legend.position = c(0.8,0.8),
        legend.box.background = element_rect(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('SPEI Optimal Correlations')

ggsave(plot, file = '/home/zhoylman/soil-moisture-validation/figs/spei_density.png')  

## seasonal correlations all time with optimal
seasonal_cor = readRDS('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/spei-correlation-matrix-seasonal.RDS')

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

seasonal_generalized_optimal = seasonal_cor %>%
  group_by(site_id, month) %>%
  do(important_cor = extract_cors_seasonal(.$cor_filtered))

final_generalized_seasonal = seasonal_generalized_optimal %>%
  group_by(month) %>%
  do(summary = bind_rows(.$important_cor)) %>%
  unnest(summary) %>%
  group_by(month, generalized_depth)%>%
  summarise(median = median(value, na.rm = T),
            upper = quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T)) %>%
  drop_na(generalized_depth) %>%
  mutate(generalized_depth = factor(generalized_depth, levels = c('Shallow', 'Middle', 'Deep')))

plot_seasonal = ggplot(data = final_generalized_seasonal, 
       aes(x = month, y = median, ymax = upper, 
           ymin = lower, color = generalized_depth,
           fill = generalized_depth))+
  geom_ribbon(alpha = 0.6) +
  geom_line(color = 'black') + 
  facet_wrap(~generalized_depth)+
  theme_bw(base_size = 16)+
  labs(x = 'Month', y = 'Correlation (r)')+
  scale_fill_manual(breaks = c("Shallow", "Middle", "Deep"),
                     values = c("forestgreen", "blue", "purple"))+
  scale_color_manual(breaks = c("Shallow", "Middle", "Deep"),
                    values = c("forestgreen", "blue", "purple"))+
  scale_x_continuous(breaks = seq(1,11,by =2))+
  theme(legend.position='none',
        strip.background = element_blank(),
        axis.ticks.x=element_blank())

ggsave(plot_seasonal, file = '/home/zhoylman/soil-moisture-validation/figs/spei_generalized_seasonal.png', height = 4, width = 8, units = 'in')  


## Seasonal best time scales
min_stations = 0

seasonal_cor = readRDS('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/spei-correlation-matrix-seasonal.RDS')
station_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/station_meta.csv') %>%
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(st_crs(4326))

montana = st_read('https://eric.clst.org/assets/wiki/uploads/Stuff/gz_2010_us_040_00_500k.json') %>%
  filter(NAME %in% c('Montana', 'North Dakota', 'South Dakota', 'Wyoming', 'Idaho'))
  #filter(NAME %in% c('Montana'))


mt_sites = st_intersection(station_meta, montana)

seasonal_cor = seasonal_cor %>% 
  filter(site_id %in% mt_sites$site_id)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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

site_specific_seasonal = seasonal_expand %>%
  group_by(site_id, generalized_depth, month) %>%
  slice(which.max(value)) %>%
  mutate(timescale_num = parse_number(timescale))

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

table = ggplot(data = final, aes(x = month, y = paste0('SPEI\n', generalized_depth, ' VWC'), fill = median_cor, label = optimal_timescale))+
  geom_tile()+
  geom_text()+
  scale_x_continuous(breaks = 1:12, limits = c(1,12))+
  labs(x = 'Month', y = NULL)+
  # scale_fill_gradient2(low = ("red"), 
  #                      mid = "white", 
  #                      high = ("blue"), 
  #                      name = NULL,#"Pearson's r",
  #                      midpoint = 0.5)+
  scale_fill_gradientn(colours = c('red','orange', 'white', 'cyan', 'blue'),
                       breaks = seq(0.20,0.65, by = 0.15),
                       name = NULL,
                       limits = c(0.2,0.65))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'bottom',
        legend.key.width= unit(2, 'cm'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(angle=0, vjust = 0.5, hjust = 0.5))+
  ggtitle('UMRB Results', subtitle = paste('Montana', 'North Dakota', 'South Dakota', 'Wyoming', 'Idaho', sep = ', '))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

table

ggsave(table, file = '/home/zhoylman/soil-moisture-validation/figs/spei_seasonal_table_umrb.png', height = 5, width = 7, units = 'in')  

seasonal_optimal_scatter = ggplot(data = final, aes(x = month, y = optimal_timescale, color = paste0('SPI\n', generalized_depth, ' VWC'), size = median_cor))+
  geom_point()+
  geom_line(size = 0.5, alpha = 0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = 'bottom',
        legend.key.width= unit(2, 'cm'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(angle=0, vjust = 0.5, hjust = 0.5),
        legend.box="vertical", legend.margin=margin())+
  labs(x = 'Month', y = 'Optimal Timescale', color = NULL, size = 'Median Correlation (r)')+
  scale_x_continuous(breaks = 1:12)
  
seasonal_optimal_scatter

ggsave(seasonal_optimal_scatter, file = '/home/zhoylman/soil-moisture-validation/figs/spi_seasonal_scatter.png', height = 5, width = 7, units = 'in')  
