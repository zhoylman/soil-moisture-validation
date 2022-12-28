library(tidyverse)
library(scales)

metric = 'spi'
metric_upper = 'SPI'

readRDS(paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlations/',metric,'-cor.RDS'))

full_seasonal_corelation_all = Filter(function(a) any(!is.na(a)), out) %>%
  purrr::map(1) %>%
  lapply(., function(x){return(x %>% mutate(generalized_depth = generalized_depth %>% as.character()))}) %>%
  bind_rows()

density_data = full_seasonal_corelation_all %>%
  filter(standardize_method == 'drought_anomaly') %>%
  group_by(generalized_depth) %>%
  do(density_x = density(.$timescale_mode)$x,
     density_y = density(.$timescale_mode)$y) %>%
  unnest(cols = c(density_x, density_y)) %>%
  filter(density_x > 0 & density_x < 730) 

mode_data = density_data %>%
  group_by(generalized_depth) %>%
  filter(density_y == max(density_y)) %>%
  mutate(simplified_timescale = plyr::round_any(density_x, 10))

density_plot = ggplot(data = density_data)+
  geom_line(aes(x = density_x, y = density_y, color = generalized_depth))+
  geom_point(data = mode_data, aes(x = density_x, y = density_y, color = generalized_depth))+
  geom_text(data = mode_data, aes(x = density_x+70, y = density_y, label = paste0(simplified_timescale, ' days')))+
  ggtitle(paste0(metric_upper, ' Optimal Timescales (May - Oct)'))+
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.72,0.8))+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  guides(color=guide_legend(title="Soil Depth", reverse = TRUE))+
  labs(x = 'Timescale (Days)', y = 'Density')+
  scale_colour_manual(values = c('blue', 'forestgreen', '#800080'))

full_monthly_corelation_all = Filter(function(a) any(!is.na(a)), out) %>%
  purrr::map(2) %>%
  bind_rows()

monthly_filtered = full_monthly_corelation_all %>%
  filter(standardize_method == 'drought_anomaly' ,
         timescale %in% mode_data$simplified_timescale) %>%
  filter(generalized_depth == mode_data$generalized_depth[1] & timescale == mode_data$simplified_timescale[1] |
           generalized_depth == mode_data$generalized_depth[2] & timescale == mode_data$simplified_timescale[2] |
           generalized_depth == mode_data$generalized_depth[3] & timescale == mode_data$simplified_timescale[3])

monthly_final = monthly_filtered %>%
  filter(n > 93) %>% 
  group_by(generalized_depth, timescale, month) %>%
  summarise(median_r = quantile(pearson_r, 0.5, na.rm = T),
            upper_r = quantile(pearson_r, 0.75, na.rm = T),
            lower_r = quantile(pearson_r, 0.25, na.rm = T)) %>%
  mutate(dummy_date = paste0('2022-',month,'-01') %>% lubridate::ymd(.),
         generalized_depth_timescale = paste0(generalized_depth, ' [', timescale, ' day ', metric_upper, ']'))

seasonal_plot = ggplot(monthly_final)+
  geom_ribbon(aes(x = dummy_date, y = median_r, ymax = upper_r, ymin = lower_r, fill = generalized_depth_timescale), show.legend = F, alpha = 0.5)+
  geom_line(aes(x = dummy_date, y = median_r), show.legend = F) +
  facet_wrap(~fct_rev(generalized_depth_timescale), nrow = 3)+
  theme_bw(base_size = 16)+
  theme(strip.background = element_rect(colour="white", fill="white"))+
  labs(x = NULL, y = 'Correlation (r)')+ 
  scale_x_date(date_labels = "%b")+
  scale_colour_manual(values = c('blue', 'forestgreen', '#800080'))+
  scale_fill_manual(values = c('blue', 'forestgreen', '#800080'))


full_plot = cowplot::plot_grid(density_plot, seasonal_plot, rel_widths = c(0.6,0.4))

ggsave(full_plot, file = paste0('/home/zhoylman/soil-moisture-validation/figs/',metric,'-cor.png'), width = 9, height = 5)
