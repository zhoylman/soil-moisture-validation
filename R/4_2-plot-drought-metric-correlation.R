library(tidyverse)
library(scales)

metric = 'spi'
metric_upper = casefold(metric, upper = TRUE)

metrics = c('spi', 'spei', 'eddi')


getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for(i in 1:length(metrics)){
  metric = metrics[i]
  metric_upper = casefold(metric, upper = TRUE)
  
  out = readRDS(paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlations/',metric,'-cor.RDS'))
  
  full_seasonal_corelation_all = Filter(function(a) any(!is.na(a)), out) %>%
    purrr::map(1) %>%
    lapply(., function(x){return(x %>% mutate(generalized_depth = generalized_depth %>% as.character()))}) %>%
    bind_rows()
  
  density_data = full_seasonal_corelation_all %>%
    filter(standardize_method == 'drought_anomaly') %>%
    group_by(generalized_depth) %>%
    do(density_x = density(.$timescale_mode, bw = 'nrd')$x,
       density_y = density(.$timescale_mode, bw = 'nrd')$y,
       density_bw = density(.$timescale_mode, bw = 'nrd')$bw,
       n = length(.$timescale_mode)) %>%
    unnest(cols = c(density_x, density_y, n, density_bw)) %>%
    filter(density_x > 0 & density_x < 730) 
  
  # test = full_seasonal_corelation_all %>%
  #   filter(standardize_method == 'drought_anomaly') %>%
  #   group_by(generalized_depth) %>%
  #   do(hist_x = hist(.$timescale_mode, bw = 'nrd'))
  
  mode_data = density_data %>%
    group_by(generalized_depth) %>%
    filter(density_y == max(density_y)) %>%
    mutate(simplified_timescale = plyr::round_any(density_x, 10),
           n = n[1])
  
  # test_split = full_seasonal_corelation_all %>%
  #   filter(standardize_method == 'drought_anomaly') %>%
  #   group_split(generalized_depth)
  #   
  # ggplot()+
  #   geom_histogram(data = test_split[[1]], aes(timescale_mode, y = ..density.., fill = generalized_depth), color = 'transparent', fill = 'blue', alpha = 0.2)+
  #   geom_histogram(data = test_split[[2]], aes(timescale_mode, y = ..density.., fill = generalized_depth), color = 'transparent', fill = 'red', alpha = 0.2)+
  #   geom_histogram(data = test_split[[3]], aes(timescale_mode, y = ..density.., fill = generalized_depth), color = 'transparent', fill = 'green', alpha = 0.2)+
  #   geom_density(data = test_split[[1]], aes(timescale_mode, color = generalized_depth), color = 'blue',)+
  #   geom_density(data = test_split[[2]], aes(timescale_mode, color = generalized_depth), color = 'red')+
  #   geom_density(data = test_split[[3]], aes(timescale_mode, color = generalized_depth), color = 'green')+
  #   geom_point(data = mode_data, aes(x = density_x, y = density_y, color = generalized_depth))+
  #   geom_text(data = mode_data, aes(x = density_x+70, y = density_y, label = paste0(simplified_timescale, ' days')))
  
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
  
  
  monthly_table_data = full_monthly_corelation_all %>%
    #extract relevant standized soil moisture metric
    filter(standardize_method == 'drought_anomaly')%>%
    #filter for at least 3 years of data
    filter(n > 93) %>%
    #compute site and probe specific optimal timescales
    group_by(month, name, site_id) %>%
    {if(metric == 'eddi') filter(., pearson_r == min(pearson_r)) else filter(., pearson_r == max(pearson_r))} %>%
    ungroup() %>%
    #generalize across sites
    group_by(month, generalized_depth) %>%
    # summarise(median_r = median(pearson_r),
    #           #median timescale?
    #           median_timescale = median(timescale) %>% plyr::round_any(., 10))
    do(density_x = density(.$timescale, bw = 'nrd')$x,
       density_y = density(.$timescale, bw = 'nrd')$y) %>%
    unnest(cols = c(density_x, density_y)) %>%
    ungroup() %>%
    group_by(month, generalized_depth) %>%
    filter(density_y == max(density_y)) %>%
    mutate(simplified_timescale = plyr::round_any(density_x, 10),
           str_match = paste0('month=', month, ',generalized_depth=',generalized_depth,',timescale=',simplified_timescale)) %>%
    select(month, generalized_depth, str_match)
  
  monthly_table_final = full_monthly_corelation_all %>%
    #extract relevant standized soil moisture metric
    filter(standardize_method == 'drought_anomaly') %>%
    #compute str_match
    mutate(str_full = paste0('month=', month, ',generalized_depth=',generalized_depth,',timescale=',timescale)) %>%
    #bind kernal density maxima to data
    left_join(., monthly_table_data, by = c('month', 'generalized_depth')) %>%
    #filter for month, generalized depth and simplified_timescale
    filter(str_full == str_match) %>%
    group_by(month, generalized_depth) %>%
    summarise(median_r = median(pearson_r, na.rm = T),
              median_timescale = median(timescale),
              max_timescale = max(timescale),
              min_timescale = min(timescale)) %>%
    mutate(dummy_date = paste0('2022-',month,'-01') %>% lubridate::ymd(.))
  
  # plot it
  density_plot = ggplot(data = density_data)+
    # geom_histogram(data = test_split[[1]], aes(timescale_mode, y = ..density.., fill = generalized_depth), 
    #                color = 'transparent', fill = 'blue', alpha = 0.1, show.legend = F)+
    # geom_histogram(data = test_split[[2]], aes(timescale_mode, y = ..density.., fill = generalized_depth), 
    #                color = 'transparent', fill = 'forestgreen', alpha = 0.1, show.legend = F)+
    # geom_histogram(data = test_split[[3]], aes(timescale_mode, y = ..density.., fill = generalized_depth), 
    #                color = 'transparent', fill = '#800080', alpha = 0.1, show.legend = F)+
    geom_line(aes(x = density_x, y = density_y, color = generalized_depth))+
    geom_point(data = mode_data, aes(x = density_x, y = density_y, color = generalized_depth))+
    geom_text(data = mode_data, aes(x = density_x+70, y = density_y, label = paste0(simplified_timescale, ' days')))+
    ggtitle(paste0(metric_upper, ' Optimal Timescales (May - Oct)'))+
    theme_bw(base_size = 14)+
    theme(legend.position = c(0.72,0.8),
          plot.title = element_text(size = 14, hjust = 0.5))+
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
    guides(color=guide_legend(title="Soil Depth", reverse = TRUE))+
    labs(x = 'Timescale (Days)', y = 'Density')+
    scale_colour_manual(values = c('blue', 'forestgreen', '#800080'))
  
  seasonal_plot = ggplot(monthly_final)+
    geom_ribbon(aes(x = dummy_date, y = median_r, ymax = upper_r, ymin = lower_r, fill = generalized_depth_timescale), show.legend = F, alpha = 0.5)+
    geom_line(aes(x = dummy_date, y = median_r), show.legend = F) +
    facet_wrap(~fct_rev(generalized_depth_timescale), nrow = 3)+
    theme_bw(base_size = 14)+
    theme(strip.background = element_rect(colour="white", fill="white"))+
    labs(x = 'Month', y = 'Correlation (r)')+ 
    scale_x_date(date_labels = "%b")+
    scale_colour_manual(values = c('blue', 'forestgreen', '#800080'))+
    scale_fill_manual(values = c('blue', 'forestgreen', '#800080'))
  
  seasonal_table = ggplot(data = monthly_table_final, aes(x = dummy_date, y = generalized_depth, fill = median_r, label = median_timescale))+
    geom_tile(width = 31)+
    geom_text()+
    labs(x = 'Month', y = NULL)+
    {if(metric == 'eddi'){
      scale_fill_gradientn(colours = rev(c('red','orange', 'white', 'cyan', 'blue')),
                           breaks = seq(-0.7,-0.1, by = 0.1),
                           name = 'Correlation (r)',
                           limits = c(-0.7,-0.1))
    } else {
      scale_fill_gradientn(colours = (c('red','orange', 'white', 'cyan', 'blue')),
                           breaks = seq(0.7, 0.1, by = -0.1),
                           name = 'Correlation (r)',
                           limits = c(0.1,0.7))
    }}+
    theme_bw(base_size = 16)+
    theme(legend.position = 'bottom',
          legend.key.width= unit(2, 'cm'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(angle=60, vjust = 0.5, hjust = 0.5))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))+
    ggtitle(paste0(metric_upper,' Monthly Correlation by Depth'))+
    theme(plot.title = element_text(size = 14, hjust = 0.5))+
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5))+
    scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = expansion(add = 4))
  
  top_row = cowplot::plot_grid(density_plot, seasonal_plot, rel_widths = c(0.6,0.4), labels = c('(a)', '(b)'))
  full_plot = cowplot::plot_grid(top_row, seasonal_table, nrow = 2, labels = c('', '(c)'))
  
  ggsave(full_plot, file = paste0('/home/zhoylman/soil-moisture-validation/figs/',metric,'-cor.png'), width = 9, height = 10)
}
