#load libs
library(tidyverse)
library(scales)

#define drought metric names
metrics = c('spi', 'spei', 'eddi')

#import station meta
stations_meta = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/beta-standardized-station-meta-6-years-min-CDF-w-mean.csv')

#define USCRN if you want to only use USCRN
uscrn = stations_meta %>%
  filter(network == 'USCRN')

#define mode function 
getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#for loop per drought metric
for(i in 1:length(metrics)){
  #define metric lower and upper case
  metric = metrics[i]
  metric_upper = casefold(metric, upper = TRUE)
  
  #import data from script 4_1
  out = readRDS(paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlations/',metric,'-6-years-min-cor-rmse-clamped-w-mean-beta-standardized.RDS'))
  
  #summarize the full season results 
  full_seasonal_corelation_all = Filter(function(a) any(!is.na(a)), out) %>%
    purrr::map(1) %>%
    lapply(., function(x){return(x %>% mutate(generalized_depth = generalized_depth %>% as.character(),
                                              standardize_method = standardize_method %>% as.character()))}) %>%
    bind_rows()
  
  #compute density data to compute optimal timescale for whole season
  density_data = full_seasonal_corelation_all %>%
    filter(standardize_method == 'drought_anomaly') %>%
    #filter(site_id %in% uscrn$site_id) %>%
    group_by(generalized_depth) %>%
    do(density_x = density(.$timescale_mode, bw = 'nrd')$x,
       density_y = density(.$timescale_mode, bw = 'nrd')$y,
       density_bw = density(.$timescale_mode, bw = 'nrd')$bw,
       n = length(.$timescale_mode)) %>%
    unnest(cols = c(density_x, density_y, n, density_bw)) %>%
    filter(density_x > 0 & density_x < 730) 
  
  #compute the mode of the KDE (this is the optimal seasonal timescale by depth)
  mode_data = density_data %>%
    group_by(generalized_depth) %>%
    filter(density_y == max(density_y)) %>%
    mutate(simplified_timescale = plyr::round_any(density_x, 10),
           n = n[1])
  
  # compute the monthly correlations (we will use this to compute average rmse)
  full_monthly_corelation_all = Filter(function(a) any(!is.na(a)), out) %>%
    purrr::map(2) %>%
    lapply(., function(x){return(x %>% mutate(generalized_depth = generalized_depth %>% as.character(),
                                              standardize_method = standardize_method %>% as.character()))}) %>%
    bind_rows()
  
  #compute full season rmse given the optimal timescale from mode data (joined with flexible below)
  rigid_season_rmse = full_monthly_corelation_all %>%
    filter(standardize_method == 'drought_anomaly',
           timescale %in% mode_data$simplified_timescale,
           month >= 5 & month <= 10) %>% 
    #filter for modal timescale by depth to compute seasonal average rmse for optimal timescale
    filter(generalized_depth == mode_data$generalized_depth[1] & 
           timescale == mode_data$simplified_timescale [1] | 
           generalized_depth == mode_data$generalized_depth[2] & 
           timescale == mode_data$simplified_timescale [2]| 
           generalized_depth == mode_data$generalized_depth[3] & 
           timescale == mode_data$simplified_timescale [3]|
           generalized_depth == mode_data$generalized_depth[4] & 
           timescale == mode_data$simplified_timescale [4]
           ) %>%
    #filter for a minimum of 30 obs to compute rmse (less than 1% of the data)
    filter(n > 30) %>%
    group_by(site_id, generalized_depth) %>%
    summarize(#median_seasonal_rigid_rmse = median(linearFit_RMSE, na.rm = T),
              median_seasonal_rigid_r = median(pearson_r, na.rm = T))
  
  #re compute monthly filtered data (same as above but for all months (not just May - Oct))
  monthly_filtered = full_monthly_corelation_all %>%
    filter(standardize_method == 'drought_anomaly' ,
           timescale %in% mode_data$simplified_timescale) %>%
    #filter(site_id %in% uscrn$site_id) %>%
    filter(generalized_depth == mode_data$generalized_depth[1] & timescale == mode_data$simplified_timescale[1] |
             generalized_depth == mode_data$generalized_depth[2] & timescale == mode_data$simplified_timescale[2] |
             generalized_depth == mode_data$generalized_depth[3] & timescale == mode_data$simplified_timescale[3] |
             generalized_depth == mode_data$generalized_depth[4] & timescale == mode_data$simplified_timescale[4])
  
  #compute montly final correlations (this is for fig 2-4(b))
  monthly_final = monthly_filtered %>%
    group_by(generalized_depth, timescale, month) %>%
    summarise(median_r = quantile(pearson_r, 0.5, na.rm = T),
              upper_r = quantile(pearson_r, 0.75, na.rm = T),
              lower_r = quantile(pearson_r, 0.25, na.rm = T)) %>%
    mutate(dummy_date = paste0('2022-',month,'-01') %>% lubridate::ymd(.),
           generalized_depth_timescale = paste0(generalized_depth, ' [', timescale, ' day ', metric_upper, ']'))
  
  #compute multi part string match for month, generalized depth, and timescale
  monthly_table_data = full_monthly_corelation_all %>%
    #extract relevant standized soil moisture metric
    filter(standardize_method == 'drought_anomaly')%>%
    #filter(site_id %in% uscrn$site_id) %>%
    #compute site and probe specific optimal timescales
    group_by(month, name, site_id) %>%
    {if(metric == 'eddi') filter(., pearson_r == min(pearson_r)) else filter(., pearson_r == max(pearson_r))} %>%
    ungroup() %>%
    #generalize across sites
    group_by(month, generalized_depth) %>%
    do(density_x = density(.$timescale, bw = 'nrd')$x,
       density_y = density(.$timescale, bw = 'nrd')$y) %>%
    unnest(cols = c(density_x, density_y)) %>%
    ungroup() %>%
    group_by(month, generalized_depth) %>%
    filter(density_y == max(density_y)) %>%
    mutate(simplified_timescale = plyr::round_any(density_x, 10),
           str_match = paste0('month=', month, ',generalized_depth=',generalized_depth,',timescale=',simplified_timescale)) %>%
    dplyr::select(month, generalized_depth, str_match)
  
  # use the above calculated string match to extract relevant information for plot
  monthly_table_final = full_monthly_corelation_all %>%
    #extract relevant standized soil moisture metric
    filter(standardize_method == 'drought_anomaly') %>%
    #filter(site_id %in% uscrn$site_id) %>%
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
  
  # use the above calculated string match to extract relevant information for rmse comparison
  monthly_flexible_rmse = full_monthly_corelation_all %>%
    #extract relevant standized soil moisture metric
    filter(standardize_method == 'drought_anomaly') %>%
    #filter(site_id %in% uscrn$site_id) %>%
    #compute str_match
    mutate(str_full = paste0('month=', month, ',generalized_depth=',generalized_depth,',timescale=',timescale)) %>%
    #bind kernal density maxima to data
    left_join(., monthly_table_data, by = c('month', 'generalized_depth')) %>%
    #filter for month, generalized depth and simplified_timescale
    filter(str_full == str_match) %>%
    #filter by months of interst
    filter(month >= 5 & month <= 10) %>%
    group_by(site_id, generalized_depth) %>%
    summarise(#median_seasonal_flexible_rmse = median(linearFit_RMSE, na.rm = T),
              median_seasonal_flexible_r = median(pearson_r, na.rm = T))
  
  #left join the two rmse comparison datasets
  #final_rmse_comparison = left_join(rigid_season_rmse, monthly_flexible_rmse, by = c('site_id', 'generalized_depth'))
  
  #write final rmse comparison csv
  #write_csv(final_rmse_comparison, paste0('~/soil-moisture-validation-data/processed/rmse-comparison/',metric,'-6-years-min-season-rmse-clamped.csv'))
  
  #rename and reorder for consistancy
  density_data = density_data %>%
    mutate(generalized_depth = ifelse(generalized_depth == "Middle (8in - 20in)",
                                      "Middle (8-20in)", generalized_depth),
           generalized_depth = fct_relevel(generalized_depth, c("Depth Averaged",
                                                                "Shallow (0-4in)",
                                                                "Middle (8-20in)",
                                                                "Deep (>20in)")))
  
  
  #reorder 
  monthly_final = monthly_final %>%
    mutate(generalized_depth = ifelse(generalized_depth == "Middle (8in - 20in)",
                                      "Middle (8-20in)", generalized_depth),
           generalized_depth = as.factor(generalized_depth),
           generalized_depth = fct_relevel(generalized_depth, c("Depth Averaged",
                                                                "Shallow (0-4in)",
                                                                "Middle (8-20in)",
                                                                "Deep (>20in)")),
           generalized_depth_timescale = paste0(generalized_depth, ' [', timescale, ' day ', metric_upper, ']') %>% as.factor())
  
  monthly_final$generalized_depth_timescale = factor(monthly_final$generalized_depth_timescale, 
                                                      levels = c(monthly_final[order(monthly_final$generalized_depth), "generalized_depth_timescale"] %>% unique())$generalized_depth_timescale)
  
  
  monthly_table_final = monthly_table_final %>%
    mutate(generalized_depth = ifelse(generalized_depth == "Middle (8in - 20in)",
                                      "Middle (8-20in)", generalized_depth),
           #max is 0.703 - for color truncate at 0.7
           median_r = ifelse(median_r > 0.7, 0.7, median_r),
           generalized_depth = fct_relevel(generalized_depth, c("Depth Averaged",
                                                                "Shallow (0-4in)",
                                                                "Middle (8-20in)",
                                                                "Deep (>20in)")))
  
  mode_data = mode_data %>%
    mutate(generalized_depth = ifelse(generalized_depth %>% as.character() == "Middle (8in - 20in)",
                                      "Middle (8-20in)", generalized_depth %>% as.character()),
           generalized_depth = as.factor(generalized_depth),
           generalized_depth = fct_relevel(generalized_depth, c("Depth Averaged",
                                                                "Shallow (0-4in)",
                                                                "Middle (8-20in)",
                                                                "Deep (>20in)")))

  
  # plot it
  density_plot = ggplot(data = density_data)+
    geom_line(aes(x = density_x, y = density_y, color = generalized_depth))+
    geom_point(data = mode_data, aes(x = density_x, y = density_y, color = generalized_depth))+
    geom_text(data = mode_data, aes(x = density_x+70, y = density_y, label = paste0(simplified_timescale, ' days')))+
    ggtitle(paste0(metric_upper, ' Optimal Timescales (May - Oct)'))+
    theme_bw(base_size = 14)+
    theme(legend.position = c(0.72,0.8),
          legend.margin=margin(r = 0.6, t = 0.3, l = 0.3, b = 0.3, unit='cm'),
          plot.title = element_text(size = 14, hjust = 0.5))+
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
    guides(color=guide_legend(title="Soil Depth", reverse = F))+
    labs(x = 'Timescale (Days)', y = 'Density')+
    scale_colour_manual(values = c('black', 'blue', 'forestgreen', '#800080'))
  
  seasonal_plot = ggplot(monthly_final)+
    geom_ribbon(aes(x = dummy_date, y = median_r, ymax = upper_r, ymin = lower_r, fill = generalized_depth_timescale), show.legend = F, alpha = 0.5)+
    geom_line(aes(x = dummy_date, y = median_r), show.legend = F) +
    facet_wrap(~(generalized_depth_timescale), nrow = 4)+
    theme_bw(base_size = 14)+
    theme(strip.background = element_rect(colour="white", fill="white"))+
    labs(x = 'Month', y = 'Correlation (r)')+ 
    scale_x_date(date_labels = "%b")+
    #scale_colour_manual(values = c('black', 'blue', 'forestgreen', '#800080'))+
    scale_fill_manual(values =  c('black', 'blue', 'forestgreen', '#800080'))
  
  seasonal_table = ggplot(data = monthly_table_final, aes(x = dummy_date, y = generalized_depth, fill = median_r, label = median_timescale))+
    geom_tile(width = 31)+
    geom_text()+
    scale_y_discrete(limits=rev)+
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
  full_plot = cowplot::plot_grid(top_row, seasonal_table, nrow = 2, rel_widths = c(0.6,0.4), labels = c('', '(c)'))
  
  ggsave(full_plot, file = paste0('/home/zhoylman/soil-moisture-validation/figs-revision1/',metric,'-6-years-min-cor-clamped-w-mean.png'), width = 9, height = 11)
}
