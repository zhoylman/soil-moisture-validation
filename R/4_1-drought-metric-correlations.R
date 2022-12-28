library(tidyverse)
library(scales)

getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

metric = 'eddi'
metric_upper = 'EDDI'

drought = read_csv(paste0('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/',metric,'-data-long-10s.csv'))%>%
  rename(drought = metric)

soil_moisture = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide.csv') %>%
  pivot_longer(., cols = -c(site_id,date)) %>%
  mutate(time = date) %>%
  select(site_id,time,name,value)

ids = unique(soil_moisture$site_id)

# i = 100
#
# tictoc::tic()
# 
# temp_data = drought %>%
#   filter(site_id == ids[i]) %>%
#   left_join(.,  soil_moisture %>%
#               filter(site_id == ids[i]),
#             by = c('site_id', 'time'))
# 
# full_seasonal_corelation = temp_data %>%
#   drop_na(value) %>%
#   mutate(month = lubridate::month(time)) %>%
#   filter(month >= 5 & month <= 10) %>%
#   group_by(name,timescale) %>%
#   do(#linearFit = lm(.$value ~ .$spi),
#      pearson = cor( .$drought, .$value,method="pearson")) %>%
#   mutate(#slope = coef(linearFit)[2],
#          pearson_r = pearson[1]) %>%
#          #linearFir_r2 = summary(linearFit)$r.squared,
#          #linearFit_p_value = summary(linearFit)$coefficients[2,4])%>%
#   group_by(name) %>%
#   filter(pearson_r == max(pearson_r)) %>%
#   mutate(depth = parse_number(name),
#          generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
#                                      ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
#                                             ifelse(depth > 20, 'Deep (>20in)', NA))),
#          standardize_method = gsub("[[:digit:]].*$", "", name),
#          standardize_method = str_sub(standardize_method, end=-2)) %>%
#   group_by(generalized_depth, standardize_method) %>%
#   summarise(median_pearson_r = median(pearson_r),
#             timescale_mode = getMode(timescale),
#             timescale_range = max(timescale)- min(timescale)) %>%
#   ungroup() %>%
#   mutate(site_id = ids[i],
#          drought_metric = metric) %>%
#   select(site_id, drought_metric, generalized_depth,
#          standardize_method, median_pearson_r, timescale_mode)
# 
# monthly_correlation = temp_data %>%
#   drop_na(value) %>%
#   mutate(month = lubridate::month(time)) %>%
#   group_by(name,timescale, month) %>%
#   do(#linearFit = lm(.$value ~ .$spi),
#     pearson = cor( .$drought, .$value,method="pearson"),
#     n = length(.$value)) %>%
#   mutate(#slope = coef(linearFit)[2],
#     pearson_r = pearson[1],
#     n = n[1]) %>%
#   #linearFir_r2 = summary(linearFit)$r.squared,
#   #linearFit_p_value = summary(linearFit)$coefficients[2,4])%>%
#   #group_by(name, month) %>%
#   #filter(pearson_r == max(pearson_r)) %>%
#   mutate(depth = parse_number(name),
#          generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
#                                      ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
#                                             ifelse(depth > 20, 'Deep (>20in)', NA))),
#          standardize_method = gsub("[[:digit:]].*$", "", name),
#          standardize_method = str_sub(standardize_method, end=-2)) %>%
#   # group_by(generalized_depth, standardize_method, month) %>%
#   # summarise(median_pearson_r = median(pearson_r),
#   #           timescale_mode = getMode(timescale),
#   #           timescale_range = max(timescale)- min(timescale)) %>%
#   # ungroup() %>%
#   mutate(site_id = ids[i],
#          drought_metric = metric) %>%
#   select(site_id, drought_metric, month, generalized_depth,
#          standardize_method, depth, name, pearson_r, timescale, n)
# 
# tictoc::toc()


library(doSNOW)
tictoc::tic()
cl = makeSOCKcluster(10)
registerDoSNOW(cl)
pb <- txtProgressBar(min=1, max=length(ids), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
out = foreach(i = 1:length(ids), .packages = c('tidyverse'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    temp_data = drought %>%
      filter(site_id == ids[i]) %>%
      left_join(.,  soil_moisture %>%
                  filter(site_id == ids[i]),
                by = c('site_id', 'time'))
    
    if(metric == 'eddi'){
      full_seasonal_corelation = temp_data %>%
        drop_na(value) %>%
        mutate(month = lubridate::month(time)) %>%
        filter(month >= 5 & month <= 10) %>%
        group_by(name,timescale) %>%
        do(#linearFit = lm(.$value ~ .$spi),
          pearson = cor( .$drought, .$value,method="pearson")) %>%
        mutate(#slope = coef(linearFit)[2],
          pearson_r = pearson[1]) %>%
        #linearFir_r2 = summary(linearFit)$r.squared,
        #linearFit_p_value = summary(linearFit)$coefficients[2,4])%>%
        group_by(name) %>%
        filter(pearson_r == min(pearson_r)) %>%
        mutate(depth = parse_number(name),
               generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                           ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                  ifelse(depth > 20, 'Deep (>20in)', NA))),
               standardize_method = gsub("[[:digit:]].*$", "", name),
               standardize_method = str_sub(standardize_method, end=-2)) %>%
        group_by(generalized_depth, standardize_method) %>%
        summarise(median_pearson_r = median(pearson_r),
                  timescale_mode = getMode(timescale),
                  timescale_range = max(timescale)- min(timescale)) %>%
        ungroup() %>%
        mutate(site_id = ids[i],
               drought_metric = metric) %>%
        select(site_id, drought_metric, generalized_depth,
               standardize_method, median_pearson_r, timescale_mode)
    } else {
      full_seasonal_corelation = temp_data %>%
        drop_na(value) %>%
        mutate(month = lubridate::month(time)) %>%
        filter(month >= 5 & month <= 10) %>%
        group_by(name,timescale) %>%
        do(#linearFit = lm(.$value ~ .$spi),
          pearson = cor( .$drought, .$value,method="pearson")) %>%
        mutate(#slope = coef(linearFit)[2],
          pearson_r = pearson[1]) %>%
        #linearFir_r2 = summary(linearFit)$r.squared,
        #linearFit_p_value = summary(linearFit)$coefficients[2,4])%>%
        group_by(name) %>%
        filter(pearson_r == max(pearson_r)) %>%
        mutate(depth = parse_number(name),
               generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                           ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                  ifelse(depth > 20, 'Deep (>20in)', NA))),
               standardize_method = gsub("[[:digit:]].*$", "", name),
               standardize_method = str_sub(standardize_method, end=-2)) %>%
        group_by(generalized_depth, standardize_method) %>%
        summarise(median_pearson_r = median(pearson_r),
                  timescale_mode = getMode(timescale),
                  timescale_range = max(timescale)- min(timescale)) %>%
        ungroup() %>%
        mutate(site_id = ids[i],
               drought_metric = metric) %>%
        select(site_id, drought_metric, generalized_depth,
               standardize_method, median_pearson_r, timescale_mode)
    }
    
    
    
    monthly_correlation = temp_data %>%
      drop_na(value) %>%
      mutate(month = lubridate::month(time)) %>%
      group_by(name,timescale, month) %>%
      do(#linearFit = lm(.$value ~ .$spi),
        pearson = cor( .$drought, .$value,method="pearson"),
        n = length(.$value)) %>%
      mutate(#slope = coef(linearFit)[2],
        pearson_r = pearson[1],
        n = n[1]) %>%
      #linearFir_r2 = summary(linearFit)$r.squared,
      #linearFit_p_value = summary(linearFit)$coefficients[2,4])%>%
      #group_by(name, month) %>%
      #filter(pearson_r == max(pearson_r)) %>%
      mutate(depth = parse_number(name),
             generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                         ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                ifelse(depth > 20, 'Deep (>20in)', NA))),
             standardize_method = gsub("[[:digit:]].*$", "", name),
             standardize_method = str_sub(standardize_method, end=-2)) %>%
      # group_by(generalized_depth, standardize_method, month) %>%
      # summarise(median_pearson_r = median(pearson_r),
      #           timescale_mode = getMode(timescale),
      #           timescale_range = max(timescale)- min(timescale)) %>%
      # ungroup() %>%
      mutate(site_id = ids[i],
             drought_metric = metric) %>%
      select(site_id, drought_metric, month, generalized_depth,
             standardize_method, depth, name, pearson_r, timescale, n)
    
    export = list(full_seasonal_corelation, monthly_correlation)
    export
  }, error = function(e){
    export = NA
    export
  })
  
} 
close(pb)
stopCluster(cl)
tictoc::toc()

saveRDS(out, paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlations/',metric,'-cor.RDS'))