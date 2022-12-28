library(tidyverse)
library(Hmisc)
library(magrittr)
library(data.table)

`%notin%` = Negate(`%in%`)

#site_of_interest = sites[178]

vwc_cor = function(site_of_interest, vwc, drought, min_n){
  #filter vwc data for site of interest
  vwc_filtered_partial = vwc %>%
    filter(site_id == site_of_interest)
  
  #filter drought for site of interest and time in vwc
  drought_filtered = drought %>%
    filter(site_id == site_of_interest,
           time %in% vwc_filtered_partial$date)
  
  #filter VWC for valid drought timescales (not enough data to calculate
  #drought metrics earlier than 2010 - 30 years required)
  vwc_filtered = vwc_filtered_partial %>%
    filter(date %in% drought_filtered$time)
  
  #compute valid depths
  vwc_cols = colnames(vwc_filtered) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'date')) %>%
    mutate(depth = parse_number(value))
  
  #compute valid timescales
  time_scales = colnames(drought_filtered) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'time')) 
  
  #check for time continuity
  if(all.equal(vwc_filtered$date, drought_filtered$time)){
    #compute correlation matrix
    cor_matrix = bind_cols(vwc_filtered %>% 
                             dplyr::select(-c(site_id, date)),
                           drought_filtered %>% 
                             dplyr::select(-c(site_id, time))) %>%
      na_if(., Inf) %>%
      na_if(., -Inf) %>%
      as.matrix() %>%
      rcorr(type = "pearson")
    
    # find significant P values of correlation matrix
    P_filter = ifelse(cor_matrix$P < 0.05, 1, NA)
    # make sure min number of observations is met
    n_filter = ifelse(cor_matrix$n > min_n, 1, NA)
    
    #extract pearson correlation matrix
    r_matrix = (cor_matrix %$%
      #filter for significance and min n            
      r * P_filter * n_filter) %>%
      as.data.frame() %>%
      rownames_to_column(var = "timescale") %>%
      as_tibble() %>%
      filter(timescale %in% time_scales$value) %>%
      dplyr::select(c(timescale, vwc_cols$value)) %>%
      mutate(site_id = site_of_interest) %>%
      relocate(site_id)
      
    return(r_matrix)
  }
  return(NA)
}

#significance and n test
p_n_filter = function(x, min_n_seasonal, time_scales, vwc_cols, site_of_interest){
  x_ = x[[1]]
  P_filter = ifelse(x_$P < 0.05, 1, NA)
  # make sure min number of observations is met
  n_filter = ifelse(x_$n > min_n_seasonal, 1, NA)
  
  #extract pearson correlation matrix
  r_matrix = (x_ %$%
                #filter for significance and min n            
                r * P_filter * n_filter)%>%
      as.data.frame() %>%
      rownames_to_column(var = "timescale") %>%
      as_tibble() %>%
      filter(timescale %in% time_scales$value) %>%
      dplyr::select(c(timescale, vwc_cols$value)) %>%
      mutate(site_id = site_of_interest) %>%
      relocate(site_id)
  
  return(list(r_matrix))
}

# sites[i], vwc, drought, min_n_seasonal
# site_of_interest = sites[i]
vwc_cor_seasonal = function(site_of_interest, vwc, drought, min_n_seasonal){
  #filter vwc data for site of interest
  vwc_filtered_partial = vwc %>%
    filter(site_id == site_of_interest)
  
  #filter drought for site of interest and time in vwc
  drought_filtered = drought %>%
    filter(site_id == site_of_interest,
           time %in% vwc_filtered_partial$date) %>%
    mutate(date = time) %>%
    relocate(date) %>%
    dplyr::select(-c('time'))
  
  #filter VWC for valid drought timescales (not enough data to calculate
  #drought metrics earlier than 2010 - 30 years required)
  vwc_filtered = vwc_filtered_partial %>%
    filter(date %in% drought_filtered$date) %>%
    select_if(~sum(!is.na(.)) >= 5)
  
  #compute valid depths
  vwc_cols = colnames(vwc_filtered) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'date')) %>%
    mutate(depth = parse_number(value))
  
  #compute valid timescales
  time_scales = colnames(drought_filtered) %>%
    as_tibble() %>%
    filter(value %notin% c('site_id', 'time')) 
  
  #check for time continuity
  if(all.equal(vwc_filtered$date, drought_filtered$date)){
    #compute valid months (4 obs needed for function, then we further filter from there)
    months_ = left_join(vwc_filtered, drought_filtered, by = c('site_id', 'date')) %>%
      mutate(month = lubridate::month(date)) %>%
      group_by(month) %>%
      summarise(n = length(date)) %>%
      filter(n > 4)
    
    #compute correlation matrix
    cor_matrix = left_join(vwc_filtered, drought_filtered, by = c('site_id', 'date')) %>%
      mutate(month = lubridate::month(date)) %>%
      filter(month %in% months_$month) %>%
      group_by(month) %>%
      na_if(., Inf) %>%
      na_if(., -Inf) %>%
      dplyr::select(-c(site_id, date)) %>%
      do(cor = as.matrix(.) %>%
           rcorr(., type = "pearson")) %>%
      dplyr::select(cor) %>%
      do(cor_filtered = p_n_filter(., min_n_seasonal, time_scales, vwc_cols, site_of_interest)) %>%
      ungroup() %>%
      mutate(month  = months_$month,
             site_id = site_of_interest) %>%
      relocate(site_id, month)
    
    return(cor_matrix)
  }
  return(NA)
}
