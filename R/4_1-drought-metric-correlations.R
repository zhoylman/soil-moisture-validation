#load libs  
library(tidyverse)
library(scales)
library(doSNOW)

#define get mode function
getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#define drought metirc to process
metric = 'eddi'
metric_upper = casefold(metric, upper = T)

#import drought data
drought = read_csv(paste0('~/soil-moisture-validation-data/processed/drought-metrics/',metric,'-data-long-10s.csv'))%>%
  rename(drought = metric) %>%
  #clamp data at -2 and 2
  mutate(drought = ifelse(drought > 2, 2, drought),
         drought = ifelse(drought < -2, -2, drought))

#import soil moisture data
soil_moisture = read_csv('~/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide-6-years-min-CDF-w-mean.csv') %>%
  pivot_longer(., cols = -c(site_id,date)) %>%
  mutate(time = date) %>%
  select(site_id,time,name,value) %>%
  #convert anomalies to standard normal if they are anomolies but leave raw alone
  mutate(value_new = ifelse(grepl("anomaly", name),  qnorm(value), value),
         #clamp data at -2 and 2
         value_new = ifelse(value_new > 2, 2, value_new),
         value_new = ifelse(value_new < -2, -2, value_new),
         value = ifelse(grepl("anomaly", name),  value_new, value)) %>%
  select(site_id,time,name,value)

#compute ids to process (site locations)
ids = unique(soil_moisture$site_id)

#set up time accounting
tictoc::tic()

#set up the multi processing scheme and cluster
cl = makeSOCKcluster(10)
registerDoSNOW(cl)

#define parameters for the 
progress_bar = function(i){
  # <<- defines variable to global enviorment
  pb <<- txtProgressBar(min=1, max=length(i), style=3)
  progress <<- function(n) setTxtProgressBar(pb, n)
  opts <<- list(progress=progress)
}

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
#out = foreach(i = c(1,100,200,300,400,500,600,700), .packages = c('tidyverse'), .options.snow= progress_bar(ids)) %dopar% {
out = foreach(i = 1:length(ids), .packages = c('tidyverse'), .options.snow= progress_bar(ids)) %dopar% {
  #clear memory
  gc()
  #set up basic error handling
  tryCatch({
    #define dataset of interest
    temp_data = drought %>%
      #filter data by id of interest
      filter(site_id == ids[i]) %>%
      #left join with soil moisture data
      left_join(.,  soil_moisture %>%
                  filter(site_id == ids[i]),
                #by site id and time
                by = c('site_id', 'time'))
    
    #if statement for eddi vs spi and spei becasue of directionality change
    #in correlation (negative for eddi, positive for spi, spei)
    if(metric == 'eddi'){
      #compute full season correlation (May 1st - Oct 31st) 
      #SPI/SPEI doesn't make sense with snow! (rationale)
      full_seasonal_corelation = temp_data %>%
        #convert Inf and -Inf to NA
        mutate(drought = ifelse(is.infinite(drought), NA, drought),
               value = ifelse(is.infinite(value), NA, value)) %>%
        #drop NAs in soil moisture responses
        drop_na(value) %>%
        #drop NAs in drought metrics
        drop_na(drought) %>%
        #compute month id
        mutate(month = lubridate::month(time)) %>%
        #filter for the period of interest
        filter(month >= 5 & month <= 10) %>%
        #define the group by scheme (name = soil moisture depth, timescale)
        group_by(name, timescale) %>%
        #compute correlation 2 ways, as a linear model and as a pearson cor
        #linear model will be used to compute RMSE
        do(
          #compute linear model
          #linearFit = lm(.$value ~ .$drought),
          #compute pearson r
          pearson = cor(.$drought, .$value,method="pearson"),
          #compute number of obs in model 
          n = length(.$value)
          ) %>%
        #extract information from each model
        mutate(
          #extract number of obs in model
          n = unlist(n),
          #extract pearson r
          pearson_r = unlist(pearson)
          #extract r2 from lm
          #linearFit_r2 = summary(linearFit)$r.squared,
          #compute RMSE from the lm
          #linearFit_RMSE = RMSE(linearFit)
          ) %>%
        #filter out no data models
        filter(n > 1) %>%
        #group by soil moisture depth and standardization method
        group_by(name) %>%
        #filter for the minimum person r (this where the eddi vs spi functions differ)
        filter(pearson_r == min(pearson_r)) %>%
        #compute some information about depth
        #compute some information about depth
        mutate(depth = parse_number(name),
               #compute the generalized depth
               generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                           ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                  ifelse(depth > 20, 'Deep (>20in)', NA))),
               #add a 'Mean soil moisture' generalized depth
               generalized_depth = ifelse(str_detect(name, 'mean'), 'Depth Averaged', generalized_depth),
               #compute the generalized method (drought, storage, raw)
               standardize_method = gsub("[[:digit:]].*$", "", name),
               standardize_method = str_sub(standardize_method, end=-2),
               standardize_method = ifelse(name == 'mean_soil_moisture', 'soil_moisture', standardize_method),
               standardize_method = ifelse(name == 'drought_anomaly_mean', 'drought_anomaly', standardize_method),
               standardize_method = ifelse(name == 'storage_anomaly_mean', 'storage_anomaly', standardize_method)) %>%
        #group_by the genrealized depth
        group_by(generalized_depth, standardize_method) %>%
        #compute summaries by generalized depth and method
        summarise(median_pearson_r = median(pearson_r),
                  #median_linearFit_r2 = median(linearFit_r2),
                  #median_linearFit_RMSE = median(linearFit_RMSE),
                  timescale_mode = getMode(timescale),
                  timescale_range = max(timescale)- min(timescale)) %>%
        #ungroup
        ungroup() %>%
        #append the site ID back into dataset
        mutate(site_id = ids[i],
               #rename drought metric
               drought_metric = metric) %>%
        #select final data
        select(site_id, drought_metric, generalized_depth,
               standardize_method, median_pearson_r, #median_linearFit_r2, median_linearFit_RMSE, 
               timescale_mode)
    } else {
      #compute full season correlation (May 1st - Oct 31st) 
      #SPI/SPEI doesn't make sense with snow! (rationale)
      full_seasonal_corelation = temp_data %>%
        #convert Inf and -Inf to NA
        mutate(drought = ifelse(is.infinite(drought), NA, drought),
               value = ifelse(is.infinite(value), NA, value)) %>%
        #drop NAs in soil moisture responses
        drop_na(value) %>%
        #drop NAs in drought metrics
        drop_na(drought) %>%
        #compute month id
        mutate(month = lubridate::month(time)) %>%
        #filter for the period of interest
        filter(month >= 5 & month <= 10) %>%
        #define the group by scheme (name = soil moisture depth, timescale)
        group_by(name, timescale) %>%
        #compute correlation 2 ways, as a linear model and as a pearson cor
        #linear model will be used to compute RMSE
        do(
          #compute linear model
          #linearFit = lm(.$value ~ .$drought),
          #compute pearson r
          pearson = cor(.$drought, .$value,method="pearson"),
          #compute number of obs in model 
          n = length(.$value)
        ) %>%
        #extract information from each model
        mutate(
          #extract number of obs in model
          n = unlist(n),
          #extract pearson r
          pearson_r = unlist(pearson)
          #extract r2 from lm
          #linearFit_r2 = summary(linearFit)$r.squared,
          #compute RMSE from the lm
          #linearFit_RMSE = RMSE(linearFit)
        ) %>%
        #filter out no data models
        filter(n > 1) %>%
        #group by soil moisture depth and standardization method
        group_by(name) %>%
        #filter for the minimum person r (this where the eddi vs spi functions differ)
        filter(pearson_r == max(pearson_r)) %>%
        #compute some information about depth
        mutate(depth = parse_number(name),
               #compute the generalized depth
               generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                           ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                  ifelse(depth > 20, 'Deep (>20in)', NA))),
               #add a 'Mean soil moisture' generalized depth
               generalized_depth = ifelse(str_detect(name, 'mean'), 'Depth Averaged', generalized_depth),
               #compute the generalized method (drought, storage, raw)
               standardize_method = gsub("[[:digit:]].*$", "", name),
               standardize_method = str_sub(standardize_method, end=-2),
               standardize_method = ifelse(name == 'mean_soil_moisture', 'soil_moisture', standardize_method),
               standardize_method = ifelse(name == 'drought_anomaly_mean', 'drought_anomaly', standardize_method),
               standardize_method = ifelse(name == 'storage_anomaly_mean', 'storage_anomaly', standardize_method)) %>%
        #group_by the genrealized depth
        group_by(generalized_depth, standardize_method) %>%
        #compute summaries by generalized depth and method
        summarise(median_pearson_r = median(pearson_r),
                  #median_linearFit_r2 = median(linearFit_r2),
                  #median_linearFit_RMSE = median(linearFit_RMSE),
                  timescale_mode = getMode(timescale),
                  timescale_range = max(timescale) - min(timescale)) %>%
        #ungroup
        ungroup() %>%
        #append the site ID back into dataset
        mutate(site_id = ids[i],
               #rename drought metric
               drought_metric = metric) %>%
        #select final data
        select(site_id, drought_metric, generalized_depth,
               standardize_method, median_pearson_r, #median_linearFit_r2, median_linearFit_RMSE, 
               timescale_mode)
    }
    
    #compute the monthly correlations
    monthly_correlation = temp_data %>%
      #convert Inf and -Inf to NA
      mutate(drought = ifelse(is.infinite(drought), NA, drought),
             value = ifelse(is.infinite(value), NA, value)) %>%
      #drop NAs in soil moisture responses
      drop_na(value) %>%
      #drop NAs in drought metrics
      drop_na(drought) %>%
      #compute month identifier
      mutate(month = lubridate::month(time)) %>%
      #group_by name (soil moisture data), timescale and month
      group_by(name, timescale, month) %>%
      #compute per group data
      do(
        #compute the linear fit
        #linearFit = lm(.$value ~ .$drought),
        #compute the pearson correlation
        pearson = cor( .$drought, .$value,method="pearson"),
        #compute the number of obs in each model based on soil moisture obs
        n = length(.$value)
        ) %>%
      #extract some information from each model
      mutate(
        #extract pearsons r
        pearson_r = unlist(pearson),
        #extract number of obs
        n = unlist(n)
        #extract r2 from lm
        #linearFit_r2 = summary(linearFit)$r.squared,
        #compute RMSE from the lm
        #linearFit_RMSE = RMSE(linearFit)
        ) %>%
      #compute some information about depth
      mutate(depth = parse_number(name),
             #compute the generalized depth
             generalized_depth = ifelse( depth <= 4, 'Shallow (0-4in)',
                                         ifelse(depth >= 8 & depth <= 20, 'Middle (8in - 20in)',
                                                ifelse(depth > 20, 'Deep (>20in)', NA))),
             #add a 'Mean soil moisture' generalized depth
             generalized_depth = ifelse(str_detect(name, 'mean'), 'Depth Averaged', generalized_depth),
             #compute the generalized method (drought, storage, raw)
             standardize_method = gsub("[[:digit:]].*$", "", name),
             standardize_method = str_sub(standardize_method, end=-2),
             standardize_method = ifelse(name == 'mean_soil_moisture', 'soil_moisture', standardize_method),
             standardize_method = ifelse(name == 'drought_anomaly_mean', 'drought_anomaly', standardize_method),
             standardize_method = ifelse(name == 'storage_anomaly_mean', 'storage_anomaly', standardize_method))  %>%
      #add back site identifier
      mutate(site_id = ids[i],
             #rename drought metric name
             drought_metric = metric) %>%
      #select variables of interest
      select(site_id, drought_metric, month, generalized_depth,
             standardize_method, depth, name, pearson_r, 
             #linearFit_r2, linearFit_RMSE, 
             timescale, n)
    
    #define export list
    export = list(full_seasonal_corelation, monthly_correlation)
    #push outside of foreach function
    rm(temp_data, full_seasonal_corelation, monthly_correlation)
    gc(); gc()
    return(export)
    rm(export)
    gc(); gc()
  }, error = function(e){
    #if there is an error, return NA
    export = NA
    #push outside of foreach function
    return(export)
    gc(); gc()
  })
} 
#close the progress bar
close(pb)
#stop cluster
stopCluster(cl)
#compute run time
tictoc::toc()

#save out final data as RDS dataset
saveRDS(out, paste0('~/soil-moisture-validation-data/processed/correlations/',metric,'-6-years-min-cor-rmse-clamped-w-mean.RDS'))

#memory management
gc(); gc()
