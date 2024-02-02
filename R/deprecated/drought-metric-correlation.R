library(tidyverse)
library(doSNOW)
source('/home/zhoylman/soil-moisture-validation/R/correlation-functions.R')

drought_metric = 'eddi'

#function to remove lists with only NA
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

vwc = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/standardized-soil-moisture/standardized-soil-moisture-data-wide.csv')
drought = read_csv(paste0('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/', drought_metric,'-data-wide-10s.csv'))
sites = unique(vwc$site_id)

min_n = 365
min_n_seasonal = 60

#Full time period correlation
tictoc::tic()
cl = makeSOCKcluster(20)
registerDoSNOW(cl)
pb = txtProgressBar(min=1, max=length(sites), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
#length(ids)
out = foreach(i = 1:length(sites), .packages = c('tidyverse', 'lubridate', 'magrittr', 'data.table', 'Hmisc'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    return = vwc_cor(sites[i], vwc, drought, min_n)
    return
  }, error = function(e){
    return = NA
    return
  })
} 
close(pb)
stopCluster(cl)
tictoc::toc()

final = out %>%
  na.omit.list() %>%
  bind_rows()

write_csv(final, paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/', drought_metric,'-correlation-matrix.csv'))

### seasonal correlation

tictoc::tic()
cl = makeSOCKcluster(20)
registerDoSNOW(cl)
pb = txtProgressBar(min=1, max=length(sites), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
#length(ids)
out_seasonal = foreach(i = 1:length(sites), .packages = c('tidyverse', 'lubridate', 'magrittr', 'data.table', 'Hmisc'), .options.snow=opts) %dopar% {
  gc()
  tryCatch({
    return = vwc_cor_seasonal(sites[i], vwc, drought, min_n_seasonal)
    return
  }, error = function(e){
    return = NA
    return
  })
} 
close(pb)
stopCluster(cl)
tictoc::toc()

final_seasonal = out_seasonal %>%
  na.omit.list() %>%
  bind_rows()

saveRDS(final_seasonal, paste0('/home/zhoylman/soil-moisture-validation-data/processed/correlation-matrix/', drought_metric, '-correlation-matrix-seasonal.RDS'))
