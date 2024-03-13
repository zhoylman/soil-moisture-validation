#download, bind and and write out gridmet data as nc files for point extraction

#load libs
library(tidyverse)
library(raster)
library(sf)

# define the years of interest
year = 1979:2021

#define download-merge function
download_gridmet = function(variable){
  system(paste0('mkdir /mnt/data1/gridmet/', variable, '_raw'))
  ## define urls
  url = paste0('https://www.northwestknowledge.net/metdata/data/', variable,'_',year,'.nc')
  dest = paste0('/mnt/data1/gridmet/',variable,'_raw/',variable,'_',year,'.nc')
  #download raw data
  for(i in 1:length(url)){
    download.file(url[i], dest[i])
  }
  # use CDO to merge Climate Data Operator (CDO) command line tool
  # https://www.isimip.org/protocol/preparing-simulation-files/cdo-help/
  command = paste0('cdo mergetime', paste0(' ', dest, collapse = ''),
                   paste0(' /mnt/data1/gridmet/gridmet_',variable,'.nc'),
                   collapse = '')

  system(command)
  #compute time 
  data = brick(paste0('/mnt/data1/gridmet/gridmet_',variable,'.nc'))
  time = data.frame(datetime = as.Date(as.numeric(substring(names(data),2)), origin="1900-01-01"))
  readr::write_csv(time, paste0('/mnt/data1/gridmet/gridmet_',variable,'_time.csv'))
}

# pr, pet
download_gridmet('pr') 