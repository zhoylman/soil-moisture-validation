#download, bind and and write out gridmet data as nc files for extraction

#load libs
library(tidyverse)
library(raster)
library(sf)

# pet
year = 1979:2021
url = paste0('https://www.northwestknowledge.net/metdata/data/pet_',year,'.nc')
dest = paste0('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/pet_raw/pet_',year,'.nc')
for(i in 1:length(url)){
  download.file(url[i], dest[i])
}

# use CDO to merge Climate Data Operator (CDO) command line tool
# https://www.isimip.org/protocol/preparing-simulation-files/cdo-help/

command = paste0('cdo mergetime', paste0(' ', dest, collapse = ''), 
                 ' /home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pet.nc',
                 collapse = '')

system(command)

#compute time 
pet = brick('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pet.nc')
time = data.frame(datetime = as.Date(as.numeric(substring(names(pet),2)), origin="1900-01-01"))
readr::write_csv(time, '/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pet_time.csv')

## precip
url = paste0('https://www.northwestknowledge.net/metdata/data/pr_',year,'.nc')
dest = paste0('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/pr_raw/pr_',year,'.nc')
for(i in 1:length(url)){
  download.file(url[i], dest[i])
}

# use CDO to merge Climate Data Operator (CDO) command line tool
# https://www.isimip.org/protocol/preparing-simulation-files/cdo-help/

command = paste0('cdo mergetime', paste0(' ', dest, collapse = ''), 
                 ' /home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pr.nc',
                 collapse = '')

system(command)

#compute time 
pr = brick('/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pr.nc')
time = data.frame(datetime = as.Date(as.numeric(substring(names(pr),2)), origin="1900-01-01"))
readr::write_csv(time, '/home/zhoylman/soil-moisture-validation-data/raw/gridmet/gridmet_pr_time.csv')