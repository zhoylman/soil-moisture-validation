library(tictoc)
library(foreach)
library(doParallel)
library(parallel)
library(dplyr)
library(data.table)
library(sf)
library(stringr)
library(magrittr)
library(spData)
library(tidyverse)
library(rlist)
library(xml2)
library(snow)
library(doSNOW)
library(progress)

#define special
`%notin%` = Negate(`%in%`)

#get us domain
data(us_states)

#compute conus
conus = us_states %>%
  filter(NAME %notin% c('Alaska', 'Hawaii')) %>%
  st_transform(st_crs(4326))

#define dunctions for getting NRCS metadata
grabNRCS.meta<-function(ntwrks="SCAN",cnvrt.elev=FALSE){
  options(timeout = 120)
  
  #QC check:
  ntwrks.cntrl<-c("SCAN","SNTL","SNTLT","SNOW","MPRC","OTHER","COOP","USGS","MSNT","BOR","CLMIND","ALL")
  ntwrks.QC<-ntwrks %in% ntwrks.cntrl
  if(any(ntwrks.QC)==FALSE){
    stop("Incorrect ntwrk name(s). Please use the controlled terms: 'SCAN','SNTL','SNTLT','SNOW','MPRC','OTHER','COOP','USGS','MSNT','BOR','CLMIND', or you may use 'ALL' to pull metadata for all NRCS networks")
  }
  if(ntwrks[1]=="ALL"){
    #grab all networks from controlled list:
    ntwrks<-ntwrks.cntrl[1:(length(ntwrks.cntrl)-1)]
  }
  #ntwrks<-"SCAN"
  ntwrks.Lcase<-tolower(ntwrks)
  #create dynamic url link using "ntwrks" vector to access the data we want:
  urls<-lapply(ntwrks.Lcase,function(x) paste0("https://wcc.sc.egov.usda.gov/nwcc/yearcount?network=",x,"&counttype=listwithdiscontinued&state="))
  #grab metadata from each of the sites:
  stationMetadata<-lapply(urls, function(x) xml2::read_html(x, ))
  
  #grab headers:
  headers.NRCS<-lapply(stationMetadata, function(x) x %>% rvest::html_nodes("th"))
  #clean the headers of HTML code:
  metaHeaders<-lapply(headers.NRCS, function(x) gsub("<.*?>", "", x))
  #grab the actual metadata from the urls:
  metadata.NRCS<-lapply(stationMetadata, function(x) x %>% rvest::html_nodes("td"))
  #find the first entry of metadata using the NRCS station type
  #searchNtwks<-paste0("\\<td>")
  searchNtwks<-paste(paste0("\\<td>",ntwrks),collapse="|")
  findLineStart<-lapply(metadata.NRCS, function(x)  grep(searchNtwks,x))
  #clean the metadata of HTML tags:
  cleanMeta<-lapply(metadata.NRCS,function(x) gsub("<.*?>", "", x))
  #loop thru and finish the cleanup and conversion to DF:
  NRCS.metadata<-list()
  for(i in 1:length(cleanMeta)){
    #subset to only keep metadata from table:
    startIndex<-findLineStart[[i]][1]
    endIndex<-findLineStart[[i]][length(findLineStart[[i]])]+length(metaHeaders[[i]])-1
    cleanMeta.sub<-cleanMeta[[i]][startIndex:endIndex]
    #split data into equal chunks
    metadata.NRCS.clean<-split(cleanMeta.sub, ceiling(seq_along(cleanMeta.sub)/length(metaHeaders[[i]])))
    #put data into df:
    NRCS.metadata[[i]]<-data.frame(do.call(rbind,metadata.NRCS.clean))
    #add headers as names:
    names(NRCS.metadata[[i]])<-metaHeaders[[i]]
    #remove excess white space within the ntwrks column:
    NRCS.metadata[[i]]$ntwk<-trimws(NRCS.metadata[[i]]$ntwk,"both")
    #partition the site name and site ID:
    NRCS.metadata[[i]]$site_id<-paste0(NRCS.metadata[[i]]$ntwk,":",gsub(".*\\((.*)\\).*", "\\1", NRCS.metadata[[i]]$site_name))
    #add the ntwrk name to the siteID for more accurate reference:
    #NRCS.metadata[[i]]$site_id
    #remove the site id from the "site.name" column
    NRCS.metadata[[i]]$site_name<-gsub("\\s*\\([^\\)]+\\)","",as.character(NRCS.metadata[[i]]$site_name))
    #convert feet to meters:
    if(cnvrt.elev==TRUE){
      #convert elevation column
      NRCS.metadata[[i]]$elev<-as.numeric(levels(unlist(NRCS.metadata[[i]]$elev)))[unlist(NRCS.metadata[[i]]$elev)]*0.3048
      #rename to elev_m
      names(NRCS.metadata[[i]])[names(NRCS.metadata[[i]]) == 'elev'] <- 'elev_m'
    }
    else{
      #change the 'elev' column name to 'elev_ft' to make avoid ambiguity
      names(NRCS.metadata[[i]])[names(NRCS.metadata[[i]]) == 'elev'] <- 'elev_ft'
    }
    
  }
  #set the names of the dataframes to the corresponding network name:
  names(NRCS.metadata)<-ntwrks
  return(NRCS.metadata)#output the metadata to user
  
}#end Function

#define dunctions for getting NRCS data
grabNRCS.elements<-function(site_id="SCAN:2221"){
  #grab full site_id code (which includes network) to assign at the end
  site_id.label<-site_id
  #site.QC.check: check to see if site_id is from one of the accepted networks:
  ntwrks.unique<-unique(gsub(":.*","",site_id))
  ntwrks.cntrl<-c("SCAN","SNTL","SNTLT","OTHER")
  ntwrks.QC<-ntwrks.unique %in% ntwrks.cntrl
  if(any(ntwrks.QC)==FALSE){
    stop("Unable to pull metadata. Metadata can only be pulled for the following NRCS networks: 'SCAN','SNTL','SNTLT', and 'OTHER'")
  }
  #strip any "ntwrk" name out of the site_ids if reading in from other function:
  site_id<-gsub(".*:","",site_id)
  #dynamically create URLs:
  urls<-lapply(site_id,function(x) paste0("https://wcc.sc.egov.usda.gov/nwcc/site?sitenum=",x))
  #grab element level metadata from each of the sites:
  site.elements<-lapply(urls, function(x) xml2::read_html(x))
  #define %>% globally or it won't work:
  `%>%` <-magrittr::`%>%`
  site.elements.sub<-lapply(site.elements, function(x) x %>% rvest::html_nodes(xpath='//*[@id="report"]/option'))
  #clean the html tags:
  elements.clean<-lapply(site.elements.sub, function(x) gsub("<.*?>", "", x))
  #assign names now - can use this later to properly name the resulting list of sites:
  names(elements.clean)<-site_id.label
  #find sites with no element level metadata:
  noElementMeta.sites<-which(lapply(elements.clean, function(x) length(x))==0)
  #remove sites with  no element level metadata:
  if(length(noElementMeta.sites)>0){
    elements.clean<-elements.clean[-noElementMeta.sites]
    if(length(elements.clean)==0){
      stop("No element level metadata for site of interest.")
    }
  }
  #subset the data based on grep
  #elementSubset<-lapply(elements.clean, function(x) x[(grep("Individual elements",x)+1):(grep("Daily",x)-1)])
  #clean the data elements:
  variables.init<-lapply(elements.clean,function(x) gsub("[[:digit:]]|\\(|\\)","",x))
  #a bit more cleaning:
  variables<-lapply(variables.init, function(x) substr(x,1,nchar(x)-3))
  #grab the dates:
  #dates<-lapply(elementSubset, function(x) gsub(".*\\((.*)\\).*", "\\1", x))
  #put information into dataframe and output:
  elementData<-lapply(Map(cbind,variables), function(x) data.frame(x))
  #set names within each dataframe:
  #elementData<-lapply(elementData, stats::setNames, nm = c("element","date.start"))
  #final Check to make sure element level data and dates actually exist:
  QC.check<-lapply(elements.clean, function(x) ((grep("Individual elements",x)+1)-(grep("Daily",x)-1)))
  QC.check.index<-which(unlist(QC.check)>0)  #this checks to see if an NRCS site doesn't exist or doesn't have element data
  #replace those that fail QC test with NA
  elementData[QC.check.index]<-NA
  #set the names of each nested dataframe to the NRCS station id:
  names(elementData)<-names(elements.clean)
  #return the list of dataframes comprising element level data:
  return(elementData)
}#end function

#pull all the meta
nrcs_sites_raw = grabNRCS.meta(ntwrks=c("ALL"))

#sort and filter meta for networks of interest
nrcs_sites = nrcs_sites_raw %>%
  bind_rows() %>%
  as_tibble() %>%
  filter(ntwk %in% c('SCAN', 'SNTL', 'SNTLT', 'OTHER'))

#rev up a cluster for parallel downloads
cl = makeCluster(c(detectCores() - 22), outfile = '/home/zhoylman/Desktop/log.txt')
registerDoParallel(cl)

#get list of all elements assosiated with each station (soil moisture, temp, etc)
nrcs_elements = foreach(i = 1:length(nrcs_sites$site_id), .packages = c('dplyr', 'magrittr')) %dopar% {
  tryCatch({
    grabNRCS.elements(nrcs_sites$site_id[i]) %$%
      get(nrcs_sites$site_id[i]) %>%
      mutate(site_id = nrcs_sites$site_id[i])
  }, error = function(e){
    return(data.frame(element = NA, date.start = NA, site_id = NA))
  })
} %>%
  bind_rows()

#determine which sites are valid (have the data we want)
valid_nrcs_sites = nrcs_elements %>%
  #remove missing data
  tidyr::drop_na() %>%
  #first filter for soil moisture data presence
  filter(x %in% c('Soil Moisture Percent', 'Soil Moisture Percent ')) %>%
  #group_by id and assign boolean and time
  group_by(site_id) %>%
  summarise(soil_moisture = T) %>%
  #now the same thing but for temperaure
  left_join(.,  nrcs_elements %>%
              tidyr::drop_na() %>%
              filter(x %in% c('Soil Temperature', 'Soil Temperature ')) %>%
              group_by(site_id) %>%
              summarise(soil_temperature = T), by = c('site_id')) %>%
  #join datasets
  left_join(., nrcs_sites %>% 
              #map data to numeric form
              select(site_id, latitude, longitude, elev_ft, ntwk, state) %>%
              mutate(latitude = latitude %>% as.numeric(),
                     longitude = longitude %>% as.numeric(),
                     elev_ft = elev_ft %>% as.numeric()),  by = 'site_id') %>%
  select(site_id, ntwk, state, latitude, longitude, elev_ft, soil_moisture, soil_temperature) %>%
  #convert to spatial
  st_as_sf(., coords = c('longitude', 'latitude')) %>%
  st_set_crs(st_crs(4326)) %>%
  #select CONUS stations
  st_intersection(., conus) %>%
  #filter for both soil moisture and soil temperature
  filter(soil_moisture == T & soil_temperature == T)
#stop cluster (memory control)
stopCluster(cl)

#it appears that NRCS limits to 3 simultaneous requests?? 
#rev up a new cluster that conforms to NRCS limits
cl = makeCluster(3)
registerDoSNOW(cl)
iterations = length(valid_nrcs_sites$site_id)

#define the parameters used for the progress bar
pb = progress_bar$new(
  format = "percent = :letter [:bar] :elapsed | eta: :eta",
  total = iterations,    # 100 
  width = 60)
progress_letter = (((1:length(valid_nrcs_sites$site_id)) /length(valid_nrcs_sites$site_id))*100) %>% round(., 2) # token reported in progress bar

# allowing progress bar to be used in foreach -----------------------------
progress = function(n){
  pb$tick(tokens = list(letter = progress_letter[n]))
} 

opts = list(progress = progress)

#pull data from valid sites
nrcs_data = foreach(i = 1:length(valid_nrcs_sites$site_id),  
                    .packages = c('tidyverse', 'dplyr', 'magrittr', 'readr', 'stringr', 'sf', 'data.table','tidyr'),
                    .options.snow = opts) %dopar%{
  tryCatch({
    #GUI builder https://wcc.sc.egov.usda.gov/reportGenerator/
    #define the urls used to pull data
    baseURL = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultiTimeSeriesGroupByStationReport/daily/start_of_period/"
    endURL = "/POR_BEGIN,POR_END/SMV:-2:value,SMV:-4:value,SMV:-8:value,SMV:-20:value,SMV:-40:value,STV:-2:value,STV:-4:value,STV:-8:value,STV:-20:value,STV:-40:value,PRCP::value,TMAX::value,TMIN::value"
    fullURL = paste0(baseURL, valid_nrcs_sites$site_id[i] %>% parse_number(), ':', valid_nrcs_sites$state[i], ':', valid_nrcs_sites$ntwk[i], endURL)
    #get the raw data and parse
    temp = RCurl::getURL(fullURL) %>%
      read_csv(., skip = 69)
    
    #compute the name conversion
    name_conversion_temp = tibble(names = names(temp),
                                  names_wo_id = names(temp) %>% str_remove(., valid_nrcs_sites$site_id[i] %>% parse_number() %>% as.character())) %>%
      mutate(partial_element = ifelse(names %like% 'Soil Temperature', 'soil_temperature_', 
                                      ifelse(names %like% 'Moisture', 'soil_moisture_', 
                                             ifelse(names %like% 'Precipitation Increment', 'precipitation_increment_in', 
                                                    ifelse(names %like% 'Air Temperature Maximum', 'max_air_temperature_degF', 
                                                           ifelse(names %like% 'Air Temperature Minimum', 'min_air_temperature_degF', NA))))),
             partial_depth = stringr::str_extract(names_wo_id, "[[:digit:]]+"),
             new_names = ifelse(names == 'Date', 'date', paste0(partial_element, partial_depth,'in')),
             new_names = str_remove(new_names, 'NAin'))
    
    #define the final tibble
    final = temp %>% rename_at(vars(name_conversion_temp$names), function(x) name_conversion_temp$new_names) %>%
      as_tibble() %>%
      pivot_longer(cols = -c(date)) %>%
      mutate(site_id = valid_nrcs_sites$site_id[i]) %>%
      select(site_id, date, name, value) %>%
      drop_na()
    
    return(final)
  }, error = function(e){
    return(NA)
  }
  )
}
#stop cluster (memory management)
stopCluster(cl)

# this is what works, NRCS API doesnt do well with multiple requests. 
#alternative approach (single requests)
out = list()

for(i in 1:length(valid_nrcs_sites$site_id)){ # and try catch
  tryCatch({
    tictoc::tic()
    #GUI builder https://wcc.sc.egov.usda.gov/reportGenerator/
    baseURL = "https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customMultiTimeSeriesGroupByStationReport/daily/start_of_period/"
    endURL = "/POR_BEGIN,POR_END/SMV:-2:value,SMV:-4:value,SMV:-8:value,SMV:-20:value,SMV:-40:value,STV:-2:value,STV:-4:value,STV:-8:value,STV:-20:value,STV:-40:value,PRCP::value,TMAX::value,TMIN::value"
    fullURL = paste0(baseURL, valid_nrcs_sites$site_id[i] %>% parse_number(), ':', valid_nrcs_sites$state[i], ':', valid_nrcs_sites$ntwk[i], endURL)
    
    #get raw data and parse
    temp = RCurl::getURL(fullURL) %>%
      read_csv(., skip = 69)
    
    #compute name conversion
    name_conversion_temp = tibble(names = names(temp),
                                  names_wo_id = names(temp) %>% str_remove(., valid_nrcs_sites$site_id[i] %>% parse_number() %>% as.character())) %>%
      mutate(partial_element = ifelse(names %like% 'Soil Temperature', 'soil_temperature_', 
                                      ifelse(names %like% 'Moisture', 'soil_moisture_', 
                                             ifelse(names %like% 'Precipitation Increment', 'precipitation_increment_in', 
                                                    ifelse(names %like% 'Air Temperature Maximum', 'max_air_temperature_degF', 
                                                           ifelse(names %like% 'Air Temperature Minimum', 'min_air_temperature_degF', NA))))),
             partial_depth = stringr::str_extract(names_wo_id, "[[:digit:]]+in"),
             partial_depth = stringr::str_extract(partial_depth, "[[:digit:]]+"),
             new_names = ifelse(names == 'Date', 'date', paste0(partial_element, partial_depth,'in')),
             new_names = str_remove(new_names, 'NAin'))
    
    #define the final tibble
    final = temp %>% rename_at(vars(name_conversion_temp$names), function(x) name_conversion_temp$new_names) %>%
      as_tibble() %>%
      pivot_longer(cols = -c(date)) %>%
      mutate(site_id = valid_nrcs_sites$site_id[i]) %>%
      select(site_id, date, name, value) %>%
      drop_na()
    
    out[[i]] = final
    final
    print(i/length(valid_nrcs_sites$site_id)*100)
    tictoc::toc()
    
  }, error = function(e){
    out[[i]] = NA
    print(paste0('error on ', i))
    tictoc::toc()
  })
}


# final data
#make it long (tidy)
data_long = out[!is.na(out)] %>% 
  bind_rows() %>%
  tidyr::drop_na(value)

write_csv(data_long, '/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-long.csv')

#write it without precip and temp
data_long_sm_only = data_long %>%
  dplyr::filter(name %notin% c('max_air_temperature_degF', 'min_air_temperature_degF', 'precipitation_increment_in'))

write_csv(data_long_sm_only, '/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-long-soil-data-only.csv')

#filter final site meta data for data that actually has soil moisture (~10 dropped or so)
valid_nrcs_sites_final = valid_nrcs_sites %>%
  filter(site_id %in% unique(data_long$site_id))

#filter and extract geospatial information
valid_nrcs_sites_final_table = valid_nrcs_sites_final %>%
  mutate(latitude = unlist(map(valid_nrcs_sites_final$geometry,2)),
         longitude = unlist(map(valid_nrcs_sites_final$geometry,1))) %>%
  select(ntwk, site_id, latitude, longitude, elev_ft) %>%
  as_tibble() %>%
  select(-geometry) %>%
  `colnames<-`(c('network', 'site_id', 'latitude', 'longitude', 'elevation_ft'))

write_csv(valid_nrcs_sites_final_table, '/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-meta.csv')

valid_nrcs_sites_final_table_soil = valid_nrcs_sites_final_table %>%
  filter(site_id %in% unique(data_long_sm_only$site_id))

write_csv(valid_nrcs_sites_final_table_soil, '/home/zhoylman/soil-moisture-validation-data/processed/nrcs-soil-moisture/nrcs-soil-moisture-metasoil-data-only.csv')

