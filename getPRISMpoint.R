# function to pull point time series of PRISM data from RCC-ACIS webservice
# MAC 05/29/24

# USAGE
# getPRISMpoint(lat, lon, startDate, endDate)
#
# Arguments
# lat...point latitude in decimal degrees
# lon...point longitude in decimal degrees
# startDate...first date for query, not before 1895-01-01 for PRISM data; YYYY-MM-DD format
# endDate...last date for query, not after current month-1; ; YYYY-MM-DD format
#
# Example
# df<-getPRISMpoint(34,-110,"1895-01-01","2023-12-31")

# import functions
import::from(RCurl, postForm)
import::from(jsonlite, fromJSON)
import::from(magrittr, "%>%")
import::from(dplyr, mutate_if)

# define function
getPRISMpoint<-function(lat,lon,startDate,endDate){
  
  ##### download and process data for location ---- 
  # set location 
  #lat=32.071 
  #lon=-110.640 
  # set dates
  #startDate<-"1895-01-01"
  #endDate<-"2023-12-31"
  
  #download daily PRISM 
  jsonQuery=paste0('{"loc":"',lon,',',lat,'","sdate":"',startDate,'","edate":"',endDate,'","grid":"21",
                                    "meta":"ll,elev","elems":[{"name":"mly_pcpn","units":"mm"},{"name":"mly_mint","units":"degreeC"},{"name":"mly_maxt","units":"degreeC"}]}')
  
  outData<-postForm("http://data.rcc-acis.org/GridData",.opts = list(postfields = jsonQuery, 
                                                                     httpheader = c('Content-Type' = 'application/json', Accept = 'application/json')))
  # json to datframe
  outData<-fromJSON(outData)
  data<-data.frame(outData$data)
  #
  
  # get metadata
  # ll<-data.frame(matrix(unlist(outData$meta), nrow=length(outData$meta), byrow=T))
  meta<-outData$meta
  meta<-do.call(cbind.data.frame, meta)

  # character to numeric conversion
  data[,2:4]<- data[,2:4] %>% mutate_if(is.character, as.numeric)
  # set colnames
  colnames(data)<-c("date","precip","minT","maxT")
  # calculate daily average temperature
  data$avgT<-(data$maxT+data$minT)/2
  
  # convert date column
  data$date<-as.Date(paste0(data$date,"-01"), format='%Y-%m-%d')
  
  # add in date elements
  data$month<-as.numeric(format(data$date, "%m"))
  data$year<-as.numeric(format(data$date, "%Y"))
  ######
  
  # create list to return
  data<-list(data,meta)
  
  # return data
  return(data)
  
}


