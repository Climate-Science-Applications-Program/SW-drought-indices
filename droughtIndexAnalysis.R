# Download data/process drought indices, analyze results
# MAC 06/02/24

# load libraries
library(SPEI)
library(dplyr)
library(ggplot2)

# load functions
source('getPRISMpoint.R')
source('calcZscore.R')


# GET DATA
# Tucson 32.12933603670422, -110.94978964298414
# Prescott 34.54491250593747, -112.46996711781452
lat<-34.54
lon<-(-112.45)
data<-getPRISMpoint(lat,lon,"1895-01-01","2023-12-31")
# extract results
meta<-data[[2]]
data<-data[[1]]

##### get sun hour estimate for P-M calculation ----
# get monthly average sunshine hours
dates<-as.data.frame(seq.Date(as.Date("2000-01-01"),as.Date("2000-12-31"), by="day")); colnames(dates)<-"date"
dates$month<-as.numeric(format(dates$date,"%m"))
dates$doy<-as.numeric(format(dates$date,"%j"))
dates$DayL<-solrad::DayLength(dates$doy,lat)
# get monthly averages
moSun<-dates %>% group_by(month) %>%
  summarise(avgHrs=mean(DayL))

# calculate sunshine hours
yrs<-seq(min(data$year),max(data$year),1)

tsun<-list()
for(i in 1:length(yrs)){
  dataYr<-subset(data, year==yrs[i])
  temp<-macroBiome::cliBrtSunDurFrcPoints(
    dataYr$avgT,
    dataYr$precip,
    lat,
    lon,
    meta$elev/3.281,
    year = yrs[i],
    aprchSIM = c("Solar123")
  )
  temp<-cbind.data.frame(moSun, t(temp))
  colnames(temp)[3]<-"sunFrac"
  temp$sunhours<-temp$avgHrs*temp$sunFrac
  temp$year<-yrs[i]
  tsun[[i]]<-temp
  print(yrs[i])
}

tsun<-do.call(rbind.data.frame, tsun)
save(tsun, file = "tsunSample_Prescott.RData")
load("C:/Users/Crimmins/OneDrive - University of Arizona/RProjects/SWDroughtIndices/tsunSample_Prescott.RData")
#####

#####
load("C:/Users/Crimmins/OneDrive - University of Arizona/RProjects/SWDroughtIndices/tsunSample_Prescott.RData")
# calculate PET water balances
# thorntwaite
tho <- thornthwaite(data$avgT, lat)
har <- hargreaves(data$minT, data$maxT, lat = lat)
harPrec <- hargreaves(data$minT, data$maxT, lat = lat, Pre=data$precip)
# penman, needs total sunshine hours
#penman<-penman(data$minT, data$maxT, tsun = tsun$sunhours, lat = lat, z = meta$elev/3.281, na.rm = TRUE)

# plot
#plot(ts(cbind(tho, har, harPrec,penman), fr = 12))
#pairs(cbind(tho, har, harPrec, penman))
pairs(cbind(tho, har))

# calc bias, rmse between methods
# quantile regression of monthly values
######

##### Plot climatologies
climo<-cbind.data.frame(data,tho,har,harPrec)
climo<-climo %>% group_by(month) %>%
                  summarise(count = n(),
                            minT=mean(minT, na.rm = TRUE),
                            maxT=mean(maxT, na.rm = TRUE),
                            avgT=mean(avgT, na.rm = TRUE),
                            precip=sum(precip, na.rm = TRUE)/n(),
                            tho=sum(tho, na.rm=TRUE)/n(),
                            har=sum(har, na.rm=TRUE)/n(),
                            harPrec=sum(harPrec, na.rm=TRUE)/n())
climo$monthAbb<-as.factor(month.abb[climo$month])
longClimo<-tidyr::pivot_longer(climo, cols = 3:9, names_to = "vars", values_to = "values")
longClimo$varType<-ifelse(longClimo$vars %in% c("minT","maxT","avgT"),"temps","waterBal")

# climograph
ggplot(longClimo, aes(month,values, color=vars))+
  geom_line()+
  facet_wrap(.~varType, ncol = 1, scales = "free_y")+
  scale_x_continuous(breaks = seq(1, 12, by = 1))+
  theme_bw()+
  ylab("deg C and mm")+
  ggtitle("Monthly Average Temperature and Water Balance Vars")

# time series plots - all months
ggplot(data, aes(year, maxT-minT, color=as.factor(month)))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()
# annual time series plots
annClimo<-cbind.data.frame(data,tho,har,harPrec)
annClimo<-annClimo %>% group_by(year) %>%
  summarise(count = n(),
            minT=mean(minT, na.rm = TRUE),
            maxT=mean(maxT, na.rm = TRUE),
            avgT=mean(avgT, na.rm = TRUE),
            #precip=sum(precip, na.rm = TRUE),
            tho=sum(precip-tho, na.rm=TRUE),
            har=sum(precip-har, na.rm=TRUE),
            harPrec=sum(precip-harPrec, na.rm=TRUE))

longAnnClimo<-tidyr::pivot_longer(annClimo, cols = 3:8, names_to = "vars", values_to = "values")
longAnnClimo$varType<-ifelse(longAnnClimo$vars %in% c("minT","maxT","avgT"),"temps","waterBal")
# plot
ggplot(longAnnClimo, aes(year,values, color=vars))+
  geom_line()+
  facet_wrap(.~varType, ncol = 1, scales = "free_y")+
  theme_bw()+
  geom_smooth(method = "lm")+
  ylab("deg C and mm")+
  ggtitle("Annual Average Temperature and Water Balance Vars")+
  theme_bw()  

#####


#####
# calculate drought indices
scales<-c(3,6,12)
refStart<-c(1991,1)
refEnd<-c(2020,12)
##### SPI
# calculate all full SPI timescales
fullSPI<-lapply(1:length(scales), function(x) spi(data$precip,scales[x])$fitted)
# clean Inf values
fullSPI<-lapply(fullSPI, function(x) replace(x, is.infinite(x),NA))
  fullSPI<-do.call(cbind, fullSPI)
  colnames(fullSPI)<-paste0("spiFull",scales)
fullSPI<-as.data.frame(fullSPI)
# calculate all shorter ref period timescales
shortSPI<-lapply(1:length(scales), function(x)
                 spi(ts(data$precip, freq=12, start=c(1895,1)),
                     scales[x], ref.start = refStart, ref.end = refEnd)$fitted)
# clean Inf values
shortSPI<-lapply(shortSPI, function(x) replace(x, is.infinite(x),NA))
  shortSPI<-do.call(cbind,shortSPI)
  colnames(shortSPI)<-paste0("spiShort",scales)
shortSPI<-as.data.frame(shortSPI)
# look at results
  pairs(cbind(shortSPI$spiShort3, fullSPI$spiFull3))

# add in diagnostics for SPI values --- median vs mean...
#####  
  
##### SPEI -thornthwaite ----
# calculate all full SPEI timescales
fullSPEI_th<-lapply(1:length(scales), function(x) spei(data$precip-tho,scales[x])$fitted)
# clean Inf values
fullSPEI_th<-lapply(fullSPEI_th, function(x) replace(x, is.infinite(x),NA))
  fullSPEI_th<-do.call(cbind, fullSPEI_th)
  colnames(fullSPEI_th)<-paste0("speiFull_th",scales)
fullSPEI_th<-as.data.frame(fullSPEI_th)
# calculate all shorter ref period timescales
shortSPEI_th<-lapply(1:length(scales), function(x)
                 spi(ts(data$precip-tho, freq=12, start=c(1895,1)),
                     scales[x], ref.start = refStart, ref.end = refEnd)$fitted)
# clean Inf values
shortSPEI_th<-lapply(shortSPEI_th, function(x) replace(x, is.infinite(x),NA))
  shortSPEI_th<-do.call(cbind,shortSPEI_th)
  colnames(shortSPEI_th)<-paste0("speiShort_th",scales)
shortSPEI_th<-as.data.frame(shortSPEI_th)
# look at results
pairs(cbind(shortSPEI_th$speiShort_th3, fullSPEI_th$speiFull_th3))
#####

##### SPEI -thornthwaite ----
# calculate all full SPEI timescales
fullSPEI_th<-lapply(1:length(scales), function(x) spei(data$precip-tho,scales[x])$fitted)
# clean Inf values
fullSPEI_th<-lapply(fullSPEI_th, function(x) replace(x, is.infinite(x),NA))
  fullSPEI_th<-do.call(cbind, fullSPEI_th)
  colnames(fullSPEI_th)<-paste0("speiFull_th",scales)
fullSPEI_th<-as.data.frame(fullSPEI_th)
# calculate all shorter ref period timescales
shortSPEI_th<-lapply(1:length(scales), function(x)
                 spei(ts(data$precip-tho, freq=12, start=c(1895,1)),
                     scales[x], ref.start = refStart, ref.end = refEnd)$fitted)
# clean Inf values
shortSPEI_th<-lapply(shortSPEI_th, function(x) replace(x, is.infinite(x),NA))
  shortSPEI_th<-do.call(cbind,shortSPEI_th)
  colnames(shortSPEI_th)<-paste0("speiShort_th",scales)
shortSPEI_th<-as.data.frame(shortSPEI_th)
# look at results
pairs(cbind(shortSPEI_th$speiShort_th3, fullSPEI_th$speiFull_th3))
#####


##### SPEI - hargreaves ----
# calculate all full SPEI timescales
fullSPEI_har<-lapply(1:length(scales), function(x) spei(data$precip-har,scales[x])$fitted)
# clean Inf values
fullSPEI_har<-lapply(fullSPEI_har, function(x) replace(x, is.infinite(x),NA))
  fullSPEI_har<-do.call(cbind, fullSPEI_har)
  colnames(fullSPEI_har)<-paste0("speiFull_har",scales)
fullSPEI_har<-as.data.frame(fullSPEI_har)
# calculate all shorter ref period timescales
shortSPEI_har<-lapply(1:length(scales), function(x)
                  spei(ts(data$precip-har, freq=12, start=c(1895,1)),
                  scales[x], ref.start = refStart, ref.end = refEnd)$fitted)
# clean Inf values
shortSPEI_har<-lapply(shortSPEI_har, function(x) replace(x, is.infinite(x),NA))
  shortSPEI_har<-do.call(cbind,shortSPEI_har)
  colnames(shortSPEI_har)<-paste0("speiShort_har",scales)
shortSPEI_har<-as.data.frame(shortSPEI_har)
# look at results
pairs(cbind(shortSPEI_har$speiShort_har3, fullSPEI_har$speiFull_har3))
#####

##### standardized variables
# avgTemp
tavgZ<-climZ_mean(scales,data$avgT,data$date,data$month)
colnames(tavgZ)[3:(length(scales)+2)]<-paste0("tavg",colnames(tavgZ)[3:(length(scales)+2)])
# minTemp
tminZ<-climZ_mean(scales,data$minT,data$date,data$month)
colnames(tminZ)[3:(length(scales)+2)]<-paste0("tmin",colnames(tminZ)[3:(length(scales)+2)])
# maxTemp
tmaxZ<-climZ_mean(scales,data$maxT,data$date,data$month)
colnames(tmaxZ)[3:(length(scales)+2)]<-paste0("tmax",colnames(tmaxZ)[3:(length(scales)+2)])
# thoPET
thoZ<-climZ_sum(scales,tho,data$date,data$month)
colnames(thoZ)[3:(length(scales)+2)]<-paste0("tho",colnames(thoZ)[3:(length(scales)+2)])
# harPET
harZ<-climZ_sum(scales,har,data$date,data$month)
colnames(harZ)[3:(length(scales)+2)]<-paste0("har",colnames(harZ)[3:(length(scales)+2)])
## water balance Z's - equivalent to SPEI without fitting 
  # # thoPET-P
  # thoZwb<-climZ_sum(scales,data$precip-tho,data$date,data$month)
  # colnames(thoZwb)[3:(length(scales)+2)]<-paste0("thoWB",colnames(thoZwb)[3:(length(scales)+2)])
  # # harPET-P
  # harZwb<-climZ_sum(scales,data$precip-har,data$date,data$month)
  # colnames(harZwb)[3:(length(scales)+2)]<-paste0("harWB",colnames(harZwb)[3:(length(scales)+2)])

#####
# analysis plots
#pairs(cbind.data.frame(fullSPI$spiFull3,fullSPEI_har$speiFull_har3,fullSPEI_th$speiFull_th3))
#GGally::ggpairs(cbind.data.frame(fullSPI$spiFull3,fullSPEI_har$speiFull_har3,fullSPEI_th$speiFull_th3))
#GGally::ggpairs(cbind.data.frame(fullSPI$spiFull6,fullSPEI_har$speiFull_har6,fullSPEI_th$speiFull_th6))
#GGally::ggpairs(cbind.data.frame(fullSPI$spiFull12,fullSPEI_har$speiFull_har12,fullSPEI_th$speiFull_th12))
#temp<-cbind.data.frame(data$month,fullSPI$spiFull3,fullSPEI_har$speiFull_har3,fullSPEI_th$speiFull_th3)
#GGally::ggpairs(temp, columns=2:4, aes(color=as.factor(data$month), alpha = 0.5))

# drought index time series plots
# 3-month scale
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull3,fullSPEI_har$speiFull_har3,fullSPEI_th$speiFull_th3, tempZ$tempZ_3)
colnames(temp)<-c("date","month","year","spi3","spei3_har","spei3_th","tavgZ")
temp<-tidyr::pivot_longer(temp, cols = 4:7, names_to = "vars", values_to = "values")

ggplot(temp, aes(date, values, color=vars))+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_x_date(date_breaks = "6 months", 
               labels=scales::date_format("%b-%Y"),
               limits = as.Date(c('2001-01-01','2003-01-01')))+
  theme_bw()

# 12-month scale
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull12,fullSPEI_har$speiFull_har12,fullSPEI_th$speiFull_th12, tempZ$tempZ_12)
colnames(temp)<-c("date","month","year","spi12","spei12_har","spei12_th","tavgZ")
temp<-tidyr::pivot_longer(temp, cols = 4:7, names_to = "vars", values_to = "values")

ggplot(temp, aes(date, values, color=vars))+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_x_date(date_breaks = "6 months", 
               labels=scales::date_format("%b-%Y"),
               limits = as.Date(c('2000-01-01','2022-01-01')))+
  theme_bw()

# 12-month scale --- SPI-SPEI difference plots
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull12,fullSPEI_har$speiFull_har12,fullSPEI_th$speiFull_th12, tempZ$tempZ_12)
colnames(temp)<-c("date","month","year","spi12","spei12_har","spei12_th","tavgZ")
  temp$diffHar<-temp$spei12_har-temp$spi12
  temp$diffTho<-temp$spei12_th-temp$spi12
  temp<-temp[,-c(4:6)]
temp<-tidyr::pivot_longer(temp, cols = 4:6, names_to = "vars", values_to = "values")

ggplot(temp, aes(date, values, color=vars))+
  geom_line()+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = -1, linetype="dashed")+
  scale_x_date(date_breaks = "6 months", 
               labels=scales::date_format("%b-%Y"),
               limits = as.Date(c('2015-01-01','2017-01-01')))+
  theme_bw()

# 12-month scale --- climatological period of record comparison
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull12,fullSPEI_har$speiFull_har12,fullSPEI_th$speiFull_th12,
                       shortSPI$spiShort12, shortSPEI_har$speiShort_har12, shortSPEI_th$speiShort_th12, tempZ$tempZ_12)
colnames(temp)<-c("date","month","year","spi12_full","spei12_har_full","spei12_th_full",
                  "spi12_short","spei12_har_short","spei12_th_short","tavgZ")
temp<-tidyr::pivot_longer(temp, cols = 4:10, names_to = "vars", values_to = "values")

ggplot(subset(temp, vars %in% c("spi12_full","spi12_short")), aes(date, values, color=vars))+
  geom_line()+
  geom_hline(yintercept = 0)+
  scale_x_date(date_breaks = "6 months", 
               labels=scales::date_format("%b-%Y"),
               limits = as.Date(c('2019-01-01','2022-01-01')))+
  theme_bw()


#####
# drought index comparison scatterplots -- 12month
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull12,fullSPEI_har$speiFull_har12,fullSPEI_th$speiFull_th12,
                       shortSPI$spiShort12, shortSPEI_har$speiShort_har12, shortSPEI_th$speiShort_th12, tempZ$tempZ_12)
colnames(temp)<-c("date","month","year","spi12_full","spei12_har_full","spei12_th_full",
                  "spi12_short","spei12_har_short","spei12_th_short","tavgZ")

ggplot(subset(temp, year>=1981), aes(spi12_full,tavgZ ,color=spei12_har_full))+
  geom_point()+
  ylim(-4,4)+
  xlim(-4,4)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_abline(slope=1,intercept = 0)+
  geom_smooth(method = "lm")+
  scale_color_distiller(type="div", palette="BrBG", direction=1)+
  theme_bw()
  
ggplot(temp, aes(spi12_full,spei12_har_full,color=tavgZ))+
  geom_point(alpha=1)+
  ylim(-4,0)+
  xlim(-4,0)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method = "lm")+
  scale_color_distiller(type="div", palette="RdBu", direction=-1)+
  theme_bw()

# drought index comparison scatterplots -- 3 month
temp<-cbind.data.frame(data$date,data$month,data$year,fullSPI$spiFull3,fullSPEI_har$speiFull_har3,fullSPEI_th$speiFull_th3,
                       shortSPI$spiShort3, shortSPEI_har$speiShort_har3, shortSPEI_th$speiShort_th3, tempZ$tempZ_3)
colnames(temp)<-c("date","month","year","spi3_full","spei3_har_full","spei3_th_full",
                  "spi3_short","spei3_har_short","spei3_th_short","tavgZ")

ggplot(subset(temp, month==7), aes(spi3_full,tavgZ ,color=spei3_har_full))+
  geom_point()+
  ylim(-4,4)+
  xlim(-4,4)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  #geom_abline(slope=1,intercept = 0)+
  geom_smooth(method = "lm")+
  scale_color_distiller(type="div", palette="BrBG", direction=1)+
  theme_bw()

ggplot(temp, aes(spi3_full,spei3_th_full,color=tavgZ))+
  geom_point(alpha=1)+
  ylim(-4,0)+
  xlim(-4,0)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline(slope=1,intercept = 0)+
  geom_smooth(method = "lm")+
  scale_color_distiller(type="div", palette="RdBu", direction=-1)+
  theme_bw()

