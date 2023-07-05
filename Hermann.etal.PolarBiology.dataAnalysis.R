

rm(list=ls())
# Top ---------------------------------------------------------------------

set.seed(42)


library(tidyverse)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(adehabitatLT)
library(ggrepel)
library(gganimate)
library(move)
library(moveVis)
library(viridis)
library(RColorBrewer)
library(lubridate)
library(argosfilter)
library(pvclust)
library(TraMineR)
library(lubridate)
library(ggmosaic)
library(glmm)
library(modelr)
library(geosphere)
library(adehabitatHR)
library(ggmap)
library(ggsn)
library(ggnewscale)
library(ggstar)
library(ggspatial)
library(maptools)
library(lme4)
library(moments)
library(magick)
library(remotes)
library(readr)
library(sp)
library(sf)
library(ggpubr)
#remotes::install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = F) #If needing RSP install
library(RSP)
library(actel)
library(latticeDensity)


#To use, add command for scale_fill_manual(values=cbPalette)
theme_set(theme_bw(base_size=25))


'%notin%'<-Negate('%in%')

tags<-read.csv("2019_Sculpin_Tagged.csv")%>%
  filter(Tag_Type=="V9PA")%>%
  mutate(Fish_detailID=ifelse(Species=="Slimy",paste0("Msco",substr(Fish_ID,3,5)),paste0("Mqua",substr(Fish_ID,3,5))))

surgeryTime<-read.csv("2019_Sculpin_Tagged.csv")%>%
  filter(Tag_Type=="V9PA")%>%
  dplyr::select(Date,Anesthesia=Anestheisa,Surgery.start=Tag_Start,Surgery.stop=Tag_Stop,Recovery,Release)%>%
  mutate(Anesthesia=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Anesthesia))),
         Surgery.start=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.start))),
         Surgery.stop=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Surgery.stop))),
         Recovery=mdy_hm(paste(Date,Recovery)),
         Release=mdy_hm(paste(Date,gsub("([a,p])"," \\1",Release))))%>%
  filter(!is.na(Anesthesia))%>%dplyr::select(-Date)%>%
  mutate(Surgery.start=if_else(as.numeric(difftime(Surgery.start,Anesthesia,units="mins"))<0,
                               Surgery.start+days(1),Surgery.start),
         Surgery.stop=if_else(as.numeric(difftime(Surgery.stop,Surgery.start,units="mins"))<0,
                              Surgery.stop+days(1),Surgery.stop),
         Recovery=if_else(as.numeric(difftime(Recovery,Surgery.stop,units="mins"))<0,
                          Recovery+days(1),Recovery),
         Release=if_else(as.numeric(difftime(Release,Recovery,units="mins"))<0,
                         Release+days(1),Release))
#Time in anesthesia before surgery starts
mean(as.numeric(difftime(surgeryTime$Surgery.start,surgeryTime$Anesthesia,units="mins")),na.rm=T)
sd(as.numeric(difftime(surgeryTime$Surgery.start,surgeryTime$Anesthesia,units="mins")),na.rm=T)/sqrt(nrow(filter(surgeryTime,!is.na(Anesthesia))))
#Time under surgery
mean(as.numeric(difftime(surgeryTime$Surgery.stop,surgeryTime$Surgery.start,units="mins")),na.rm=T)
sd(as.numeric(difftime(surgeryTime$Surgery.stop,surgeryTime$Surgery.start,units="mins")),na.rm=T)/sqrt(nrow(filter(surgeryTime,!is.na(Surgery.stop))))
#Time for recovery
mean(as.numeric(difftime(surgeryTime$Release,surgeryTime$Recovery,units="mins")),na.rm=T)
range(as.numeric(difftime(surgeryTime$Release,surgeryTime$Recovery,units="mins")),na.rm=T)




# Acceleration and Pressure Exploration -----------------------------------


#Trying to match up the As with the Ps
V9APs<-read_csv("allsculpin.AP.data.csv",col_types = cols(Date.and.Time = col_datetime(format = "%Y-%m-%d %H:%M:%S"), 
                                                          Fish_ID = col_character(), Receiver_ID = col_character(), 
                                                          Release_date = col_date(format = "%Y-%m-%d"), 
                                                          Sensor.Measure = col_double(), Sensor.Value = col_double(), 
                                                          Tag_ID = col_character(), Year = col_character()))%>%
  group_by(Fish_ID)%>%
  mutate(time.delay_s=difftime(Date.and.Time,lag(Date.and.Time),units="secs"),
         pSensor.p=ifelse(Tag_Type=="V9A"&lag(Tag_Type)=="V9P"&time.delay_s<180,lag(Sensor.Measure),NA),
         pSensor.1=ifelse(Tag_Type=="V9A",ifelse(lag(Tag_Type)=="V9P"&time.delay_s<180,lag(Sensor.Measure),NA),NA),
         pSensor.2=ifelse(Tag_Type=="V9A",ifelse(lead(Tag_Type)=="V9P"&lead(time.delay_s)<180,lead(Sensor.Measure),NA),NA),
         pSensor.1=ifelse(pSensor.1<0,0,pSensor.1),
         pSensor.2=ifelse(pSensor.2<0,0,pSensor.2),
         season=ifelse(Date.and.Time>"2019-11-09","Winter","Summer"),
         Fish_ID=gsub("\\.5","",Fish_ID))

V9APs$pSensor<-ifelse(is.na(V9APs$pSensor.1),V9APs$pSensor.2,
                      ifelse(is.na(V9APs$pSensor.2),V9APs$pSensor.1,(V9APs$pSensor.1+V9APs$pSensor.2)/2))
#Only use this if we're doing the "gold" standard
V9APs$pSensor<-ifelse(is.na(V9APs$pSensor.p),NA,V9APs$pSensor.p)
##
V9APs$pSensor<-ifelse(V9APs$pSensor<0,0,V9APs$pSensor)

#How many does that cover??
nrow(filter(V9APs,!is.na(pSensor)))/nrow(filter(V9APs,Tag_Type=="V9A"))
#57.5% are the perfect pressure, they came from the pressure immediately at the start of the acceleration sampling window (180 lag and 0 lead)
#74.2% for either of the two P preceding the A (either the exact time or 180 sec away; 360 lag and 0 lead)
#78.1% within 180 seconds (the previous and/or next transmission was detected with a P value; 180 lag and lead)
#86.7% for any of the nearest P to that when A starts (at most 360 sec away, also options 0 and 180 sec away; 360 sec lag and 180 lead)
#94.4% within 12 minutes of the acceleration transmission (4 transmission times; 720 sec lag and lead)


V9APs%>%group_by(Fish_ID)%>%summarise(N=n())
V9APs%>%group_by(season)%>%summarise(N=n(),p=N/nrow(V9APs))


#Keep only the winter detections
V9APs<-filter(V9APs,season=="Winter")
V9APs%>%group_by(Fish_ID)%>%summarise(N=n())
V9APs%>%group_by(Station.name)%>%summarise(N=n())
(20891+13446)/nrow(V9APs)
V9APs%>%group_by(Fish_ID)%>%summarise(N=n_distinct(Station.name))%>%ungroup()%>%summarise(m=mean(N),s=sd(N)/sqrt(length(N)),r=range(N))
V9APs%>%group_by(Array)%>%summarise(N=n())
V9APs%>%group_by(Fish_ID)%>%summarise(N=n_distinct(Array))%>%ungroup()%>%summarise(m=mean(N),s=sd(N)/sqrt(length(N)),r=range(N))

V9APs%>%group_by(Tag_Type)%>%summarise(N=n(),p=n()/nrow(V9APs))

#How far off are those from the estimates?
summary(c(abs(V9APs$pSensor-V9APs$pSensor.1),abs(V9APs$pSensor-V9APs$pSensor.2)))
#V9APs$off1<-abs(V9APs$pSensor-V9APs$pSensor.1)
#V9APs$off2<-abs(V9APs$pSensor-V9APs$pSensor.2)

#Related to the distance off the bottom
V9APs<-V9APs%>%
  mutate(depth=ifelse(Tag_Type=="V9P",ifelse(Sensor.Measure<0,0,Sensor.Measure),pSensor), #positive depths
         receiverProp=ifelse(Tag_Type=="V9P",ifelse(Sensor.Measure<0,0,Sensor.Measure)/-Depth_m,NA), 
         minOff_Bottom=ifelse((depth*-1)-minDepth<=0,0,(depth*-1)-minDepth)) #Positive numbers are off bottom, if negative then assumed zero




#How fast are they changing from one to the next around the acceleration (up to 12 minutes before and 12 minutes after)
V9APs<-V9APs%>%
  ungroup()%>%
  mutate(diffpSensor=pSensor.2-pSensor.1,
         directionDepth=ifelse(diffpSensor<0,"Down",ifelse(diffpSensor>0,"Up","Stable")),
         rateDepth_ms=diffpSensor/(as.numeric(time.delay_s)+as.numeric(lead(time.delay_s))),
         accelCat=ifelse(Tag_Type=="V9P",NA,ifelse(Sensor.Measure>=3.4,"Burst",
                                                   ifelse(Sensor.Measure<=0.06,"Stationary",
                                                          ifelse(Sensor.Measure>0.06&Sensor.Measure<=0.78,"Low Activity","Swimming")))),
         accelCat=factor(accelCat,levels=c("Stationary","Low Activity","Swimming","Burst")))%>%
  left_join(dplyr::select(tags%>%mutate(Tag_ID=as.character(Tag_ID)),Tag_ID,Fish_detailID),by=c("Fish_ID"="Tag_ID"))

#How many of the As can I get a depth rate for?
nrow(filter(V9APs,!is.na(rateDepth_ms)))/nrow(filter(V9APs,Tag_Type=="V9A"))
#36.7% of accelerations have pressure measures at the next transmission both before AND after
#65.2% of accelerations have pressure measures within 12 minutes both before AND after (so I can calculate rate of change)


accurateIDs<-read.csv("allSculpininWindsor (Samples Found 03-16-2022)+Gen Processing Order.csv")%>%
  mutate(Sample.ID=gsub("^F","FH",ID))
table(accurateIDs$Species_GeneticID)

accurateTags<-left_join(tags,accurateIDs,by=c("Fish_ID"="Sample.ID"))



# Mapping with depth usage -------------------------------------------
spatial<-read.csv("spatial.csv")


#Doing it by distances
tremblay<-st_read("tremblay_and_Islands.shp")
tremblay<-fortify(tremblay)

tremblay<-filter(tremblay,OBJECTID=="6")

#Overlaying the detection heatmap on bathymetry
spatialBuffer<-st_as_sf(filter(spatial,Station.name!="Camp"),coords=c("Longitude","Latitude"))
st_crs(spatialBuffer)<-CRS("+proj=longlat +datum=WGS84")  
spatialBuffer <- st_transform(spatialBuffer, CRS("+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84"))
spatialBuffer <- st_buffer(spatialBuffer, 130)
spatialBuffer <- as_Spatial(spatialBuffer)
spatialBuffer <- spTransform(spatialBuffer,CRS("+proj=longlat +datum=WGS84"))

plotBuffer<-fortify(spatialBuffer,region="Station.name")
plotBuffer<-left_join(plotBuffer,detections%>%group_by(Station.name)%>%summarise(N=n()),by=c("id"="Station.name"))
plotBuffer[is.na(plotBuffer)]<-0


bathy<-raster("gebco_2020_n72.59971618652344_s72.22480773925781_w-81.42654418945312_e-80.47622680664062.nc",varname="elevation")
bathyDF<-as.data.frame(bathy[[1]],xy=T)
colnames(bathyDF)<-c("Long","Lat","z")
bathyDF$z<-ifelse(bathyDF$z>1,NA,bathyDF$z)
bathyDF$absZ<-abs(bathyDF$z)
bathyDF<-filter(bathyDF,z>-270)

bathyPoly<-rasterToPolygons(bathy)
bathyPoly<-fortify(bathyPoly,region="Elevation.relative.to.sea.level")
bathyPoly$id<-as.numeric(bathyPoly$id)
bathyPolyT<-filter(bathyPoly,id<=0 & !(long>-80.95&lat<72.35) & 
                     !(long>-80.7 &lat<72.45) & !long>-80.64) #Keep only the Tremblay area (below 0 and not over in the next area)
bathyPolyW<-filter(bathyPoly,id<=0 & id>-270)
bathyPolyL<-filter(bathyPoly,id>0)


load("googleTR.rda")


water.col <- colorRampPalette(c("purple4", "royalblue4", "royalblue1", "skyblue3", "lightblue3", "white"))

plotV9APs<-V9APs%>%
  filter(Tag_Type=="V9P")%>%
  mutate(depthColor=case_when(Sensor.Measure<10~"a",
                              Sensor.Measure>=10 & Sensor.Measure<50~"b",
                              Sensor.Measure>=50~"c"))%>%
  dplyr::select(Station.name,Longitude,Latitude,depthColor)%>%
  #arrange(Station.name,depthColor)%>%
  group_by(Station.name,depthColor)%>%
  mutate(N=n())%>%
  slice(which(row_number() %% 10 == 1))
plotV9APs2<-distinct(plotV9APs)%>%arrange(desc(N))

stamen<-get_stamenmap(bbox=bbox(bathy),maptype="terrain")
ggmap(stamen)


#Complete detections Map
bathyMap<-image_graph(width=4156.25,height=3906.25,res=300)
ggmap(stamen)+
  #Bathy background
  geom_contour_filled(data=filter(bathyPolyT,lat<75.502),aes(long,lat,z=-1*id),
                      color="black",breaks=c(-1000,0,10,50,100,250,300),alpha=0.6)+
  scale_fill_manual(values=rev(water.col(5)),name="Depth (m)",
                    labels=c("<10 m","10 - 50 m","50 - 100 m","100 - 250 m",">250 m"))+
  geom_star(data=filter(spatial,Type=="Release"),aes(Longitude,Latitude),fill="firebrick2",size=11)+ #Campsite
  geom_point(aes(x=-81.10518, y=72.35267),size=8,shape=22,fill="seagreen3")+ #Fyke Net
  new_scale_fill()+
  #Array
  #geom_point(data=filter(spatial,Type!="Release"),aes(Longitude,Latitude),color="grey5",size=0.85)+
  geom_polygon(data=plotBuffer,aes(long,lat,group=group),fill=rgb(0,0,0,alpha=0.5))+
  #Detections+Depths
  geom_point(data=plotV9APs2,
             aes(Longitude,Latitude,fill=depthColor,size=N),
             shape=21)+
  scale_fill_manual(values=rev(water.col(5))[1:3],guide="none")+
  scale_size_continuous(range=c(3,11),guide="none")+
  #Gate labels
  geom_label(aes(x=-80.822,y=72.433,label="A"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-80.858,y=72.397,label="B"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-80.940,y=72.355,label="C"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.039,y=72.329,label="D"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.078,y=72.308,label="E"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.119,y=72.277,label="F"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  geom_label(aes(x=-81.188,y=72.246,label="G"),color="white",fill="black",size=18,vjust=0.25,hjust=0.5)+
  #"Legend" of shape types
  annotate("rect",xmin=-81.41,xmax=-81.215,ymin=72.313,ymax=72.359,fill="white",alpha=0.85,color="black",size=1)+
  geom_text(aes(x=-81.4,y=72.353,label="Location"),hjust=0,size=10)+
  geom_text(aes(x=-81.37,y=72.341,label="Campsite"),hjust=0,size=8)+
  geom_text(aes(x=-81.37,y=72.331,label="Fyke Net"),hjust=0,size=8)+
  geom_text(aes(x=-81.37,y=72.321,label="Receivers"),hjust=0,size=8)+
  geom_star(aes(x=-81.39,y=72.34),size=8,fill="firebrick2")+ #Camp
  geom_point(aes(x=-81.39, y=72.33),size=8,shape=22,fill="seagreen3")+ #Fyke Net
  geom_point(aes(x=-81.39, y=72.32),size=8,shape=21,fill="white")+ #All the receiver points
  geom_point(aes(x=-81.39, y=72.32),size=2)+ #All the receiver points
  #Other aesthetics
  theme(legend.position=c(0.015,0.63),legend.justification=0,legend.title=element_text(size=30),
        legend.background=element_rect(fill=rgb(1,1,1,alpha=0.88),color="black",size=1),
        legend.key.height=unit(26,"points"),legend.text=element_text(size=25),legend.margin=margin(10,20,20,15),
        axis.title=element_text(size=32),axis.text=element_text(size=25),
        axis.ticks=element_line(size=2),axis.ticks.length=unit(0.2,"cm"),
        plot.margin=margin(60,100,20,20), panel.border = element_rect(fill=NA,color="black",size=1.25))+
  scale_y_continuous(expand=expansion(0),breaks=c(72.25,72.3,72.35,72.4,72.45,72.5),
                     labels=c("72.25°N","72.30°N","72.35°N","72.40°N","72.45°N","72.50°N"),name="Latitude")+
  scale_x_continuous(expand=expansion(0),breaks=c(-81.4,-81.2,-81,-80.8,-80.6),
                     labels=c("-81.4°W","-81.2°W","-81.0°W","-80.8°W","-80.6°W"),name="Longitude")+
  scalebar(location="bottomleft",anchor=c(x=-81.25,y=72.491),dist=5,dist_unit="km",transform=T,height=0.025,
           x.min=-81.49,x.max=-80.61,y.min=72.236,y.max=72.501,st.dist=0.03,st.size=13,st.bottom=T,st.color="black")+
  north(x.min=-81.49,x.max=-80.61,y.min=72.24,y.max=72.501,
        anchor=c(x=-81.4,y=72.5),location="topleft",symbol=12)+
  annotate("text",x=-81.357,y=72.47,label="N",color="black",size=18,fontface="bold")+
  coord_cartesian(ylim=c(72.236,72.501),clip="on")
dev.off()


####Inset####
load("insetMap_googleMap.R")
countries<-map_data("world")%>%filter(region%in%c("Canada","USA"))
countriesLabels<-countries%>%
  group_by(region,subregion)%>%summarise(long=mean(long),lat=mean(lat))
canada<-readOGR("lpr_000b16a_e.shp")
canada <- spTransform(canada, CRSobj = CRS("+proj=longlat +ellps=WGS84"))

canadaDF<-fortify(canada)
canadaDF<-canadaDF%>%
  mutate(id=factor(id),
         idName=case_when(id=="0" ~ canada@data$PRENAME[1],
                          id=="1" ~ canada@data$PRENAME[2],
                          id=="2" ~ canada@data$PRENAME[3],
                          id=="3" ~ canada@data$PRENAME[4],
                          id=="4" ~ canada@data$PRENAME[5],
                          id=="5" ~ canada@data$PRENAME[6],
                          id=="6" ~ canada@data$PRENAME[7],
                          id=="7" ~ canada@data$PRENAME[8],
                          id=="8" ~ canada@data$PRENAME[9],
                          id=="9" ~ canada@data$PRENAME[10],
                          id=="10" ~ canada@data$PRENAME[11],
                          id=="11" ~ canada@data$PRENAME[12],
                          T ~ canada@data$PRENAME[13]))
canadaMap<-canadaDF%>%
  group_by(idName)%>%
  mutate(labLong=case_when(idName=="Ontario"~mean(long)+9,
                           idName=="Nunavut"~mean(long)-5,
                           TRUE~mean(long)),
         labLat=ifelse(idName=="Ontario",mean(lat)+15,mean(lat)-1))%>%
  dplyr::select(labLong,labLat,idName)%>%unique()%>%
  bind_rows(data.frame(labLong=-72.5,labLat=69,idName="Baffin Island"))

inset<-image_graph(width=1500,height=1500,res=300,bg="transparent")
ggmap(insetMap)+
  geom_polygon(data=canadaDF,aes(long,lat,group=group),fill="white",color="black",size=0.5,alpha=0.2)+
  geom_text(data=canadaMap%>%filter(idName%in%c("Nunavut","Manitoba","Ontario","Quebec")),
            aes(labLong,labLat,label=idName),color="white",size=6)+
  #scale_x_continuous(limits=c(-115,-45))+
  #scale_y_continuous(limits=c(40,74))+
  theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
        plot.background = element_blank(), panel.border=element_rect(fill=NA,color="black",size=2))+
  scalebar(location="bottomleft",dist=1000,dist_unit="km",transform=T,st.color="white",st.bottom=F,st.dist=0.05,
           x.min=-115,x.max=-45,y.min=40,y.max=74)+
  #geom_point(aes(-70.934167,43.135),size=4,color="red") #Added to get UNH in right place
  geom_rect(aes(xmin=-81.89,xmax=-79.6,ymin=71.7,ymax=72.7),fill="transparent",color="firebrick1",size=2)+ #Added to surround Tremblay
  coord_map(xlim=c(-95,-60),ylim=c(60,75))

dev.off()

fullMap<-image_composite(bathyMap,inset,offset="+2388+2054") #Fit right into the corner of the plot, r margin 20, plot width 1250, inset should be 450x450
#fullMap<-image_composite(bathyMap,inset,offset="+800+710") #Fit outside the map, make the r margin 100 and the plot width 1330, inset should be 525x525
fullMap




# #Histograms of AP data (in manu) --------------------------------------------------


#Acceleration histogram
oSkew<-filter(V9APs,Tag_Type=="V9A")%>%
  group_by(Fish_detailID)%>%
  summarise(s=skewness(Sensor.Measure))
plot<-filter(V9APs,Tag_Type=="V9A")%>%
  ggplot(aes(Sensor.Measure))+
  geom_histogram(boundary=0,binwidth=0.06,fill="grey80",color="black")+
  scale_x_continuous(name=expression(paste("ODBA (ms"^"-2"*")")),expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,3.5,by=0.5))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~Fish_detailID,nrow=2,scales="free_y")
oSkew<-filter(V9APs,Tag_Type=="V9A")%>%
  group_by(Fish_detailID)%>%
  summarise(s=skewness(Sensor.Measure))%>%
  left_join(ggplot_build(plot)$data[[1]]%>%
              group_by(PANEL)%>%
              summarise(max=max(count))%>%
              mutate(Fish_detailID=rev(unique(V9APs$Fish_detailID))))
tiff("figure2a_odbahist.tiff",width=4687.5,height=4687.5,res=300)
ggplot(filter(V9APs,Tag_Type=="V9A"),aes(Sensor.Measure))+
  geom_histogram(boundary=0,binwidth=0.06,fill="grey80",color="black")+
  geom_vline(xintercept=0.06,lty=2,size=1)+
  geom_text(data=oSkew,aes(x=2.5,y=max*0.9,label=paste("Skewness = ",round(s,3))),inherit.aes=F,size=8)+
  scale_x_continuous(name=expression(paste("ODBA (ms"^"-2"*")")),expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,3.5,by=0.5))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~Fish_detailID,nrow=3,scales="free_y")
dev.off()

#Depth histogram
dSkew<-filter(V9APs,Tag_Type=="V9P")%>%
  group_by(Fish_detailID)%>%
  summarise(s=skewness(depth))
plot<-filter(V9APs,Tag_Type=="V9P")%>%
  ggplot(aes(depth))+
  geom_histogram(boundary=0,binwidth=1,fill="grey80",color="black")+
  scale_x_continuous(name="Depth (m)",expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,55,by=10))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~Fish_detailID,nrow=2,scales="free_y")
dSkew<-filter(V9APs,Tag_Type=="V9P")%>%
  group_by(Fish_detailID)%>%
  summarise(s=skewness(depth))%>%
  left_join(ggplot_build(plot)$data[[1]]%>%
              group_by(PANEL)%>%
              summarise(max=max(count))%>%
              mutate(Fish_detailID=rev(unique(V9APs$Fish_detailID))))
tiff("figure2b_depthhist.tiff",width=4687.5,height=4687.5,res=300)
ggplot(filter(V9APs,Tag_Type=="V9P"),aes(depth))+
  geom_histogram(boundary=0,binwidth=1,fill="grey80",color="black")+
  geom_text(data=dSkew,aes(x=c(39,13,39,39,39),y=max*0.9,label=paste("Skewness = ",round(s,3))),inherit.aes=F,size=8)+
  scale_x_continuous(name="Depth (m)",expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,55,by=10))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~Fish_detailID,nrow=3,scales="free_y")
dev.off()

#Summaries of that data
summary(filter(V9APs,Tag_Type=="V9A")$Sensor.Measure)
sd(filter(V9APs,Tag_Type=="V9A")$Sensor.Measure)/sqrt(nrow(filter(V9APs,Tag_Type=="V9A")))
nrow(filter(V9APs,Tag_Type=="V9A" & Sensor.Measure>3.46))/nrow(filter(V9APs,Tag_Type=="V9A"))*100
skewness(filter(V9APs,Tag_Type=="V9A")$Sensor.Measure)

summary(filter(V9APs,Tag_Type=="V9P")$Sensor.Measure)
sd(filter(V9APs,Tag_Type=="V9P")$Sensor.Measure)/sqrt(nrow(filter(V9APs,Tag_Type=="V9P")))
skewness(filter(V9APs,Tag_Type=="V9P")$Sensor.Measure)

#By individual--odds are depth and evens are acceleration
V9APs%>%group_by(Tag_ID)%>%summarise(s=skewness(Sensor.Measure))
sd(c(2.42,2.92,2.58,2.91,2.78))/sqrt(5)
sd(c(4.04,3.47,2.28,-1.02,0.757))/sqrt(5)

summary(filter(V9APs,Tag_Type=="V9P")$receiverProp)
skewness(filter(V9APs,Tag_Type=="V9P")$receiverProp)



#Aggregate fish and facet by time
omonthlySkew<-filter(V9APs,Tag_Type=="V9A")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  group_by(month)%>%
  summarise(s=skewness(Sensor.Measure))
plot<-filter(V9APs,Tag_Type=="V9A")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  ggplot(aes(Sensor.Measure))+
  geom_histogram(boundary=0,binwidth=0.06,fill="grey80",color="black")+
  scale_x_continuous(name=expression(paste("ODBA (ms"^"-2"*")")),expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,3.5,by=0.5))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~month,nrow=2,scales="free_y")
omonthlySkew<-filter(V9APs,Tag_Type=="V9A")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  group_by(month)%>%
  summarise(s=skewness(Sensor.Measure))%>%
  left_join(ggplot_build(plot)$data[[1]]%>%
              group_by(PANEL)%>%
              summarise(max=max(count))%>%
              mutate(month=c("November","December","January","February","March","April")))%>%
  mutate(month=factor(month,levels=c("November","December","January","February","March","April")))
filter(V9APs,Tag_Type=="V9A")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  ggplot(aes(Sensor.Measure))+
  geom_histogram(boundary=0,binwidth=0.06,fill="grey80",color="black")+
  geom_text(data=omonthlySkew,aes(x=2.5,y=max*0.9,label=paste("Skewness = ",round(s,3))),inherit.aes=F,size=8)+
  geom_vline(xintercept=0.06,lty=2,size=1)+
  scale_x_continuous(name=expression(paste("ODBA (ms"^"-2"*")")),expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,3.5,by=0.5))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~month,nrow=2,scales="free_y")




#Skewness for this
dmonthlySkew<-filter(V9APs,Tag_Type=="V9P")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  group_by(month)%>%
  summarise(s=skewness(depth))
plot<-filter(V9APs,Tag_Type=="V9P")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  ggplot(aes(depth))+
  geom_histogram(boundary=0,binwidth=1,fill="grey80",color="black")+
  geom_text(data=dmonthlySkew,aes(x=35,y=1000,label=paste("Skewness = ",round(s,3))),inherit.aes=F)+
  scale_x_continuous(name="Depth (m)",expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,55,by=10))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~month,nrow=2,scales="free_y")
dmonthlySkew<-filter(V9APs,Tag_Type=="V9P")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  group_by(month)%>%
  summarise(s=skewness(depth))%>%
  left_join(ggplot_build(plot)$data[[1]]%>%
              group_by(PANEL)%>%
              summarise(max=max(count))%>%
              mutate(month=c("November","December","January","February","March","April")))%>%
  mutate(month=factor(month,levels=c("November","December","January","February","March","April")))
filter(V9APs,Tag_Type=="V9P")%>%
  mutate(month=factor(month(Date.and.Time,label=T,abbr=F),levels=c("November","December","January","February","March","April")))%>%
  ggplot(aes(depth))+
  geom_histogram(boundary=0,binwidth=1,fill="grey80",color="black")+
  geom_text(data=dmonthlySkew,aes(x=39,y=max*0.9,label=paste("Skewness = ",round(s,3))),inherit.aes=F,size=8)+
  scale_x_continuous(name="Depth (m)",expand=expansion(mult=c(0.05,0.05)),
                     breaks=seq(0,55,by=10))+
  scale_y_sqrt(name="Count",expand=expansion(mult=c(0,0.1)))+
  theme(panel.grid=element_blank(),panel.background=element_rect(color="black",size=1.5),panel.spacing=unit(0.5,"lines"),
        axis.ticks=element_line(size=2),axis.text=element_text(size=33),axis.title=element_text(size=50),strip.text=element_text(size=45))+
  facet_wrap(~month,nrow=2,scales="free_y")






# Sequences of AP data ----------------------------------------------------


#Inertia Method (in manu) ##### 

#A new method for assigning sequences that builds inertia
#If you have had multiple detections in a row, you can miss some
#But you can't miss more than you've gotten
#So you get a point for each detection you hit, and lose one for each detection you miss
inertiaSeqs<-V9APs%>%
  arrange(Fish_detailID,Date.and.Time)%>%
  mutate(inertia=0,
         seqNum=0)%>%
  dplyr::select(Fish_detailID,Date.and.Time,season,Tag_Type,Sensor.Measure,Mass,Station.name,receiverDepth=Depth_m,
                Latitude,Longitude,time.delay_s,inertia,seqNum)
n=0 
for (i in 2:nrow(inertiaSeqs)) {
  if (inertiaSeqs$Fish_detailID[i]!=inertiaSeqs$Fish_detailID[i-1]) { #If it's a new fish, start back down at zero
    inertiaSeqs$inertia[i]=0
    n=n+1
  } else if (inertiaSeqs$time.delay_s[i]<100) {
    inertiaSeqs$inertia[i]=inertiaSeqs$inertia[i-1]
  } else if (inertiaSeqs$time.delay_s[i]<=180) { #If the last 
    inertiaSeqs$inertia[i]=inertiaSeqs$inertia[i-1]+1
  } else if (inertiaSeqs$inertia[i-1]>0 & inertiaSeqs$time.delay_s[i]>180) {
    inertiaSeqs$inertia[i]=inertiaSeqs$inertia[i-1]-floor((as.numeric(floor(inertiaSeqs$time.delay_s[i]/180)))^1.1)
    if (inertiaSeqs$inertia[i]<=0) {
      inertiaSeqs$inertia[i]=0
      n=n+1
    }
  } else {
    inertiaSeqs$inertia[i]=0
    n=n+1
  }
  inertiaSeqs$seqNum[i]=n
}
inertiaSeqs<-inertiaSeqs%>%
  group_by(Fish_detailID)%>%
  mutate(seqNum=seqNum-min(seqNum)+1)

inertiaSeqSummary<-inertiaSeqs%>%
  group_by(Fish_detailID,seqNum)%>%
  summarise(N=n(),N_A=sum(Tag_Type=="V9A"),N_P=sum(Tag_Type=="V9P"),maxMissingTime=max(time.delay_s),
            duration_hrs=as.numeric(difftime(max(Date.and.Time),min(Date.and.Time),units="hours")))


#### Looking at those sequences longer than 8 hours ####
longInert<-filter(inertiaSeqSummary,duration_hrs>8)%>%dplyr::select(Fish_detailID,seqNum)

inertlongSeqs<-inertiaSeqs%>%
  right_join(longInert)%>%mutate(Sensor.Measure=ifelse(Tag_Type=="V9A",Sensor.Measure,Sensor.Measure*-1))%>%
  mutate(g = c(0, cumsum(diff(time.delay_s) > 180)),
         sg=paste0(Tag_Type,seqNum,g))

#Example plot for how inertia works
inertlongSeqs%>%filter(seqNum=="343")%>%
  ggplot()+
  geom_line(aes(x=Date.and.Time,inertia))
8^1.1






#### Plotting long sequences ####


seqPlot1<-ggplot(filter(inertlongSeqs, Fish_detailID=="Mqua034"))+
  geom_ribbon(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua034"),
              aes(Date.and.Time,ymin=Sensor.Measure-1.7,ymax=Sensor.Measure+1.7,group=sg),alpha=0.2)+
  geom_line(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua034"),
            aes(Date.and.Time,Sensor.Measure,color=Tag_Type,group=sg))+
  geom_line(data=filter(inertlongSeqs,Tag_Type=="V9A" & Fish_detailID=="Mqua034"),
            aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type,group=sg))+
  geom_point(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua034"),
             aes(Date.and.Time,Sensor.Measure,color=Tag_Type))+
  geom_point(data=filter(inertlongSeqs,Tag_Type=="V9A" & Fish_detailID=="Mqua034"),
             aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type))+
  scale_color_viridis_d(option="A",begin=0.2,end=0.7)+
  scale_y_continuous(name="Depth (m)",
                     sec.axis=sec_axis(trans=~./5,name=expression(paste("ODBA (ms"^"-2"*")"))))+
  scale_x_datetime(name="",date_labels = "%H:%M")+
  theme(axis.title.x=element_text(size=45),axis.text.x=element_text(size=15),
        axis.title.y.left = element_text(color="white",size=45),
        axis.title.y.right = element_text(color="white",size=45),
        axis.text.y.left = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.text.y.right = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        axis.ticks.y.left = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.ticks.y.right = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(45,0,0,0),
        legend.position="none",plot.background = element_rect(fill="transparent"))+
  facet_wrap(.~paste0("Seq: ",seqNum),scales="free",ncol=4) 

seqPlot2<-ggplot(filter(inertlongSeqs, Fish_detailID=="Mqua028"))+
  geom_ribbon(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua028"),
              aes(Date.and.Time,ymin=Sensor.Measure-1.7,ymax=Sensor.Measure+1.7,group=sg),alpha=0.2)+
  geom_line(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua028"),
            aes(Date.and.Time,Sensor.Measure,color=Tag_Type,group=sg))+
  geom_line(data=filter(inertlongSeqs,Tag_Type=="V9A" & Fish_detailID=="Mqua028"),
            aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type,group=sg))+
  geom_point(data=filter(inertlongSeqs,Tag_Type=="V9P" & Fish_detailID=="Mqua028"),
             aes(Date.and.Time,Sensor.Measure,color=Tag_Type))+
  geom_point(data=filter(inertlongSeqs,Tag_Type=="V9A" & Fish_detailID=="Mqua028"),
             aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type))+
  scale_color_viridis_d(option="A",begin=0.2,end=0.7)+
  scale_y_continuous(name="Depth (m)",
                     sec.axis=sec_axis(trans=~./5,name=expression(paste("ODBA (ms"^"-2"*")"))))+
  scale_x_datetime(name="",date_labels = "%H:%M")+
  theme(axis.title.x=element_text(size=45),axis.text.x=element_text(size=15),
        axis.title.y.left = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[2],size=45),
        axis.title.y.right = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[1],size=45),
        axis.text.y.left = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.text.y.right = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        axis.ticks.y.left = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.ticks.y.right = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position="none",plot.background = element_rect(fill="transparent"))+
  facet_wrap(.~paste0("Seq: ",seqNum),scales="free",ncol=4) 

seqPlot3<-ggplot(bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                      Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                           filter(inertlongSeqs, Fish_detailID=="Mqua026")))+
  geom_ribbon(data=bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                        Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                             filter(inertlongSeqs, Fish_detailID=="Mqua026" & Tag_Type=="V9P")),
              aes(Date.and.Time,ymin=Sensor.Measure-1.7,ymax=Sensor.Measure+1.7,group=sg),alpha=0.2)+
  geom_line(data=bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                      Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                           filter(inertlongSeqs, Fish_detailID=="Mqua026" & Tag_Type=="V9P")),
            aes(Date.and.Time,Sensor.Measure,color=Tag_Type,group=sg))+
  geom_line(data=bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                      Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                           filter(inertlongSeqs, Fish_detailID=="Mqua026" & Tag_Type=="V9A")),
            aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type,group=sg))+
  geom_point(data=bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                       Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                            filter(inertlongSeqs, Fish_detailID=="Mqua026" & Tag_Type=="V9P")),
             aes(Date.and.Time,Sensor.Measure,color=Tag_Type))+
  geom_point(data=bind_rows(data.frame(seqNum=c(98,99),Date.and.Time=c(ymd_hm("2020-10-10 12:00","2020-10-10 12:00")),
                                       Sensor.Measure=c(1,1),sg=c("dummy","dummy")),
                            filter(inertlongSeqs, Fish_detailID=="Mqua026" & Tag_Type=="V9A")),
             aes(Date.and.Time,Sensor.Measure*5,color=Tag_Type))+
  scale_color_viridis_d(option="A",begin=0.2,end=0.7)+
  scale_y_continuous(name="Depth (m)",
                     sec.axis=sec_axis(trans=~./5,name=expression(paste("ODBA (ms"^"-2"*")"))))+
  scale_x_datetime(name="Time",date_labels = "%H:%M")+
  theme(axis.title.x=element_text(size=45),axis.text.x=element_text(size=15),
        axis.title.y.left = element_text(color="white",size=45),
        axis.title.y.right = element_text(color="white",size=45),
        axis.text.y.left = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.text.y.right = element_text(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        axis.ticks.y.left = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[2]),
        axis.ticks.y.right = element_line(color=viridis(2,option="A",begin=0.2,end=0.7)[1]),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",plot.background = element_rect(fill="transparent"))+
  facet_wrap(.~paste0("Seq: ",seqNum),scales="free",ncol=4) 

tiff("figure3_sequences.tiff",width=7500,height=7500,res=300)
ggarrange(seqPlot1,seqPlot2,seqPlot3,nrow=3,heights=c(2,4,1),
          labels=c("Fish ID: Mqua034","Fish ID: Mqua028","Fish ID: Mqua026"),vjust=c(1.2,0,0),hjust=-0.25,font.label=list(size=45))
dev.off()

#### Frequency of bursts 
inertiaBursts<-inertiaSeqs%>%
  filter(Tag_Type=="V9A")%>%
  mutate(burst=ifelse(Sensor.Measure<=0.06,"rest","burst"),
         burstSeq=1)%>%
  dplyr::select(Date.and.Time,Fish_detailID,Station.name,Sensor.Measure,seqNum,burst,burstSeq)
n=1
for (i in 2:nrow(inertiaBursts)) {
  if ((inertiaBursts$burst[i]=="burst" & inertiaBursts$burst[i-1]=="rest") | 
      (inertiaBursts$burst[i]=="rest" & inertiaBursts$burst[i-1]=="burst") |
      inertiaBursts$seqNum[i]!=inertiaBursts$seqNum[i-1]) {
    n=n+1
  }
  inertiaBursts$burstSeq[i]<-n
}


#Calculate the duration of all the burst (and rest) periods
inertiaBursts<-inertiaBursts%>%
  group_by(Fish_detailID,seqNum,burstID=paste(burst,burstSeq,sep="-"))%>%
  mutate(startTime=min(Date.and.Time),
         endTime=max(Date.and.Time),
         duration_hrs=as.numeric(difftime(endTime,startTime,units="hours")),
         maxBurst=ifelse(burst=="burst",max(Sensor.Measure),NA),
         meanBurst=ifelse(burst=="burst",mean(Sensor.Measure),NA))

inertiaSeqBurst<-left_join(inertiaSeqs,inertiaBursts)%>%
  group_by(Fish_detailID,seqNum)%>%
  mutate(seqStartTime=min(Date.and.Time),
         seqEndTime=max(Date.and.Time),
         seqDuration=seqEndTime-seqStartTime,
         seqDuration_hrs=as.numeric(seqDuration)/60/60+(180/60/60))


inertiaSeqBurstSummary<-inertiaSeqBurst%>%
  filter(Tag_Type=="V9A")%>%
  group_by(Fish_detailID,seqNum,seqStartTime,seqDuration_hrs)%>%
  summarise(meanODBA=mean(Sensor.Measure))%>%distinct()
inertiaSeqDepthSummary<-inertiaSeqBurst%>%
  filter(Tag_Type=="V9P")%>%
  group_by(Fish_detailID,seqNum,seqStartTime,seqDuration_hrs)%>%
  summarise(meanDepth=mean(Sensor.Measure), 
            minDepth=min(Sensor.Measure),
            maxDepth=max(Sensor.Measure),
            rangeDepth=maxDepth-minDepth)%>%distinct()

inertiaSeqSummary<-full_join(inertiaSeqBurstSummary,inertiaSeqDepthSummary)%>%
  mutate(dDays=yday(seqStartTime)+
           (hour(seqStartTime)/24)+
           (minute(seqStartTime)/(24*60))+
           (second(seqStartTime)/(24*60*60)),
         dDays_diff=as.numeric(difftime(seqStartTime,ymd_hms("2019-11-09 00:00:00"),units="days")),
         atRest=ifelse(meanODBA>0.06,0,1),
         atSurface=ifelse(meanDepth>1.7,0,1))%>%
  arrange(Fish_detailID,seqNum)

nrow(filter(inertiaSeqSummary,rangeDepth<1 & seqDuration_hrs>0.05))/nrow(filter(inertiaSeqSummary,rangeDepth>=0 & seqDuration_hrs>0.05))





# Frequency and Regularity of Bursts --------------------------------------


timebetweenbursts<-filter(inertiaSeqBurst,burst=="burst")%>%
  ungroup()%>%
  dplyr::select(Fish_detailID,seqNum,burstID,startTime,endTime)%>%distinct()%>%
  group_by(Fish_detailID,seqNum)%>%
  summarise(time=difftime(startTime,lag(endTime),units="mins")) #Includes times within the same burst, so ignore negatives, and when there are no previous bursts in a sequences, so ignore NAs
meanTTB<-timebetweenbursts%>%
  group_by(Fish_detailID)%>%
  summarise(time=mean(time,na.rm=T))
mean(timebetweenbursts$time,na.rm=T)
sd(timebetweenbursts$time,na.rm=T)/sqrt(sum(timebetweenbursts$time>0,na.rm=T))

#Frequency of bursts
burstNum<-filter(inertiaSeqBurst,burst=="burst")%>%
  ungroup()%>%
  dplyr::select(Fish_detailID,burstID,duration_hrs,Date.and.Time)%>%distinct()%>%
  group_by(Fish_detailID)%>%
  summarise(N_bursts=n_distinct(burstID),
            duration_burst=sum(duration_hrs),
            N_days=as.numeric(max(Date.and.Time)-min(Date.and.Time)),
            burst.per.day=N_bursts/N_days,
            bursttime.per.day=duration_burst/N_days)



# Fitting GLMM for sequences (simple) -------------------------------------

#The max days (length of season) until ice off starts next year
maxdDays_diff=as.numeric(difftime(ymd_hms("2020-07-02 00:00:00"),ymd_hms("2019-11-09 00:00:00"),units="days"))


hist(log(inertiaSeqSummary$meanDepth+1.7))
hist((inertiaSeqSummary$meanODBA+0.001))

#Depth
summary(seqDepth.lN<-glmer((meanDepth+1.7)~scale(dDays_diff)+(scale(dDays_diff)||Fish_detailID),
                           weights=seqDuration_hrs,
                           data=filter(inertiaSeqSummary,meanDepth>=-1.7), #Only when not at the surface
                           family=gaussian(link="log"))) #Solves convergence issues resulting in negative predictions
summary(seqDepth.null<-glmer((meanDepth+1.7)~(1|Fish_detailID),
                             weights=seqDuration_hrs,
                             data=filter(inertiaSeqSummary,meanDepth>=-1.7), #Only when not at the surface
                             family=gaussian(link="log"))) #Solves convergence issues resulting in negative predictions
summary(seqDepth.lN)$AIC[1]-summary(seqDepth.null)$AIC[1]

#ODBA
summary(seqODBA.lN<-glmer((meanODBA+0.001)~scale(dDays_diff)+(scale(dDays_diff)||Fish_detailID),
                          weights=seqDuration_hrs,
                          data=inertiaSeqSummary, #Only when not at the surface
                          family=Gamma(link="log"))) #Solves convergence issues resulting in negative predictions
summary(seqODBA.null<-glmer((meanODBA+0.001)~(1|Fish_detailID),
                            weights=seqDuration_hrs,
                            data=inertiaSeqSummary, #Only when not at the surface
                            family=Gamma(link="log"))) #Solves convergence issues resulting in negative predictions
summary(seqODBA.lN)$AIC[1]-summary(seqODBA.null)$AIC[1]

#how many do I have to throw out?
nrow(filter(inertiaSeqSummary,meanDepth<=-1.7))
nrow(filter(inertiaSeqSummary,meanDepth<=-1.7))/nrow(inertiaSeqSummary)

#####Summarising output#####

###### Depth ######
#Natural scale intercepts
exp(coef(seqDepth.lN)$Fish_detailID[1]-coef(seqDepth.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff))
#Natural scale slopes
Dslopes<-as.data.frame(coef(seqDepth.lN)$Fish_detailID[2]/sd(inertiaSeqSummary$dDays_diff))
#Natural scale maximum rises over the whole period (170 days, diff of intercept and value at 170 days)
exp(coef(seqDepth.lN)$Fish_detailID[1]-coef(seqDepth.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff))-
  exp(coef(seqDepth.lN)$Fish_detailID[1]-coef(seqDepth.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff)+
        coef(seqDepth.lN)$Fish_detailID[2]/sd(inertiaSeqSummary$dDays_diff)*170)
#Actual predicted values in the range of data
fitDepth<-bind_cols(filter(inertiaSeqSummary,meanDepth>=-1.7),fit_resp=predict(seqDepth.lN,type="response"))
#predicted changes in depth over observation for that individual
sumfitDepth<-fitDepth%>%
  group_by(Fish_detailID)%>%summarise(start=max(fit_resp),
                                      final=min(fit_resp),
                                      r=final-start,
                                      dr=max(dDays_diff)-min(dDays_diff),
                                      s=r/dr)
mean(sumfitDepth$r)
sd(sumfitDepth$r)/sqrt(nrow(sumfitDepth))
mean(sumfitDepth$s)*100
sd(sumfitDepth$s)/sqrt(nrow(sumfitDepth))*100

###### ODBA ######
#Natural scale intercepts
exp(coef(seqODBA.lN)$Fish_detailID[1]-coef(seqODBA.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff))
#Natural scale slopes
Aslopes<-as.data.frame(coef(seqODBA.lN)$Fish_detailID[2]/sd(inertiaSeqSummary$dDays_diff))
#Natural scale maximum rises over the whole period (170 days, diff of intercept and value at 170 days)
exp(coef(seqODBA.lN)$Fish_detailID[1]-coef(seqODBA.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff)+
      coef(seqODBA.lN)$Fish_detailID[2]/sd(inertiaSeqSummary$dDays_diff)*170)-
  exp(coef(seqODBA.lN)$Fish_detailID[1]-coef(seqODBA.lN)$Fish_detailID[2]*mean(inertiaSeqSummary$dDays_diff)/sd(inertiaSeqSummary$dDays_diff))

#Actual predicted values in the range of data
fitODBA<-bind_cols(filter(inertiaSeqSummary,!is.na(atRest)),fit_resp=predict(seqODBA.lN,type="response"))
#predicted changes in depth over observation for that individual
odbaChecks<-fitODBA%>%
  group_by(Fish_detailID)%>%summarise(start=min(fit_resp),
                                      final=max(fit_resp),
                                      r=max(fit_resp)-min(fit_resp),
                                      dr=max(dDays_diff)-min(dDays_diff),
                                      s=r/dr)
mean(odbaChecks$r)
sd(odbaChecks$r)

###### Correlation in Slopes ######
scor=cor(Dslopes$`scale(dDays_diff)`,Aslopes$`scale(dDays_diff)`)
tiff(filename="figureS3_cor.tiff",height=2000,width=2666,res=300)
bind_cols(Aslopes,Dslopes)%>%
  mutate(Fish_ID=row.names(.))%>%
  rename(slope_ODBA=`scale(dDays_diff)...1`,slope_Depth=`scale(dDays_diff)...2`)%>%
  ggplot()+
  geom_point(aes(slope_ODBA,slope_Depth),size=5)+
  geom_text_repel(aes(slope_ODBA,slope_Depth,label=Fish_ID),size=5,
                  box.padding=0.65,min.segment.length=1,seed=1)+
  geom_text(aes(0.007,-0.002,label=paste("Pearson's rho = ",round(scor,2))),hjust=1,size=5)+
  xlab(expression(paste("Slope in ODBA (ms"^"-2"*"/day)")))+
  ylab("Slope in Depth (m/day)")
dev.off()



# Plotting GLMM Output ----------------------------------------------------

plotDepth<-ggplot(filter(fitDepth,meanDepth>=-1.7)%>%mutate(Fish_detailID=factor(Fish_detailID,levels=c("Mqua024","Mqua028","Mqua034","Mqua026","Mqua029"))))+
  geom_point(aes(dDays_diff,meanDepth+1.7,size=seqDuration_hrs,fill=Fish_detailID),
             shape=21,alpha=0.7,show.legend=F)+
  geom_line(aes(dDays_diff,fit_resp,group=Fish_detailID),size=3)+
  geom_line(aes(dDays_diff,fit_resp,group=Fish_detailID,color=Fish_detailID),size=2.25)+
  scale_y_reverse(expand=expansion(mult=c(0.05,0.05)),name="Depth (m)",
                  breaks=seq(0,max(fitDepth$meanDepth)+10,by=10))+
  scale_x_continuous(expand=expansion(mult=c(0.01,0)),name="Days into Ice-Covered Season",
                     limits=c(0,maxdDays_diff))+
  scale_size_continuous(range=c(2,12))+
  scale_color_viridis_d(option="C",name="Fish\nID",direction=-1)+
  scale_fill_viridis_d(option="C",name="Fish ID",direction=-1)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
        axis.title.y=element_text(margin=margin(0,9,0,0)),
        legend.position=c(0.75,0.25),legend.direction="horizontal",legend.background=element_rect(color="black"))+
  guides(color=guide_legend(nrow=3,override.aes=list(size=10)))
tiff(filename="figure4a_depthglmm.tiff",width=3500,height=2000,res=300)
plotDepth
dev.off()

#Gamma
plotODBA<-ggplot()+
  geom_point(data=filter(fitODBA,meanODBA>0)%>%mutate(Fish_detailID=factor(Fish_detailID,levels=c("Mqua024","Mqua028","Mqua034","Mqua026","Mqua029"))),
             aes(dDays_diff,meanODBA,size=seqDuration_hrs,fill=Fish_detailID),
             shape=21,alpha=0.7,show.legend=F)+
  geom_hline(yintercept=0.06,lty=2,size=1)+
  geom_line(data=fitODBA,aes(dDays_diff,fit_resp,group=Fish_detailID),size=3)+
  geom_line(data=fitODBA,aes(dDays_diff,fit_resp,color=Fish_detailID,group=Fish_detailID),size=2.25)+
  scale_y_continuous(expand=expansion(add=c(0.05,0.05)),name=expression(paste("ODBA (ms"^"-2"*")")),
                     breaks=c(0.001,1,2,3),labels=c(0,1,2,3))+
  scale_x_continuous(expand=expansion(mult=c(0.01,0)),name="Days Since Ice-Covered Season Start",
                     limits=c(0,maxdDays_diff))+
  scale_size_continuous(range=c(2,12))+
  scale_color_viridis_d(option="C",name="Fish ID",direction=-1)+
  scale_fill_viridis_d(option="C",name="Fish ID",direction=-1)+
  theme(legend.position="none")
tiff(filename="figure4b_odbaglmm.tiff",width=3500,height=2000,res=300)
plotODBA
dev.off()




