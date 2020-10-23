library(sf)
library(here)
library(readxl)
library(psych)
library(raster)
library(stringr)
library(lwgeom)
library(rworldmap)
library(rnaturalearth)
library(rmapshaper)
library(ggspatial)
library(ggsflabel)
library(tidyverse)

####################
##load landcover
####################
burn<-read_csv(here::here("LandCoverAttributes","PercentBurned.csv"))%>%
  rename(burn.p=PERCENTAGE)%>%
  dplyr::select(Name, burn.p)%>%
  mutate(Name=case_when(Name%in%"Hay River Lowlands"~"Fort Providence HRL",
                        TRUE~Name))

disturb<-read_csv(here::here("LandCoverAttributes","PercentDisturbance.csv"))%>%
  rename(disturb.p=direct.disturb)%>%
  dplyr::select(Name, disturb.p)%>%
  mutate(Name=case_when(Name%in%"Hay River Lowlands"~"Fort Providence HRL",
                        TRUE~Name))%>%
  rbind(tibble(
    Name=c("Whati (TASR Impact)","Jean Marie River"),
    disturb.p=c(0.015,0.031)))

wetl<-read_csv(here::here("LandCoverAttributes","PercentWetland.csv"))%>%
  rename(wetl.p=PERCENTAGE)%>%
  dplyr::select(Name, WET_FINAL, wetl.p)%>%
  spread(WET_FINAL,wetl.p)%>%
  select(-4)%>%
  rename(upland.p=`0`,
         wetland.p=`1`)%>%
  mutate(Name=case_when(Name%in%"Hay River Lowlands"~"Fort Providence HRL",
                        TRUE~Name))

####################
##load wolf data
####################
wf <- st_read(here::here("WolfSurveyShapes", "WolfSurveyAreas_ForPaper_FinalNoTweeds.shp"))

##add in new NWT data
wf <- wf%>%
  rbind(st_read(here::here("WolfSurveyShapes", "NWT_NewCensusBlocks_29092020.shp"))%>%
          st_transform(crs=st_crs(wf))%>%
          select(Name,area_km2=Km2,YearSurvey,WolfDensit))

plot(wf["WolfDensit"])

##pull out study area names
sa <- wf$Name



####################
##load moose data
####################
moose<-read.csv(here::here("MooseDensities.csv"))


##load moose shps for AB
mshp <- st_read(here::here("WMUs_Moose", "Alberta_WMUs_Moose.shp"))%>%
  st_transform("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")%>%
  select(WMU)%>%
  filter(WMU %in% c(529,512,519,517))

twdshp <- st_read(here::here("WMUs_Moose", "TweedsMooseSurveyDiss.shp"))%>%
  st_transform(st_crs(mshp)[2]$proj4string)

#plot(st_geometry(mshp))
#plot(st_geometry(twdshp))

##get weighted area means ready
moose$p.area<-1
moose$Survey.Area <- as.character(moose$Survey.Area)
###get wolf data in same proj
wf <- wf%>%st_transform(st_crs(mshp)[2]$proj4string)

wf%>%
  filter(Name%in%"Cold Lake")%>%
  st_intersection(mshp)%>%
  mutate(p=st_area(.)/sum(st_area(.)))


moose[moose$Survey.Area%in%"Cold Lake AB 529","p.area"] <-0.0649
moose[moose$Survey.Area%in%"Cold Lake AB 519","p.area"] <-0.0562
moose[moose$Survey.Area%in%"Cold Lake AB 517","p.area"] <-0.6090
moose[moose$Survey.Area%in%"Cold Lake AB 512","p.area"] <-0.2697

##edit names
moose[moose$Survey.Area%in%"Cold Lake AB 529","Survey.Area"] <-"Cold Lake"
moose[moose$Survey.Area%in%"Cold Lake AB 519","Survey.Area"] <-"Cold Lake"
moose[moose$Survey.Area%in%"Cold Lake AB 517","Survey.Area"] <-"Cold Lake"
moose[moose$Survey.Area%in%"Cold Lake AB 512","Survey.Area"] <-"Cold Lake"




wf%>%
  filter(Name%in%"Tweedsmuir")%>%
  st_intersection(twdshp)%>%
  mutate(p=st_area(.)/sum(st_area(.)))


moose[moose$Survey.Area%in%"Tweedsmuir Survey Area North","p.area"] <-0.477
moose[moose$Survey.Area%in%"Tweedsmuir Survey Area South","p.area"] <-0.522


##edit names
moose[moose$Survey.Area%in%"Tweedsmuir Survey Area North","Survey.Area"] <-"Tweedsmuir"
moose[moose$Survey.Area%in%"Tweedsmuir Survey Area South","Survey.Area"] <-"Tweedsmuir"

###Apply weighting and collapse
moose <- moose%>%
  mutate(d=Moose.Density*100*p.area)%>%
  group_by(Survey.Area)%>%
  summarise(Moose.Density=sum(d), Moose_Year_cl=mean(Moose_Year_cl))

##################
##load caribou data
##################
b.rang <- st_read(here::here("Boreal_Caribou_Range_Boundaries_AsOfJune62012", "Boreal_Caribou_Range_Boundaries_AsOfJune62012.shp"))%>%
  st_transform(st_crs(mshp)[2]$proj4string)


b.ab <- read_csv(here::here("DemographyForms", "2019_07_CLAR_CLAB_CLSK_ESAR_YATE_Demographic_Update_Preliminary.csv"))%>%
  rename("lambda" = "Lambda_Mean")%>%
  filter(!lambda %in% "#VALUE!")%>%
  rename(survival=`Surv Mean`,
         reproduction=Recruitment_Mean,
         n_survival=`n-KM`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.sk <- read_csv(here::here("DemographyForms", "Demography_form_SaskatchewanBorealShield.csv"))%>%
  rename("lambda" = "Lambda(S/(1-Rrm))")%>%
  rename(survival=S,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.tweed <- read_csv(here::here("DemographyForms", "DemographyForm NEBC_Tweeds.csv"))%>%
  rename("lambda" = "Lambda(S/(1-Rrm))")%>%
  rename(survival=S,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.nwt <- read_csv(here::here("DemographyForms", "Demography_NWT_NicUpdate.csv"))%>%
  rename("lambda" = "Lambda(S/(1-Rrm))")%>%
  filter(!lambda %in% "#VALUE!")%>%
  rename(survival=`S (km est)`,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.nwt2 <- read_csv(here::here("DemographyForms", "Demography_NWT_New.csv"))%>%
  rename("lambda" = "Lambda(S/(1-Rrm))")%>%
  filter(!lambda %in% "#VALUE!")%>%
  rename(survival=`S (km est)`,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.dem <- rbind(b.ab, b.sk, b.tweed, b.nwt, b.nwt2)%>%
  drop_na(Area)

b.dem$Year <- as.numeric(b.dem$Year)
b.dem$lambda <- as.numeric(b.dem$lambda)
b.dem$reproduction <- as.numeric(b.dem$reproduction)
b.dem$survival <- as.numeric(b.dem$survival)


##load bou lookup table
b.lookup <- read_csv(here::here("DemographyForms", "bou_lookup.csv"))
b.lookup$lambda <- NA_real_
b.lookup$surv <- NA_real_
b.lookup$repro <- NA_real_
b.lookup$weight_sum <- NA_real_


for(i in 1:nrow(b.lookup)){
  a <- b.lookup[i,]
  b <- b.dem%>%filter(Area %in% a$bousheet[1] &
                        Year %in% c((a$Year[1]-2),
                                    (a$Year[1]-1),
                                    (a$Year[1])))
  
  b.lookup[i, "weight_sum"] <- sum(b$n_survival,na.rm=TRUE)
  b.lookup[i, "lambda"] <- weighted.mean(x=b$lambda,w=b$n_survival/b.lookup[i, "weight_sum"]%>%pull(weight_sum),na.rm=TRUE)
  b.lookup[i, "surv"] <- weighted.mean(x=b$survival,w=b$n_survival/b.lookup[i, "weight_sum"]%>%pull(weight_sum),na.rm=TRUE)
  b.lookup[i, "repro"] <- weighted.mean(x=b$reproduction,w=b$n_survival/b.lookup[i, "weight_sum"]%>%pull(weight_sum),na.rm=TRUE)
  
}

b.lookup$lambda <- as.numeric(b.lookup$lambda)

##average the Jean Marie River blocks by space
JMR.av <- weighted.mean(x=c(b.lookup[b.lookup$bousheet=="Jean Marie River South","lambda"][[1]],b.lookup[b.lookup$bousheet=="Jean Marie River North","lambda"][[1]]),
                                                                                w=c(b.lookup[b.lookup$bousheet=="Jean Marie River South","w"][[1]],b.lookup[b.lookup$bousheet=="Jean Marie River North","w"][[1]]))

b.lookup[b.lookup$bousheet=="Jean Marie River South","lambda"] <- JMR.av
b.lookup[b.lookup$bousheet=="Jean Marie River North","lambda"] <- JMR.av

b.lam <- b.lookup%>%
  group_by(wolf)%>%
  summarise(lambda=weighted.mean(x=lambda,w=w,na.rm=TRUE),
            survival=weighted.mean(x=surv,w=w,na.rm=TRUE),
            reproduction=weighted.mean(x=repro,w=w,na.rm=TRUE),
            weight_sum=weighted.mean(x=weight_sum,w=w,na.rm=TRUE))%>%
  rename(Name=wolf)


##################
##load delta vegetation index
##################
rastlist <- list.files(path = here::here("dvi_annual_500m"), pattern='.tif', all.files=TRUE, full.names=TRUE)

allrasters <- stack(rastlist)

water <- raster(here::here("water.tif"))%>%
  projectRaster(allrasters[[1]])%>%
  is.na()


values(water)[values(water)==0]<-2
values(water)<- values(water)-1
values(water)[is.na(values(water))]<-0
plot(water)

allrasters<- allrasters*water


##fix up layers
for(i in 1:nlayers(allrasters)){
  values(allrasters[[i]])[values(allrasters[[i]])<0] <- 1
  values(allrasters[[i]])[values(allrasters[[i]])==0]<-NA
}



plot(allrasters)


moose <- subset(moose, !Survey.Area %in% "Tweedsmuir")%>%
  rbind(tibble(
    Survey.Area=c("Whati (TASR Impact)","Jean Marie River"),
    Moose.Density=c(1.1,4.5),
    Moose_Year_cl=2019))

wf <- wf%>%
  mutate(Name=case_when(Name%in%c("Jean Marie River South", "Jean Marie River North")~"Jean Marie River",
         TRUE~as.character(Name)))%>%
    group_by(Name)%>%
    summarise(WolfDensit=mean(WolfDensit))

plot(wf)

for(i in 1:nrow(moose)){
  start.year <- moose[i, "Moose_Year_cl"]-2003
  
  a <- allrasters[[start.year[[1]]:(start.year[[1]]+4)]]
  
  b <- raster::mask(a,as(wf[ wf$Name%in%moose[i, "Survey.Area"],]%>%st_zm(),"Spatial"))
    
  moose[i, "LAI"] <- mean(c(mean(values(b[[1]]), na.rm=TRUE),
                            mean(values(b[[2]]), na.rm=TRUE),
                            mean(values(b[[3]]), na.rm=TRUE),
                            mean(values(b[[4]]), na.rm=TRUE),
                            mean(values(b[[5]]), na.rm=TRUE)))
  print(i)
}
                    


###create final clean data

###WHY DOESNT THE NEW NWT data join?? ADD manually
disturb$Name[15]==b.lam$Name[c(12)] ##this cannot be false, but is.
##add again
b.lam$Name[c(12)] <- "Jean Marie River"
b.lam$Name[c(15)] <-"Whati (TASR Impact)"



final <- disturb%>%
  left_join(moose%>%dplyr::rename(Name=Survey.Area)%>%dplyr::select(Name, Moose.Density, LAI), by="Name")%>%
  left_join(wf%>%as_tibble()%>%dplyr::select(Name,WolfDensit), by="Name")%>%
  left_join(b.lam, by="Name")%>%
  rename(dEVI=LAI)%>%
  filter(Name!="Tweedsmuir")  ##remove Tweedsmuir- non-boreal


write.csv(final, here::here("final.csv"))





###calculate poly overlap for bou spatial weighting in lookup table
wf%>%filter(Name%in% "Fort Resolution Reference")%>%
  st_intersection(b.rang)%>%
  mutate(area=st_area(.))%>%
  mutate(w=area/sum(area))%>%
  select(HERD, w)

wf%>%filter(Name%in% "Cold Lake")%>%
  st_intersection(b.rang)%>%
  mutate(area=st_area(.))%>%
  mutate(w=area/sum(area))%>%
  select(HERD, w)

wf%>%filter(Name%in% "Clarke")%>%
  st_intersection(b.rang)%>%
  mutate(area=st_area(.))%>%
  mutate(w=area/sum(area))%>%
  select(HERD, w)


wf%>%filter(Name%in% "Jean Marie River South")%>%st_area()/ (wf%>%filter(Name%in% "Jean Marie River South")%>%st_area()+wf%>%filter(Name%in% "Jean Marie River North")%>%st_area())

############
###MAP
############


bou.ranges <- b.rang%>%
  filter(HERD%in%c("Chinchaga","Calendar","Yates","Cold Lake","Saskatchewan Boreal Plains","Saskatchewan Boreal Shield","East Side Athabasca River" ,"Snake-Sahtahneh"))%>%
  select(HERD)%>%
  rbind(st_read(here::here("Boreal_Caribou_Range_Boundaries_AsOfJune62012", "Boreal_Caribou_Revised_Study_Areas_2018.shp"))%>%
          st_transform(st_crs(b.rang))%>%
          filter(Name%in%c("Hay River Lowlands","Dehcho South","Dehcho North", "North Slave"))%>%
          select(HERD=Name))%>%
  rbind(st_read(here::here("Boreal_Caribou_Range_Boundaries_AsOfJune62012", "Boreal_Caribou_Revised_Study_Areas_2018.shp"))%>%
          st_transform(st_crs(b.rang))%>%
          filter(Name%in%c("Pine Point / Buffalo Lake"))%>%
          st_buffer(20000)%>%
          select(HERD=Name))


##Map 
##make smaller
us <- ne_countries(country="United States of America", scale='medium',returnclass = 'sf')%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
can <- read_sf(here::here("canada", "canada.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

bor<-read_sf(here::here("boreal", "NABoreal.shp"))%>%
  st_make_valid()%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  filter(TYPE%in%"BOREAL")%>%
  group_by(TYPE) %>% 
  summarise(m = sum(HA)) %>% 
  st_cast()

bor<-rmapshaper::ms_simplify(input = as(bor, 'Spatial'))


bor <- read_sf(here::here("terr-ecoregions-TNC", "tnc_terr_ecoregions.shp"))%>%
  filter(WWF_MHTNAM%in%"Boreal Forests/Taiga")%>%
  st_as_sf()%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  group_by() %>% 
  summarise() %>% 
  st_cast()




world <- getMap(resolution = "high")
world <- st_as_sf(world)
##prep extents
extent = matrix(c(-2.4E6 ,0.8E6 , -2.4E6 , 2.9E6 , 0.2E6  , 2.9E6, 0.2E6 ,0.8E6 , -2.4E6 ,0.8E6 ),ncol=2, byrow=TRUE)
pts = list(extent)
pl1 = st_polygon(pts)
pl1 <- st_sfc(pl1, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  st_transform("+proj=laea +lat_0=40 +lon_0=-100 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")

inset <- ggplot() +
  geom_sf(data = world,size=0.1) +
  geom_sf(data=bor%>%st_transform("+proj=laea +lat_0=52 +lon_0=-100 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "), fill="forestgreen", col=NA, alpha=0.6)+
  geom_sf(data = pl1, fill=NA, color="red", size = 0.4) +
  coord_sf(crs = "+proj=laea +lat_0=40 +lon_0=-100 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")+
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed", size = 0.1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(fill = NA, color = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = c(0.65,0.075),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.background = element_blank())


ggplot() +
  geom_sf(data = us, fill ="grey95") +
  geom_sf(data = can, fill ="grey95") +
  geom_sf(data=b.rang, fill="grey75", col=NA)+
  geom_sf(data=bou.ranges, fill="grey35", col=NA)+
  geom_sf(data=wf%>%filter(!Name%in%c("Tweedsmuir"))%>%left_join(disturb, by="Name"), aes(fill=disturb.p), color="black")+
  ggsflabel::geom_sf_text(data=can%>%filter(!PROV%in%c("NT", "SK")),aes(label = PROV), color="black")+
  ggsflabel::geom_sf_text(data=can%>%filter(PROV%in%"NT"),aes(label = PROV), color="black", nudge_y=-100000, nudge_x=150000)+
  ggsflabel::geom_sf_text(data=can%>%filter(PROV%in%"SK"),aes(label = PROV), color="black", nudge_y=-200000)+
  ggsflabel::geom_sf_text_repel(data=wf%>%filter(!Name%in%c("Tweedsmuir"))%>%mutate(label=1:14),aes(label = label), color="black",min.segment.length=0.01, force=10)+
  scale_fill_viridis_c(guide = guide_colourbar(direction = "horizontal"),
                       name = "Habitat alteration (%)")+
  coord_sf(xlim =c(-2.6E6,0.2E6), ylim = c(0.8E6, 2.9E6),
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  theme(panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = c(0.65,0.075),
        legend.background = element_blank())+
  annotation_scale(location = "bl", width_hint = 0.25)+
  annotation_custom(ggplotGrob(inset), xmin =-2.99E6, xmax = -1.79E6, ymin = 0.83E6, ymax = 1.53E6)



ggsave("/Users/clayton.lamb/Google Drive/Documents/University/Work/Serrouya_BouPathway/borealcaribou-pathanalysis/plots/map.png", width=5, height=4, units="in")


