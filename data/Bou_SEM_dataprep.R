library(sf)
library(here)
library(readxl)
library(psych)
library(tidyverse)
####################
##load landcover
####################
burn<-read_csv(here::here("LandCoverAttributes","PercentBurned.csv"))%>%
  rename(burn.p=PERCENTAGE)%>%
  select(Name, burn.p)

disturb<-read_csv(here::here("LandCoverAttributes","PercentDisturbance.csv"))%>%
  rename(disturb.p=PERCENTAGE)%>%
  select(Name, disturb.p)

wetl<-read_csv(here::here("LandCoverAttributes","PercentWetland.csv"))%>%
  rename(wetl.p=PERCENTAGE)%>%
  select(Name, WET_FINAL, wetl.p)%>%
  spread(WET_FINAL,wetl.p)%>%
  select(-4)%>%
  rename(upland.p=`0`,
         wetland.p=`1`)

ndvi <- read_csv(here::here("NDVI.csv"))%>%
  rename(Name=`Survey Area`)%>%
  group_by(Name)%>%
    summarise(dNDVI=mean(dNDVI),
              NDVIsummer=mean(NDVIsummer))

####################
##load wolf data
####################
wf <-st_read(here::here("WolfSurveyShapes", "WolfSurveys_ForClayton.shp"))
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

###Apply weighting and collaps
moose <- moose%>%
  mutate(d=Moose.Density*100*p.area)%>%
  group_by(Survey.Area)%>%
  summarise(Moose.Density=sum(d))

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
  rename("lambda" = "Lambda (S/(1-R/2))")%>%
  rename(survival=S,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.tweed <- read_csv(here::here("DemographyForms", "DemographyForm NEBC_Tweeds.csv"))%>%
  rename("lambda" = "Lambda (S/(1-R/2))")%>%
  rename(survival=S,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.nwt <- read_csv(here::here("DemographyForms", "Demography_NWT_NicUpdate.csv"))%>%
  rename("lambda" = "Lambda (S/(1-R/2))")%>%
  filter(!lambda %in% "#VALUE!")%>%
  rename(survival=`S (km est)`,
         reproduction=`R/2 (female calves:cow)`)%>%
  select("Area", "Year","lambda", "survival", "reproduction", "n_survival")

b.dem <- rbind(b.ab, b.sk, b.tweed, b.nwt)%>%
  drop_na(Area)

b.dem$Year <- as.numeric(b.dem$Year)
b.dem$lambda <- as.numeric(b.dem$lambda)
b.dem$reproduction <- as.numeric(b.dem$reproduction)
b.dem$survival <- as.numeric(b.dem$survival)


##load bou lookup table
b.lookup <- read_csv(here::here("DemographyForms", "bou_lookup.csv"))
b.lookup$surv <- NA
b.lookup$repro <- NA
b.lookup$weight_sum <- NA


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


b.lam <- b.lookup%>%
  group_by(wolf)%>%
  summarise(lambda=weighted.mean(x=lambda,w=w,na.rm=TRUE),
            survival=weighted.mean(x=surv,w=w,na.rm=TRUE),
            reproduction=weighted.mean(x=repro,w=w,na.rm=TRUE),
            weight_sum=weighted.mean(x=weight_sum,w=w,na.rm=TRUE),)%>%
  rename(Name=wolf)



###create final clean data
final <- burn%>%
  left_join(ndvi, by="Name")%>%
  left_join(disturb, by="Name")%>%
  left_join(wetl, by="Name")%>%
  left_join(moose%>%rename(Name=Survey.Area)%>%select(Name, Moose.Density), by="Name")%>%
  left_join(wf%>%as_tibble()%>%select(Name,WolfDensit), by="Name")%>%
  left_join(b.lam, by="Name")

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




############
###MAP
############
library(rnaturalearth)
library(rmapshaper)
##Map 
##make smaller
us <- ne_countries(country="United States of America", scale='medium',returnclass = 'sf')%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
can <- read_sf(here::here("canada", "canada.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

bor<-read_sf(here::here("boreal", "NABoreal.shp"))%>%
  st_transform("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")%>%
  filter(TYPE%in%"BOREAL")%>%
  group_by(TYPE) %>% 
  summarise(m = sum(HA)) %>% 
  st_cast()%>%
  rmapshaper::ms_simplify(input = as(bor, 'Spatial')) %>%
  st_as_sf()

ggplot() +
  geom_sf(data = us, fill ="grey95") +
  geom_sf(data = can, fill ="grey95") +
  geom_sf(data=bor2, fill="forestgreen", col=NA, alpha=0.6)+
  geom_sf(data=wf%>%filter(Name!="Tweedsmuir"), color="black", fill="grey50")+
  #scale_fill_distiller(palette = "Spectral")+
  #scale_fill_viridis_c(option="magma")+
  coord_sf(xlim =c(-2.4E6,0.2E6), ylim = c(0.8E6, 2.9E6),
           crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") +
  theme(panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue"), 
        panel.border = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        legend.position = "top")
