Clayton T. Lamb
14 October, 2020

\#\#Load Data, Functions and Cleanup Data

``` r
library(raster)
library(sf)
library(mapview)
library(here)
library(lme4)
library(sjPlot)
library(tidyverse)

##load cutblocks
cb<- st_read(here::here("data", "cutblocks", "cutblocks.shp"))
```

    ## Reading layer `cutblocks' from data source `/Users/clayton.lamb/Google Drive/Documents/University/Work/Serrouya_BouPathway/borealcaribou-pathanalysis/data/cutblocks/cutblocks.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 25843 features and 2 fields
    ## geometry type:  MULTIPOLYGON
    ## dimension:      XY
    ## bbox:           xmin: 1129500 ymin: 1361500 xmax: 1973681 ymax: 1702110
    ## epsg (SRID):    NA
    ## proj4string:    +proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +units=m +no_defs

``` r
##load dVI
rastlist <- list.files(path = here::here("data",  "dvi_annual_500m"), pattern='.tif', all.files=TRUE, full.names=TRUE)

order <- str_split(list.files(path = here::here("data",  "dvi_annual_500m"), pattern='.tif', all.files=TRUE),
            ".tif", simplify=TRUE)[,1]%>%
    as.numeric()

new.order <- tibble(reord=1:length(order),
       current=order)%>%
  arrange(current)
  
stack <- stack(rastlist[new.order$reord])

##load water mask
water <- raster(here::here("data", "water.tif"))

##remove water from dVI
stack<- stack*water

##fix up layers
for(i in 1:nlayers(stack)){
  values(stack[[i]])[values(stack[[i]])<0 & !is.na(values(stack[[i]]))] <- 0
}

##rename
names(stack) <- paste0("X", 2000:2019)

##plot
plot(stack)
```

![](README_files/figure-gfm/Load%20Data-1.png)<!-- -->

\#\#filter cutblocks to years of interest, and those of appropriate size

``` r
cb.year <-  cb%>%
  filter(HARVEST_YE %in% c(2003:2016))%>%
  mutate(area=st_area(.)%>%as.numeric())%>%
  filter(area>250000)

mean(st_area(cb.year))
```

    ## 615503.2 [m^2]

``` r
quantile(st_area(cb.year), 0.05)
```

    ## 260983.2 [m^2]

``` r
quantile(st_area(cb.year), 0.95)
```

    ## 1644879 [m^2]

\#\#Map

``` r
can <- st_read(here::here("data", "canada", "canada.shp"))%>%
  st_transform(st_crs(cb.year))
```

    ## Reading layer `canada' from data source `/Users/clayton.lamb/Google Drive/Documents/University/Work/Serrouya_BouPathway/borealcaribou-pathanalysis/data/canada/canada.shp' using driver `ESRI Shapefile'
    ## Simple feature collection with 13 features and 6 fields
    ## geometry type:  MULTIPOLYGON
    ## dimension:      XY
    ## bbox:           xmin: -141.0021 ymin: 41.68797 xmax: -52.61917 ymax: 83.11506
    ## epsg (SRID):    4326
    ## proj4string:    +proj=longlat +datum=WGS84 +no_defs

``` r
plot(st_geometry(cb.year),  border = 'grey', axes = TRUE)
plot(st_geometry(can), pch = 3, col = 'grey', add = TRUE)
plot(st_geometry(cb.year),  border = 'red', axes = TRUE, add = TRUE)
```

![](README_files/figure-gfm/map-1.png)<!-- -->

\#\#extract dVI to cublocks

``` r
a <- raster::extract(stack,as(cb.year,"Spatial"), fun=mean, na.rm=TRUE, method="simple",df=TRUE)
    
df <- cbind(as.data.frame(cb.year),a)%>%
  as.data.frame()%>%
  dplyr::select(-geometry)%>%
  gather(year,dvi,-HARVEST_YE,-ID,-prov)%>%
  mutate(year=str_sub(year,2,5)%>%as.numeric)%>%
  mutate(time=year-HARVEST_YE)
```

\#\#Plot Raw Data

``` r
precut <- df%>%
  drop_na(time)%>%
  ungroup()%>%
  filter(time%in%c(-10:-1))%>%
  group_by(ID)%>%
  summarise(baseline=mean(dvi, na.rm=TRUE))
  
contrast <- df%>%
    drop_na(time)%>%
  left_join(precut, by="ID")%>%
  drop_na(baseline)%>%
  ungroup()%>%
  group_by(ID)%>%
  mutate(change=((dvi-baseline)/baseline)*100)%>%
  ungroup()


 contrast%>%
  filter(time>=-10 & time<18)%>%
ggplot(aes(x=time, y=change))+
  geom_vline(xintercept = 0, linetype="dotted")+
  geom_point(col="red", alpha=0.01)+
  theme_bw()+
  xlab("time since cut (years)")+
  ylab("change from pre-cut dVI (%)")+
   facet_grid(.~prov)
```

![](README_files/figure-gfm/plot-1.png)<!-- -->

``` r
 contrast%>%
   filter(time>=-10 & time<18)%>%
   ggplot(aes(x=time, y=change))+
   geom_vline(xintercept = 0, linetype="dotted")+
   geom_point(col="red", alpha=0.01)+
   theme_bw()+
   xlab("Time since cut (years)")+
   ylab("Change from pre-cut vegetation index (%)")
```

![](README_files/figure-gfm/plot-2.png)<!-- -->

\#\#Test for significant effect post logging

``` r
##prep data
model.data <-
  contrast%>%
  mutate(period=case_when(time<0~"Pre",
                          time>0~"Post"))%>%
  mutate(period=fct_relevel(period, "Pre", "Post"))%>%
  drop_na(period, time)

##run model with dVI values
m1 <- lmer(dvi ~ period + (1|ID), data=model.data)
summary(m1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: dvi ~ period + (1 | ID)
    ##    Data: model.data
    ## 
    ## REML criterion at convergence: 863512
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.8448 -0.5884  0.0040  0.6068  5.4604 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 242069   492.0   
    ##  Residual             237661   487.5   
    ## Number of obs: 56164, groups:  ID, 2956
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 1889.062      9.610  196.58
    ## periodPost   359.517      4.612   77.95
    ## 
    ## Correlation of Fixed Effects:
    ##            (Intr)
    ## periodPost -0.260

``` r
##plot model
plot_model(m1,
  sort.est = TRUE,
  title = "",
type="pred",
  axis.title=c("deltaVI"),
ci.lvl=0.95)
```

    ## $period

![](README_files/figure-gfm/mixed%20model-1.png)<!-- -->

``` r
model.data%>%
  group_by(period)%>%
  summarise(mean(dvi))

##using the change
m2 <- lmer(change ~ period + (1|ID), data=model.data)
summary(m2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: change ~ period + (1 | ID)
    ##    Data: model.data
    ## 
    ## REML criterion at convergence: 543106.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -8.1880 -0.5746 -0.0070  0.5841  6.4688 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 299.6    17.31   
    ##  Residual             831.9    28.84   
    ## Number of obs: 56164, groups:  ID, 2956
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   1.6182     0.3708   4.364
    ## periodPost   23.7137     0.2701  87.791
    ## 
    ## Correlation of Fixed Effects:
    ##            (Intr)
    ## periodPost -0.394

``` r
plot_model(m2,
  sort.est = TRUE,
  title = "",
type="pred",
  axis.title=c(" % change in deltaVI"),
ci.lvl=0.95)
```

    ## $period

![](README_files/figure-gfm/mixed%20model-2.png)<!-- -->

\#\#Bootstrap results for a time-specific plot

``` r
boot.dat <- data.frame()
for(i in 1:500){
dat <-  contrast%>%
  group_by(ID)%>%
  sample_frac(1,replace = TRUE)%>%
  filter(time>=-10 & time<14)%>%
  drop_na(time)

mod <- dat%>%
  lmer(dvi ~ as.character(time) + (1|ID), data=.)

dat$pred <-predict(mod)



boot.dat <- dat%>%
  group_by(time)%>%
  summarise(av=mean(dvi,na.rm=TRUE))%>%
  mutate(iter=i)%>%
  rbind(boot.dat)
}



boot.dat%>%
   group_by(time)%>%
   summarise(mean=mean(av, na.rm=TRUE),
             lower=quantile(av, 0.05),
             upper=quantile(av,0.95))%>%
   ggplot(aes(x=time, y=mean))+
   geom_vline(xintercept = 0, linetype="dotted")+
   geom_errorbar(aes(ymin = lower, ymax = upper), width=0.01, alpha=0.2)+
   geom_point(col="red")+
   theme_bw()+
   xlab("Time since cut (years)")+
   ylab("Change from pre-cut vegetation index (%)")
```

![](README_files/figure-gfm/bootstrap-1.png)<!-- -->

``` r
precut.boot <- boot.dat%>%
  ungroup()%>%
  filter(time%in%c(-10:-1))%>%
  group_by(iter)%>%
  summarise(baseline=mean(av, na.rm=TRUE))
  
contrast.boot <- boot.dat%>%
  ungroup()%>%
  left_join(precut.boot, by="iter")%>%
  mutate(change=((av-baseline)/baseline)*100)

  
   contrast.boot%>%
     group_by(time)%>%
   summarise(mean=mean(change, na.rm=TRUE),
             lower=quantile(change, 0.05),
             upper=quantile(change,0.95))%>%
   ggplot(aes(x=time, y=mean))+
   geom_vline(xintercept = 0, linetype="dotted")+
   geom_errorbar(aes(ymin = lower, ymax = upper), width=0.01, alpha=0.2)+
   geom_point(col="red")+
   theme_bw()+
   theme(panel.grid.minor = element_blank())+
   xlab("Time since cut (years)")+
   ylab("Change from pre-cut vegetation index (%)")
```

![](README_files/figure-gfm/bootstrap-2.png)<!-- -->

``` r
    ggsave(here::here("plots", "time_since_cut_dEVI.png"),width=6, height=4)
```
