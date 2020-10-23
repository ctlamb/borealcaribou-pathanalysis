Caribou Path Analysis
================
Clayton T. Lamb
22 October, 2020

\#\#Load Data, Functions and Cleanup Data

``` r
##packages
library(here)
library(lavaan)
library(semPlot)
require(cowplot)
library(corrplot)
library(ggraph)
library(igraph)
library(QuantPsyc)
library(ggpubr)
library(MuMIn)
library(knitr)
library(piecewiseSEM)
library(tidyverse)

##data
df <- read.csv(here::here("data", "final.csv"))

##transform to instantaneous rate of growth (r)
df$caribou.lambda <- log(df$lambda)

##set seed for any bootstrapping
set.seed(2019)
```

\#\#Find intersections

``` r
###Wolf-Caribou lambda<0 intersection
wf_int <-predict(lm(WolfDensit~caribou.lambda+I(caribou.lambda^2)+I(caribou.lambda^3), data=df), newdata=data.frame(caribou.lambda=0))

##plot
df$predicted <- predict(lm(WolfDensit~caribou.lambda+I(caribou.lambda^2)+I(caribou.lambda^3), data=df), newdata=df)
ggplot(df%>%arrange(predicted))+
  geom_point(aes(x=WolfDensit, y=caribou.lambda))+
  geom_line(aes(y=caribou.lambda, x=predicted))+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept=wf_int, linetype="dashed", col="red")+
  xlab(expression(wolf~(n/1000~km^2)))+
  ylab("caribou pop. growth (r)")
```

![](README_files/figure-gfm/Find%20intersections-1.png)<!-- -->

``` r
##population declines when > 1.8 wolves/1000 sq.km


##Moose-Wolf intersection

##1.8 wolves/1000 sq.km generally reached when moose are greater than 2.8/100 sq.km
ms_int <- predict(lm(Moose.Density~WolfDensit+I(WolfDensit^2), data=df), newdata=data.frame(WolfDensit=wf_int))

##plot
df$predicted.moose <- predict(lm(Moose.Density~WolfDensit+I(WolfDensit^2), data=df), newdata=df)
ggplot(df)+
  geom_point(aes(y=WolfDensit, x=Moose.Density))+
  geom_line(aes(y=WolfDensit, x=predicted.moose))+
  geom_vline(xintercept=ms_int, linetype="dashed", col="red")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  xlab(expression(moose~(n/100~km^2)))+
  ylab(expression(wolf~(n/1000~km^2)))
```

![](README_files/figure-gfm/Find%20intersections-2.png)<!-- -->

``` r
#bootstrap errors
int.data <- data.frame()
for(i in 1:1000){
df.i <- df%>%sample_frac(1, replace=TRUE)

a <- data.frame(i=i,
                type=c("wolf", "moose"),
                val=c(predict(lm(WolfDensit~caribou.lambda+I(caribou.lambda^2)+I(caribou.lambda^3), data=df.i), newdata=data.frame(caribou.lambda=0)),
                      predict(lm(Moose.Density~WolfDensit+I(WolfDensit^2), data=df.i), newdata=data.frame(WolfDensit=1.9))))
              
int.data <- rbind(int.data,a)
}

int.data%>%
  group_by(type)%>%
  summarize(median=median(val)%>%round(2),
            upper=quantile(val,0.95)%>%round(2),
            lower=quantile(val,0.05)%>%round(2))%>%
  print()
```

\#\#transformations to linear

``` r
###transform caribou lambda
df$caribou.lambda.t <--(exp(-10*df$caribou.lambda))

##backtransform w/  log(-df$caribou.lambda.t)/-10

##make sure the rest remain linear 
a <- ggplot(df, aes(x=WolfDensit, y=caribou.lambda.t))+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

b <- ggplot(df, aes(x=Moose.Density, y=caribou.lambda.t))+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

c <- ggplot(df, aes(x=dEVI, y=caribou.lambda.t))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

ggarrange(a,b,c,
          ncol = 3, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-gfm/Transform-1.png)<!-- -->

``` r
###transform dEVI
df$dEVI.t <-exp(0.005*df$dEVI)
##backtransform w/  log(df$dEVI.t)/0.005
##make sure the rest remain linear
d <- ggplot(df, aes(y=Moose.Density, x=dEVI.t))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

e <- ggplot(df, aes(y=WolfDensit, x=dEVI.t))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

f <- ggplot(df, aes(y=caribou.lambda.t, x=dEVI.t))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

ggarrange(d,e,f,
          ncol = 3, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-gfm/Transform-2.png)<!-- -->

\#\#Plot raw data

``` r
a <-ggplot(df, aes(x=disturb.p, y=dEVI))+
    geom_line(data=df%>%mutate(pred=predict(lm(dEVI.t~disturb.p,data=.))%>%log()/.005), aes(x=disturb.p, y=pred), colour="grey")+
  geom_point()+
  theme_bw()+
  ylab("Vegetation index")+
  xlab("Habitat alteration")+
  theme(panel.grid.minor = element_blank())

b <-ggplot(df, aes(x=dEVI, y=Moose.Density))+
    geom_line(data=df%>%mutate(pred=predict(lm(Moose.Density~dEVI.t, data=.))), aes(x=dEVI, y=pred), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  xlab("Vegetation index")+
  ylab(expression(Moose~(ind./100~km^2)))+
  scale_y_continuous( limits = c(0,NA))+
  scale_x_continuous(limits = c(750,NA))

c <-ggplot(df, aes(x=Moose.Density, y=WolfDensit))+
    geom_line(data=df%>%mutate(pred=predict(lm(WolfDensit~Moose.Density, data=.))), aes(x=Moose.Density, y=pred), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  xlab(expression(Moose~(ind./100~km^2)))+
  ylab(expression(Wolf~(ind./1000~km^2)))+
  scale_x_continuous(limits = c(0,NA))

d <-ggplot(df, aes(x=WolfDensit, y=exp(caribou.lambda)))+
  geom_line(data=df%>%mutate(pred=-predict(lm(caribou.lambda.t~WolfDensit, data=.))%>%log()/-10), aes(x=WolfDensit, y=exp(pred)), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, linetype="dashed")+
  xlab(expression(Wolf~(ind./1000~km^2)))+
  ylab("Caribou pop. growth")

r1a <-ggplot(df, aes(x=disturb.p, y=exp(caribou.lambda)))+
  geom_line(data=df%>%mutate(pred=-predict(lm(caribou.lambda.t~disturb.p, data=.))%>%log()/-10), aes(x=disturb.p, y=exp(pred)), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, linetype="dashed")+
  xlab("Habitat Alteration")+
  ylab("Caribou pop. growth")


r1b <-ggplot(df, aes(x=dEVI, y=exp(caribou.lambda)))+
 geom_line(data=df%>%mutate(pred=-predict(lm(caribou.lambda.t~dEVI.t, data=.))%>%log()/-10), aes(x=log(dEVI.t)/.005, y=exp(pred)), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  geom_hline(yintercept=1, linetype="dashed")+
  xlab("Vegetation index")+
  ylab("Caribou pop. growth")+
  scale_x_continuous(limits = c(750,NA))

r1c <-ggplot(df, aes(x=disturb.p, y=WolfDensit))+
    geom_line(data=df%>%mutate(pred=predict(lm(WolfDensit~disturb.p, data=.))), aes(x=disturb.p, y=pred), colour="grey")+
  geom_point()+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  ylab(expression(Wolf~(ind./1000~km^2)))+
  xlab("Habitat Alteration")


ggarrange(a,b,c,d,r1a,r1b,r1c,nrow=2,ncol=4, labels ="AUTO")
```

![](README_files/figure-gfm/plot%20raw%20data-1.png)<!-- -->

``` r
ggsave(here::here("plots","univar.png"), width=7, height=5, units="in")
ggarrange(b,c,e,nrow=1,ncol=3, labels ="AUTO")
```

![](README_files/figure-gfm/plot%20raw%20data-2.png)<!-- -->

``` r
ggsave(here::here("plots","univar2.png"), width=7, height=2.5, units="in")


f <- ggplot(df, aes(x=WolfDensit, y=survival))+
  geom_point()+
  theme_bw()+
  ylab("Survival")+
  theme(panel.grid.minor = element_blank())+
  xlab(expression(wolf~(n/1000~km^2)))

g <- ggplot(df, aes(x=WolfDensit, y=reproduction))+
  geom_point()+
  theme_bw()+
  ylab("Reproduction")+
  theme(panel.grid.minor = element_blank())+
  xlab(expression(wolf~(n/1000~km^2)))



ggarrange(f,g,
          ncol = 2, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-gfm/plot%20raw%20data-3.png)<!-- -->

``` r
ggsave(here::here("plots","vitalrate_wolf.png"), width=6, height=2.7, units="in")
```

\#\#Correlation marix

``` r
M <- cor(df%>%
           dplyr::select(disturb.p,dEVI.t, Moose.Density, WolfDensit,caribou.lambda.t)%>%
           rename(`Habitat alteration`=disturb.p,
                  `Vegetation index`=dEVI.t,
                  `Moose density`=Moose.Density,
                  `Wolf density`=WolfDensit,
                  `Caribou pop. growth`=caribou.lambda.t), use="complete.obs")

corrplot(M, method = "number", type = "upper", order = "hclust")
```

![](README_files/figure-gfm/corr-1.png)<!-- -->

\#\#D-Separation analysis

``` r
##lay out models
##green=dEVI.t which is a greeness (food) index
##ha= habitat alteration, which is what we call disturb.p here as well
mA <- "green>moose>wolf, ha>caribou"
mB <- "green>moose>wolf>caribou, ha>wolf"
mC <- "green>moose>wolf, green>caribou, moose>caribou, ha"
mD <- "green>moose>wolf>caribou, ha"
mE <- "green>moose>wolf>caribou, ha>caribou"
mF <- "green>moose>wolf, green>caribou, ha"

modelA <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density, df),
               lm(caribou.lambda.t ~ disturb.p, df))

modelB <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density + disturb.p, df),
               lm(caribou.lambda.t ~ WolfDensit, df))

modelC <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density, df),
               lm(caribou.lambda.t ~ Moose.Density+dEVI.t, df))%>%
  update(disturb.p ~ 1)

modelD <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density, df),
               lm(caribou.lambda.t ~ WolfDensit, df))%>%
  update(disturb.p ~ 1)

modelE <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density, df),
               lm(caribou.lambda.t ~ WolfDensit +disturb.p, df))

modelF <- psem(lm(Moose.Density ~ dEVI.t, df),
               lm(WolfDensit ~ Moose.Density, df),
               lm(caribou.lambda.t ~ dEVI.t, df))%>%
  update(disturb.p ~ 1)





##summarize
data.frame(model=c("A","B","C","D","E","F"),
                      description=c(mA,mB,mC,mD,mE,mF),
                      p=round(c(summary(modelA, .progressBar = F)$Cstat$P.Value,
                           summary(modelB, .progressBar = F)$Cstat$P.Value,
                           summary(modelC, .progressBar = F)$Cstat$P.Value,
                           summary(modelD, .progressBar = F)$Cstat$P.Value,
                           summary(modelE, .progressBar = F)$Cstat$P.Value,
                           summary(modelF, .progressBar = F)$Cstat$P.Value
                           ),3),
                      K=c(summary(modelA, .progressBar = F)$IC$K,
                           summary(modelB, .progressBar = F)$IC$K,
                           summary(modelC, .progressBar = F)$IC$K,
                           summary(modelD, .progressBar = F)$IC$K,
                           summary(modelE, .progressBar = F)$IC$K,
                           summary(modelF, .progressBar = F)$IC$K
                           ),
                      AICc=round(c(summary(modelA, .progressBar = F)$IC$AICc,
                           summary(modelB, .progressBar = F)$IC$AICc,
                           summary(modelC, .progressBar = F)$IC$AICc,
                           summary(modelD, .progressBar = F)$IC$AICc,
                           summary(modelE, .progressBar = F)$IC$AICc,
                           summary(modelF, .progressBar = F)$IC$AICc
                           ),2))%>%
  mutate(dAICc=AICc-min(AICc))%>%
  arrange(dAICc)%>%
  as_tibble()%>%
  kable()
```

| model | description                                            |     p |  K |   AICc | dAICc |
| :---- | :----------------------------------------------------- | ----: | -: | -----: | ----: |
| D     | green\>moose\>wolf\>caribou, ha                        | 0.527 |  9 | 101.58 |  0.00 |
| B     | green\>moose\>wolf\>caribou, ha\>wolf                  | 0.461 | 10 | 138.94 | 37.36 |
| E     | green\>moose\>wolf\>caribou, ha\>caribou               | 0.450 | 10 | 139.50 | 37.92 |
| A     | green\>moose\>wolf, ha\>caribou                        | 0.036 |  9 | 140.34 | 38.76 |
| F     | green\>moose\>wolf, green\>caribou, ha                 | 0.013 |  9 | 151.77 | 50.19 |
| C     | green\>moose\>wolf, green\>caribou, moose\>caribou, ha | 0.043 | 10 | 180.91 | 79.33 |

``` r
###Is there another path (F), that was excluded but was maybe statistically important?
#lm(caribou.lambda.t~ disturb.p + WolfDensit, data=df)%>%summary() 
##no, wolf density remains significantly negative (p=0.0006), disturb.p has no effect (p=0.58)

##final model selection table with p<0.05 removed
aic.tab <- data.frame(model=c("A","B","C","D","E","F"),
                      description=c(mA,mB,mC,mD,mE,mF),
                      p=round(c(summary(modelA, .progressBar = F)$Cstat$P.Value,
                           summary(modelB, .progressBar = F)$Cstat$P.Value,
                           summary(modelC, .progressBar = F)$Cstat$P.Value,
                           summary(modelD, .progressBar = F)$Cstat$P.Value,
                           summary(modelE, .progressBar = F)$Cstat$P.Value,
                           summary(modelF, .progressBar = F)$Cstat$P.Value
                           ),3),
                      K=c(summary(modelA, .progressBar = F)$IC$K,
                           summary(modelB, .progressBar = F)$IC$K,
                           summary(modelC, .progressBar = F)$IC$K,
                           summary(modelD, .progressBar = F)$IC$K,
                           summary(modelE, .progressBar = F)$IC$K,
                           summary(modelF, .progressBar = F)$IC$K
                           ),
                      AICc=round(c(summary(modelA, .progressBar = F)$IC$AICc,
                           summary(modelB, .progressBar = F)$IC$AICc,
                           summary(modelC, .progressBar = F)$IC$AICc,
                           summary(modelD, .progressBar = F)$IC$AICc,
                           summary(modelE, .progressBar = F)$IC$AICc,
                           summary(modelF, .progressBar = F)$IC$AICc
                           ),2))%>%
  filter(p>0.05)%>%
  mutate(dAICc=AICc-min(AICc))%>%
  arrange(dAICc)%>%
  as_tibble()

aic.tab%>%
  write_csv(here::here("tables","aicc.csv"))

aic.tab%>%
  kable()
```

| model | description                              |     p |  K |   AICc | dAICc |
| :---- | :--------------------------------------- | ----: | -: | -----: | ----: |
| D     | green\>moose\>wolf\>caribou, ha          | 0.527 |  9 | 101.58 |  0.00 |
| B     | green\>moose\>wolf\>caribou, ha\>wolf    | 0.461 | 10 | 138.94 | 37.36 |
| E     | green\>moose\>wolf\>caribou, ha\>caribou | 0.450 | 10 | 139.50 | 37.92 |

\#\#bootstrap D-Separation analysis

``` r
####DSEP boot
mod.sel.compile.raw <- data.frame()
mod.sel.compile <- data.frame()
#len <- data.frame()
for(i in 1:1000){
  
  df.i <- df%>%sample_frac(1, replace=TRUE)
  while (length(unique(df.i$Name))<=4)
  {
    df.i <- df%>%sample_frac(1, replace=TRUE)
  }
#len <-rbind(len,data.frame(len=length(unique(df.i$Name))))


modelA <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density, df.i),
               lm(caribou.lambda.t ~ disturb.p, df.i))

modelB <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density + disturb.p, df.i),
               lm(caribou.lambda.t ~ WolfDensit, df.i))

modelC <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density, df.i),
               lm(caribou.lambda.t ~ Moose.Density+dEVI.t, df.i))%>%
  update(disturb.p ~ 1)

modelD <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density, df.i),
               lm(caribou.lambda.t ~ WolfDensit, df.i))%>%
  update(disturb.p ~ 1)

modelE <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density, df.i),
               lm(caribou.lambda.t ~ WolfDensit +disturb.p, df.i))

modelF <- psem(lm(Moose.Density ~ dEVI.t, df.i),
               lm(WolfDensit ~ Moose.Density, df.i),
               lm(caribou.lambda.t ~ dEVI.t, df.i))%>%
  update(disturb.p ~ 1)



##summarize
mod.sel <- data.frame(model=c("A","B","C","D","E","F"),
                      description=c(mA,mB,mC,mD,mE,mF),
                      p=round(c(summary(modelA, .progressBar = F)$Cstat$P.Value,
                           summary(modelB, .progressBar = F)$Cstat$P.Value,
                           summary(modelC, .progressBar = F)$Cstat$P.Value,
                           summary(modelD, .progressBar = F)$Cstat$P.Value,
                           summary(modelE, .progressBar = F)$Cstat$P.Value,
                           summary(modelF, .progressBar = F)$Cstat$P.Value
                           ),3),
                      K=c(summary(modelA, .progressBar = F)$IC$K,
                           summary(modelB, .progressBar = F)$IC$K,
                           summary(modelC, .progressBar = F)$IC$K,
                           summary(modelD, .progressBar = F)$IC$K,
                           summary(modelE, .progressBar = F)$IC$K,
                           summary(modelF, .progressBar = F)$IC$K
                           ),
                      AICc=round(c(summary(modelA, .progressBar = F)$IC$AICc,
                           summary(modelB, .progressBar = F)$IC$AICc,
                           summary(modelC, .progressBar = F)$IC$AICc,
                           summary(modelD, .progressBar = F)$IC$AICc,
                           summary(modelE, .progressBar = F)$IC$AICc,
                           summary(modelF, .progressBar = F)$IC$AICc
                           ),2))%>%
  mutate(dAICc=AICc-min(AICc),
         iter=i)%>%
  arrange(dAICc)

mod.sel.compile.raw <- rbind(mod.sel.compile.raw, mod.sel)

mod.sel.compile <- rbind(mod.sel.compile, mod.sel%>%filter(p>0.05)%>%mutate(dAICc=AICc-min(AICc)))
}
# mod.sel.compile.raw%>%
#   group_by(description)%>%
#   summarise(daic=mean(dAICc))
# 
# mod.sel.compile.raw%>%
#   group_by(description)%>%
#   summarise(p=mean(p))



##proportion of bootstrap samples where each model was top model (dAIC=0)
mod.sel.compile%>%
  filter(dAICc==0 & p>0.05)%>%
  count(description)%>%
  mutate(prop=((n/sum(n))*100)%>%round(1))%>%
  dplyr::select(-n)%>%
  arrange(-prop)%>%
  as_tibble()%>%
  kable()
```

| description                                            | prop |
| :----------------------------------------------------- | ---: |
| green\>moose\>wolf\>caribou, ha                        | 83.4 |
| green\>moose\>wolf\>caribou, ha\>wolf                  | 10.2 |
| green\>moose\>wolf\>caribou, ha\>caribou               |  2.2 |
| green\>moose\>wolf, green\>caribou, ha                 |  1.8 |
| green\>moose\>wolf, green\>caribou, moose\>caribou, ha |  1.4 |
| green\>moose\>wolf, ha\>caribou                        |  1.0 |

\#\#Plot paths

``` r
dag.dat <- data.frame()
mod.sel.compile.path <- mod.sel.compile%>%
  filter(dAICc==0 & p>0.05)
for(i in 1:nrow(mod.sel.compile.path)){
  a <- mod.sel.compile.path[i,]
  
  if(a$description %in% mA){
    b <- data.frame(from=c("green","moose","habitat alteration"),
                    to=c("moose","wolf","caribou"),
                    rel=c("+","+","-"),
                    i=i,
                    mod="A")
  }
  
  if(a$description %in% mB){
    b <- data.frame(from=c("green","moose","wolf","habitat alteration"),
                    to=c("moose","wolf","caribou","wolf"),
                    rel=c("+","+","-","+"),
                    i=i,
                    mod="B")
  }
  
  if(a$description %in% mC){
    b <- data.frame(from=c("green","moose","green","moose","habitat alteration"),
                    to=c("moose","wolf","caribou","caribou","habitat alteration"),
                    rel=c("+","+","+","-","-"),
                    i=i,
                    mod="C")
  }
  
  if(a$description %in% mD){
    b <- data.frame(from=c("green","moose","wolf","habitat alteration"),
                    to=c("moose","wolf","caribou","habitat alteration"),
                    rel=c("+","+","-","-"),
                    i=i,
                    mod="D")
  }
  
  if(a$description %in% mE){
    b <- data.frame(from=c("green","moose","wolf","habitat alteration"),
                    to=c("moose","wolf","caribou","caribou"),
                    rel=c("+","+","-","-"),
                    i=i,
                    mod="E")
  }
  
  if(a$description %in% mF){
    b <- data.frame(from=c("green","moose","green","habitat alteration"),
                    to=c("moose","wolf","caribou","habitat alteration"),
                    rel=c("+","+","+","-"),
                    i=i,
                    mod="F")
  }
  
 
  dag.dat <- rbind(dag.dat,b)
  
}

dag.dat<- dag.dat%>%
  mutate(from=case_when(from=="green"~"vegetation",
         TRUE~as.character(from)),
         to=case_when(to=="green"~"vegetation",
                      TRUE~as.character(to)),
         paste=paste0(from,to))
  


g <- graph_from_data_frame(dag.dat[1:nrow(mod.sel.compile),], vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat[1:nrow(mod.sel.compile),]$from, c("vegetation", "moose", "wolf", "caribou"))
to <- match(dag.dat[1:nrow(mod.sel.compile),]$to, c("vegetation", "moose", "wolf", "caribou"))

manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

set_graph_style(plot_margin = margin(1,1,1,1))
a <- ggraph(manual_layout) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.02, tension=0.9, 
                   position=position_jitter(width = 0.03, height = 0.03),
                   n=2) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()
a
```

![](README_files/figure-gfm/Plot%20paths-1.png)<!-- -->

``` r
###
m1a <- lm(Moose.Density~dEVI.t, data=df%>%mutate(dEVI.t=(dEVI.t-min(dEVI.t))/(max(dEVI.t)-min(dEVI.t)),
                                               Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                               WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                               caribou.lambda.t=(caribou.lambda.t-min(caribou.lambda.t))/(max(caribou.lambda.t)-min(caribou.lambda.t)))
)

summary(m1a)$r.squared
```

    ## [1] 0.4359054

``` r
summary(m1a)$coefficients[2,1]
```

    ## [1] 0.7367453

``` r
m1b <- lm(WolfDensit~Moose.Density, data=df%>%mutate(dEVI.t=(dEVI.t-min(dEVI.t))/(max(dEVI.t)-min(dEVI.t)),
                                               Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                               WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                               caribou.lambda.t=(caribou.lambda.t-min(caribou.lambda.t))/(max(caribou.lambda.t)-min(caribou.lambda.t)))
)
summary(m1b)$r.squared
```

    ## [1] 0.7687768

``` r
summary(m1b)$coefficients[2,1]
```

    ## [1] 0.7864172

``` r
m1c <- lm(caribou.lambda.t~WolfDensit, data=df%>%mutate(dEVI.t=(dEVI.t-min(dEVI.t))/(max(dEVI.t)-min(dEVI.t)),
                                                    Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                                    WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                                    caribou.lambda.t=(caribou.lambda.t-min(caribou.lambda.t))/(max(caribou.lambda.t)-min(caribou.lambda.t)))
)
summary(m1c)$r.squared
```

    ## [1] 0.7121397

``` r
summary(m1c)$coefficients[2,1]
```

    ## [1] -0.9406191

``` r
##from % models selected in bootstrap
dag.dat.top <-data.frame(from=c("vegetation","moose","wolf"),
           to=c("moose","wolf","caribou"),
           Direction=c("+","+","-"))
dag.dat.top
dag.dat.top$strength <- NA
dag.dat.top$strength[1] <- summary(m1a)$r.squared%>%round(2)
dag.dat.top$strength[2] <- summary(m1b)$r.squared%>%round(2)
dag.dat.top$strength[3] <- summary(m1c)$r.squared%>%round(2)
dag.dat.top$strength.lab <- NA
dag.dat.top$strength.lab[1] <- paste0(summary(m1a)$r.squared%>%round(2), " (", summary(m1a)$coefficients[2,1]%>%round(2), ")")
dag.dat.top$strength.lab[2] <- paste0(summary(m1b)$r.squared%>%round(2), " (", summary(m1b)$coefficients[2,1]%>%round(2), ")")
dag.dat.top$strength.lab[3] <- paste0(summary(m1c)$r.squared%>%round(2), " (", summary(m1c)$coefficients[2,1]%>%round(2), ")")



g <- graph_from_data_frame(dag.dat.top,
                           vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

b <- ggraph(manual_layout) + 
  geom_edge_link(aes(colour = Direction,label=strength.lab), 
                 width=2, 
                 angle_calc = 'along',
                 label_dodge = unit(3, 'mm'),
                 arrow = arrow(length = unit(3, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()
b
```

![](README_files/figure-gfm/Plot%20paths-2.png)<!-- -->

``` r
ggarrange(a,b,
          ncol = 2, nrow = 1,
          widths=c(1,1),
          labels="AUTO")
```

![](README_files/figure-gfm/Plot%20paths-3.png)<!-- -->

``` r
ggsave(here::here("plots","Fig3.png"), width=9, height=3.5, units="in")
```

\#\#DAGS

``` r
g <- graph_from_data_frame(data.frame(from=c("vegetation","moose", "habitat alteration"),
                                      to=c("moose","wolf","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

a <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()




g <- graph_from_data_frame(data.frame(from=c("vegetation","moose","wolf", "habitat alteration"),
                                      to=c("moose","wolf","caribou","wolf")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

b <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()



g <- graph_from_data_frame(data.frame(from=c("vegetation","moose","vegetation", "moose"),
                                      to=c("moose","wolf","caribou","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

c <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()


g <- graph_from_data_frame(data.frame(from=c("vegetation","moose","wolf"),
                                      to=c("moose","wolf","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

d <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()


g <- graph_from_data_frame(data.frame(from=c("vegetation","moose","wolf", "habitat alteration"),
                                      to=c("moose","wolf","caribou","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

e <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()


g <- graph_from_data_frame(data.frame(from=c("vegetation","moose","vegetation"),
                                      to=c("moose","wolf","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "nicely")
manual_layout$x <- c(0, 0.2, 0.8,1,0.9)
manual_layout$y <- c(0, 0.45, 0.55,1,0.2)

f <-ggraph(manual_layout) + 
  geom_edge_link(arrow = arrow(length = unit(1, 'mm')), 
                 end_cap = circle(5, 'mm')) + 
  coord_fixed()+
  geom_node_point(size = 5, color="grey")+
  geom_node_text(aes(filter=name%in%c("disturbance","vegetation",  "caribou", "habitat alteration"), label=name) ,hjust="inward",angle=0)+
  geom_node_text(aes(filter=name%in%c("moose"), label=name) ,hjust=1.2,vjust=0.05,angle=-20)+
  geom_node_text(aes(filter=name%in%c("wolf"), label=name) ,hjust=1.45,vjust=0.1,angle=-20)+
  theme_graph()

ggarrange(a,b,c,d,e,f,
          ncol = 3, nrow = 2,
          labels="AUTO")
```

![](README_files/figure-gfm/dags-1.png)<!-- -->

``` r
ggsave(here::here("plots","Fig2.png"), width=9, height=6, units="in")
```
