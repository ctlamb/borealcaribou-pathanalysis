Caribou Path Analysis
================
Clayton T. Lamb
06 May, 2020

Load Data, Functions and Cleanup Data
-------------------------------------

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
library(usethis)
library(tidyverse)

##data
df <- read.csv(here::here("data", "final.csv"))%>%
  filter(Name!="Tweedsmuir") ##remove Tweedsmuir- non-boreal

##transform to instantaneous rate of growth (r)
df$caribou.lambda <- log(df$lambda)

##set seed for any bootstrapping
set.seed(2019)
```

Plot raw data
-------------

``` r
a <-ggplot(df, aes(x=disturb.p, y=LAI))+
  geom_point()+
  theme_bw()

b <-ggplot(df, aes(x=LAI, y=Moose.Density))+
  geom_point()+
  theme_bw()+
  xlab("Vegetation index")+
  ylab(expression(Moose~(ind./100~km^2)))

c <-ggplot(df, aes(x=Moose.Density, y=WolfDensit))+
  geom_point()+
  theme_bw()+
  xlab(expression(Moose~(ind./100~km^2)))+
  ylab(expression(Wolf~(ind./1000~km^2)))

d <-ggplot(df, aes(x=WolfDensit, y=caribou.lambda))+
  geom_point()+
  theme_bw()+
  #geom_vline(xintercept=6.5)+
  geom_hline(yintercept=0, linetype="dashed")+
  xlab(expression(Wolf~(ind./1000~km^2)))+
  ylab("Caribou inst. pop. growth (r)")

e <-ggplot(df, aes(x=WolfDensit, y=lambda))+
  geom_point()+
  theme_bw()+
  #geom_vline(xintercept=6.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  xlab(expression(Wolf~(ind./1000~km^2)))+
  ylab("Caribou pop. growth")


ggarrange(a,b,c,d,nrow=2,ncol=2, labels ="AUTO")
```

![](README_files/figure-markdown_github/plot%20raw%20data-1.png)

``` r
ggsave(here::here("plots","univar.png"), width=7, height=2.5, units="in")
ggarrange(b,c,e,nrow=1,ncol=3, labels ="AUTO")
```

![](README_files/figure-markdown_github/plot%20raw%20data-2.png)

``` r
ggsave(here::here("plots","univar2.png"), width=7, height=2.5, units="in")


f <- ggplot(df, aes(x=WolfDensit, y=survival))+
  geom_point()+
  theme_bw()+
  xlab(expression(wolf~(n/1000~km^2)))

g <- ggplot(df, aes(x=WolfDensit, y=reproduction))+
  geom_point()+
  theme_bw()+
  xlab(expression(wolf~(n/1000~km^2)))



ggarrange(f,g,
          ncol = 2, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-markdown_github/plot%20raw%20data-3.png)

``` r
ggsave(here::here("plots","vitalrate_wolf.png"), width=6, height=2.7, units="in")
```

Find intersections
------------------

``` r
###Wolf-Caribou intersection

##1.9 wolves/1000 sq.km generally reached when moose are greater than 2.9/100 sq.km
predict(lm(WolfDensit~caribou.lambda+I(caribou.lambda^2)+I(caribou.lambda^3), data=df), newdata=data.frame(caribou.lambda=0))
```

    ##        1 
    ## 1.904718

``` r
##plot
df$predicted <- predict(lm(WolfDensit~caribou.lambda+I(caribou.lambda^2)+I(caribou.lambda^3), data=df), newdata=df)
ggplot(df)+
  geom_point(aes(x=WolfDensit, y=caribou.lambda))+
  geom_line(aes(y=caribou.lambda, x=predicted))+
  theme_bw()+
  geom_vline(xintercept=1.9, linetype="dashed", col="red")+
  xlab(expression(wolf~(n/1000~km^2)))+
  ylab("caribou pop. growth (r)")
```

![](README_files/figure-markdown_github/Find%20intersections-1.png)

``` r
##population declines when > 1.9 wolves/1000 sq.km


##Moose-Wolf intersection

##1.9 wolves/1000 sq.km generally reached when moose are greater than 2.9/100 sq.km
predict(lm(Moose.Density~WolfDensit+I(WolfDensit^2), data=df), newdata=data.frame(WolfDensit=1.9))
```

    ##        1 
    ## 2.951979

``` r
##plot
df$predicted.moose <- predict(lm(Moose.Density~WolfDensit+I(WolfDensit^2), data=df), newdata=df)
ggplot(df)+
  geom_point(aes(y=WolfDensit, x=Moose.Density))+
  geom_line(aes(y=WolfDensit, x=predicted.moose))+
    geom_vline(xintercept=2.95, linetype="dashed", col="red")+
  theme_bw()+
  xlab(expression(moose~(n/100~km^2)))+
  ylab(expression(wolf~(n/1000~km^2)))
```

![](README_files/figure-markdown_github/Find%20intersections-2.png)

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

transformations to linear
-------------------------

``` r
###transform caribou lambda
df$caribou.lambda <--(exp(-10*df$caribou.lambda))

##make sure the rest remain linear 
a <- ggplot(df, aes(x=WolfDensit, y=caribou.lambda))+
  geom_point()+
  theme_bw()

b <- ggplot(df, aes(x=Moose.Density, y=caribou.lambda))+
  geom_point()+
  theme_bw()

c <- ggplot(df, aes(x=LAI, y=caribou.lambda))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()

ggarrange(a,b,c,
          ncol = 3, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-markdown_github/Transform-1.png)

``` r
###transform LAI
df$LAI <-exp(0.005*df$LAI)

##make sure the rest remain linear
d <- ggplot(df, aes(y=Moose.Density, x=LAI))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()

e <- ggplot(df, aes(y=WolfDensit, x=LAI))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()

f <- ggplot(df, aes(y=caribou.lambda, x=LAI))+
  geom_point()+
    xlab("Vegetation index")+
  theme_bw()

ggarrange(d,e,f,
          ncol = 3, nrow = 1,
          labels="AUTO")
```

![](README_files/figure-markdown_github/Transform-2.png)

Correlation marix
-----------------

``` r
M <- cor(df%>%
           select(disturb.p,LAI, Moose.Density, WolfDensit,caribou.lambda)%>%
           rename(`Habitat alteration`=disturb.p,
                  `Vegetation index`=LAI,
                  `Moose density`=Moose.Density,
                  `Wolf density`=WolfDensit,
                  `Caribou pop. growth`=caribou.lambda), use="complete.obs")

corrplot(M, method = "number", type = "upper", order = "hclust")
```

![](README_files/figure-markdown_github/corr-1.png)

D-Separation analysis
---------------------

``` r
###How to count parameters:
##see Shipley 2013 https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-0976.1

#A
##green>moose>wolf, ha>caribou
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{disturb.p}
### (caribou.lambda,Moose.Density)|{disturb.p, LAI}
### (caribou.lambda,WolfDensit)|{disturb.p, Moose.Density}
### (WolfDensit, disturb.p) | {Moose.Density}
### (Moose.Density, disturb.p) | {LAI}
KA <- 8

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df),
  lm(caribou.lambda~LAI + disturb.p, data=df),
  lm(caribou.lambda~Moose.Density + disturb.p + LAI, data=df),
  lm(caribou.lambda~WolfDensit + disturb.p + Moose.Density, data=df),
  lm(WolfDensit~disturb.p + Moose.Density, data=df),
  lm(Moose.Density~disturb.p + LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pA <- p.value
fcA <- Fisher.c+((2*KA)*(nrow(df)/(nrow(df)-KA-1)))
mA <- "green>moose>wolf, ha>caribou"


#B
##green>moose>wolf>caribou, ha>wolf
##independence statements
### (WolfDensit,LAI)|{Moose.Density, disturb.p}
### (caribou.lambda,LAI)|{WolfDensit}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI}
### (caribou.lambda,disturb.p)|{WolfDensit}
### (Moose.Density, disturb.p) | {LAI}
KB <- 9

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density + disturb.p, data=df),
  lm(caribou.lambda~LAI + WolfDensit, data=df),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI, data=df),
  lm(caribou.lambda~disturb.p + WolfDensit, data=df),
  lm(Moose.Density~disturb.p + LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pB <- p.value
fcB <- Fisher.c+((2*KB)*(nrow(df)/(nrow(df)-KB-1)))
mB <- "green>moose>wolf>caribou, ha>wolf"


##C
##green>moose>wolf, green>caribou, moose>caribou, ha
##independence statements
### (WolfDensit,caribou.lambda)|{Moose.Density, LAI}
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,disturb.p) | {Moose.Density, LAI}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}
KC <- 9

indepence.test <- list(
  lm(WolfDensit~ caribou.lambda + Moose.Density + LAI, data=df),
  lm(WolfDensit~ LAI + Moose.Density, data=df),
  lm(caribou.lambda~ disturb.p + Moose.Density + LAI, data=df),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df),
  lm(Moose.Density~ disturb.p +  LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pC <- p.value
fcC <- Fisher.c+((2*KC)*(nrow(df)/(nrow(df)-KC-1)))
mC <- "green>moose>wolf, green>caribou, moose>caribou, ha"



#D
##test path 1 green>moose>wolf>caribou, ha nowhere
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{WolfDensit}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI}
### (caribou.lambda,disturb.p) | {WolfDensit}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}
KD <- 8
  
indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df),
  lm(caribou.lambda~LAI + WolfDensit, data=df),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI, data=df),
  lm(caribou.lambda~ disturb.p + WolfDensit, data=df),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df),
  lm(Moose.Density~ disturb.p +  LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pD <- p.value
fcD <- Fisher.c+((2*KD)*(nrow(df)/(nrow(df)-KD-1)))
mD <- "green>moose>wolf>caribou, ha"



#E
##green>moose>wolf>caribou, ha>caribou
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{WolfDensit, disturb.p}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI, disturb.p}
### (WolfDensit, disturb.p) | {Moose.Density}
### (Moose.Density, disturb.p) | {LAI}
KE <- 9

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df),
  lm(caribou.lambda~LAI + WolfDensit + disturb.p, data=df),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI + disturb.p, data=df),
  lm(WolfDensit~disturb.p + Moose.Density, data=df),
  lm(Moose.Density~disturb.p + LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pE <- p.value
fcE <- Fisher.c+((2*KE)*(nrow(df)/(nrow(df)-KE-1)))
mE <- "green>moose>wolf>caribou, ha>caribou"


##F
##green>moose>wolf, green>caribou
##independence statements
### (caribou.lambda,Moose.Density)|{LAI}
### (caribou.lambda,Wolf.Densit)|{Moose.Density, LAI}
### (caribou.lambda,disturb.p) | {LAI}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}

KF <- 8

indepence.test <- list(
  lm(caribou.lambda~ Moose.Density + LAI, data=df),
  lm(caribou.lambda~ WolfDensit + Moose.Density + LAI, data=df),
  lm(caribou.lambda~ disturb.p + LAI, data=df),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df),
  lm(Moose.Density~ disturb.p +  LAI, data=df)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pF <- p.value
fcF <- Fisher.c+((2*KF)*(nrow(df)/(nrow(df)-KF-1)))
mF <- "green>moose>wolf, green>caribou, ha"

##summarize
data.frame(model=c("A","B","C","D","E","F"),
                      description=c(mA,mB,mC,mD,mE,mF),
                      K=c(KA,KB,KC,KD,KE,KF),
                      p=round(c(pA,pB,pC,pD,pE,pF),3),
                      AICc=round(c(fcA,fcB,fcC,fcD,fcE,fcF),2))%>%
  mutate(dAICc=AICc-min(AICc))%>%
  arrange(dAICc)%>%
  write_csv(here::here("tables", "aic.csv"))

###Are the two other paths B,E, that were excluded statistically important?
#lm(caribou.lambda~ disturb.p + WolfDensit, data=df)%>%summary() 
##no, wolf density remains significantly negative (p=0.008), disturb.p has no effect (p=0.82)

#lm(WolfDensit~ disturb.p + Moose.Density, data=df)%>%summary() 
##no, moose density remains significantly positive (p=0.003), disturb.p has no effect (p=0.76)

read_csv(here::here("tables", "aic.csv"))%>%print()
```

bootstrap D-Separation analysis
-------------------------------

``` r
####DSEP boot
mod.sel.compile.raw <- data.frame()
len <- data.frame()
for(i in 1:1000){
  
  df.i <- df%>%sample_frac(1, replace=TRUE)
  while (length(unique(df.i$Name))<=2)
  {
    df.i <- df%>%sample_frac(1, replace=TRUE)
  }
len <-rbind(len,data.frame(len=length(unique(df.i$Name))))


#A
##green>moose>wolf, ha>caribou
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{disturb.p}
### (caribou.lambda,Moose.Density)|{disturb.p, LAI}
### (caribou.lambda,WolfDensit)|{disturb.p, Moose.Density}
### (WolfDensit, disturb.p) | {Moose.Density}
### (Moose.Density, disturb.p) | {LAI}
KA <- 8

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df.i),
  lm(caribou.lambda~LAI + disturb.p, data=df.i),
  lm(caribou.lambda~Moose.Density + disturb.p + LAI, data=df.i),
  lm(caribou.lambda~WolfDensit + disturb.p + Moose.Density, data=df.i),
  lm(WolfDensit~disturb.p + Moose.Density, data=df.i),
  lm(Moose.Density~disturb.p + LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pA <- p.value
fcA <- Fisher.c+((2*KA)*(nrow(df.i)/(nrow(df.i)-KA-1)))
mA <- "green>moose>wolf, ha>caribou"


#B
##green>moose>wolf>caribou, ha>wolf
##independence statements
### (WolfDensit,LAI)|{Moose.Density, disturb.p}
### (caribou.lambda,LAI)|{WolfDensit}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI}
### (caribou.lambda,disturb.p)|{WolfDensit}
### (Moose.Density, disturb.p) | {LAI}
KB <- 9

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density + disturb.p, data=df.i),
  lm(caribou.lambda~LAI + WolfDensit, data=df.i),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI, data=df.i),
  lm(caribou.lambda~disturb.p + WolfDensit, data=df.i),
  lm(Moose.Density~disturb.p + LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pB <- p.value
fcB <- Fisher.c+((2*KB)*(nrow(df.i)/(nrow(df.i)-KB-1)))
mB <- "green>moose>wolf>caribou, ha>wolf"


##C
##green>moose>wolf, green>caribou, moose>caribou, ha
##independence statements
### (WolfDensit,caribou.lambda)|{Moose.Density, LAI}
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,disturb.p) | {Moose.Density, LAI}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}
KC <- 9

indepence.test <- list(
  lm(WolfDensit~ caribou.lambda + Moose.Density + LAI, data=df.i),
  lm(WolfDensit~ LAI + Moose.Density, data=df.i),
  lm(caribou.lambda~ disturb.p + Moose.Density + LAI, data=df.i),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df.i),
  lm(Moose.Density~ disturb.p +  LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pC <- p.value
fcC <- Fisher.c+((2*KC)*(nrow(df.i)/(nrow(df.i)-KC-1)))
mC <- "green>moose>wolf, green>caribou, moose>caribou, ha"



#D
##test path 1 green>moose>wolf>caribou, ha nowhere
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{WolfDensit}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI}
### (caribou.lambda,disturb.p) | {WolfDensit}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}
KD <- 8

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df.i),
  lm(caribou.lambda~LAI + WolfDensit, data=df.i),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI, data=df.i),
  lm(caribou.lambda~ disturb.p + WolfDensit, data=df.i),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df.i),
  lm(Moose.Density~ disturb.p +  LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pD <- p.value
fcD <- Fisher.c+((2*KD)*(nrow(df.i)/(nrow(df.i)-KD-1)))
mD <- "green>moose>wolf>caribou, ha"



#E
##green>moose>wolf>caribou, ha>caribou
##independence statements
### (WolfDensit,LAI)|{Moose.Density}
### (caribou.lambda,LAI)|{WolfDensit, disturb.p}
### (caribou.lambda,Moose.Density)|{WolfDensit, LAI, disturb.p}
### (WolfDensit, disturb.p) | {Moose.Density}
### (Moose.Density, disturb.p) | {LAI}
KE <- 9

indepence.test <- list(
  lm(WolfDensit~ LAI + Moose.Density, data=df.i),
  lm(caribou.lambda~LAI + WolfDensit + disturb.p, data=df.i),
  lm(caribou.lambda~Moose.Density + WolfDensit + LAI + disturb.p, data=df.i),
  lm(WolfDensit~disturb.p + Moose.Density, data=df.i),
  lm(Moose.Density~disturb.p + LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pE <- p.value
fcE <- Fisher.c+((2*KE)*(nrow(df.i)/(nrow(df.i)-KE-1)))
mE <- "green>moose>wolf>caribou, ha>caribou"


##F
##green>moose>wolf, green>caribou
##independence statements
### (caribou.lambda,Moose.Density)|{LAI}
### (caribou.lambda,Wolf.Densit)|{Moose.Density, LAI}
### (caribou.lambda,disturb.p) | {LAI}
### (WolfDensit,disturb.p) | {Moose.Density}
### (Moose.Density,disturb.p) | {LAI}

KF <- 8

indepence.test <- list(
  lm(caribou.lambda~ Moose.Density + LAI, data=df.i),
  lm(caribou.lambda~ WolfDensit + Moose.Density + LAI, data=df.i),
  lm(caribou.lambda~ disturb.p + LAI, data=df.i),
  lm(WolfDensit~ disturb.p + Moose.Density, data=df.i),
  lm(Moose.Density~ disturb.p +  LAI, data=df.i)
)
P <- sapply(indepence.test, function(x) coef(summary(x))[2,4]) 
Fisher.c <- -2*sum(log(P))
d <- 2*length(P)
p.value <- 1 - pchisq(Fisher.c, d)
pF <- p.value
fcF <- Fisher.c+((2*KF)*(nrow(df.i)/(nrow(df.i)-KF-1)))
mF <- "green>moose>wolf, green>caribou, ha"


##summarize
mod.sel <- data.frame(model=c("A","B","C","D","E","F"),
                      description=c(mA,mB,mC,mD,mE,mF),
                      K=c(KA,KB,KC,KD,KE,KF),
                      p=round(c(pA,pB,pC,pD,pE,pF),3),
                      AICc=round(c(fcA,fcB,fcC,fcD,fcE,fcF),2))%>%
  mutate(dAICc=AICc-min(AICc))%>%
  arrange(dAICc)
mod.sel.compile.raw <- rbind(mod.sel.compile.raw, mod.sel)

}
# mod.sel.compile.raw%>%
#   group_by(description)%>%
#   summarise(daic=mean(dAICc))
# 
# mod.sel.compile.raw%>%
#   group_by(description)%>%
#   summarise(p=mean(p))

mod.sel.compile <- mod.sel.compile.raw%>%
  filter(dAICc==0 & p>.05)

##proportion of bootstrap samples where each model was top model (dAIC=0)
mod.sel.compile.raw%>%
  filter(dAICc==0 & p>0.05)%>%
  count(description)%>%
  mutate(prop=((n/sum(n))*100)%>%round(1))%>%
  print()
```

Plot paths
----------

``` r
dag.dat <- data.frame()
for(i in 1:nrow(mod.sel.compile)){
  a <- mod.sel.compile[i,]
  
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
  


g <- graph_from_data_frame(dag.dat[1:1000,], vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat[1:1000,]$from, c("vegetation", "moose", "wolf", "caribou"))
to <- match(dag.dat[1:1000,]$to, c("vegetation", "moose", "wolf", "caribou"))
manual_layout <- create_layout(graph = g,
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
set_graph_style(plot_margin = margin(1,1,1,1))
a <- ggraph(manual_layout) + 
  geom_conn_bundle(data = get_con(from = from, to = to), alpha = 0.03, tension=0.9, 
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

![](README_files/figure-markdown_github/Plot%20paths-1.png)

``` r
###
m1a <- lm(Moose.Density~LAI, data=df%>%mutate(LAI=(LAI-min(LAI))/(max(LAI)-min(LAI)),
                                               Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                               WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                               caribou.lambda=(caribou.lambda-min(caribou.lambda))/(max(caribou.lambda)-min(caribou.lambda)))
)

summary(m1a)$r.squared
```

    ## [1] 0.5115998

``` r
m1b <- lm(WolfDensit~Moose.Density, data=df%>%mutate(LAI=(LAI-min(LAI))/(max(LAI)-min(LAI)),
                                               Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                               WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                               caribou.lambda=(caribou.lambda-min(caribou.lambda))/(max(caribou.lambda)-min(caribou.lambda)))
)
summary(m1b)$r.squared
```

    ## [1] 0.7801225

``` r
m1c <- lm(caribou.lambda~WolfDensit, data=df%>%mutate(LAI=(LAI-min(LAI))/(max(LAI)-min(LAI)),
                                                    Moose.Density=(Moose.Density-min(Moose.Density))/(max(Moose.Density)-min(Moose.Density)),
                                                    WolfDensit=(WolfDensit-min(WolfDensit))/(max(WolfDensit)-min(WolfDensit)),
                                                    caribou.lambda=(caribou.lambda-min(caribou.lambda))/(max(caribou.lambda)-min(caribou.lambda)))
)
summary(m1c)$r.squared
```

    ## [1] 0.8036942

``` r
##from % models selected in bootstrap
dag.dat.top <-data.frame(from=c("vegetation","moose","wolf"),
           to=c("moose","wolf","caribou"),
           direction=c("+","+","-"))
dag.dat.top
```

    ##         from      to direction
    ## 1 vegetation   moose         +
    ## 2      moose    wolf         +
    ## 3       wolf caribou         -

``` r
dag.dat.top$strength <- NA
dag.dat.top$strength[1] <- summary(m1a)$r.squared
dag.dat.top$strength[2] <- summary(m1b)$r.squared
dag.dat.top$strength[3] <- summary(m1c)$r.squared


g <- graph_from_data_frame(dag.dat.top, vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat.top$from, c("vegetation", "moose", "wolf", "caribou"))
to <- match(dag.dat.top$to, c("vegetation", "moose", "wolf", "caribou"))
manual_layout <- create_layout(graph = g,
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
b <- ggraph(manual_layout) + 
  geom_edge_link(aes(colour = direction,label=round(strength,2)), 
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

![](README_files/figure-markdown_github/Plot%20paths-2.png)

``` r
ggarrange(a,b,
          ncol = 2, nrow = 1,
          widths=c(1,1),
          labels="AUTO")
```

![](README_files/figure-markdown_github/Plot%20paths-3.png)

``` r
ggsave(here::here("plots","Fig3.png"), width=9, height=3.5, units="in")
```

DAGS
----

``` r
g <- graph_from_data_frame(data.frame(from=c("vegetation","moose", "habitat alteration"),
                                      to=c("moose","wolf","caribou")), vertices = c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
from <- match(dag.dat$from, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
to <- match(dag.dat$to, c("vegetation", "moose", "wolf", "caribou", "habitat alteration"))
manual_layout <- create_layout(graph = g,
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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
                               layout = "manual", node.positions = data.frame(x = c(0, 0.2, 0.8,1,0.9),
                                                                              y = c(0, 0.45, 0.55,1,0.2)))
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

![](README_files/figure-markdown_github/dags-1.png)

``` r
ggsave(here::here("plots","Fig2.png"), width=9, height=6, units="in")
```
