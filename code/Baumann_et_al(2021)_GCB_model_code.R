# Baumann et al (2021) GCB Mixed Effects Model Code #
#By: Justin Baumann

#####################################################

#load packages

library(tidyverse)
library(ggplot2)
library(lme4)
library(car) #For checking multicollinearity, using performance instead
library(performance) #for checking R2 or ICC, and comparing model performances
library(ggeffects)
#library(cowplot)
#library(egg)
#library(ggpubr)
library(patchwork)
library(ggrepel)
library(ggsci)
library(parameters) #get model params
library(effects)
library(broom)
library(devtools)
#devtools::install_github("m-clark/visibly")
library(visibly)
library(rgdal)
library(raster)
library(spdep)

######################################################

#data prep

#read in data
a=read.csv('coral_recov_master_2021.csv') 
a1<-a %>% mutate_if(is.character, as.factor)
nrow(a1)
head(a1)
str(a1)

########## Build recovery datasets#########

##building a recovery dataset
#remove Arabian Gulf sample (n=1)
levels(a1$region)
a2<-subset(a1, region != "Arabian Gulf")
droplevels(a2$region)

#remove disturbance= disease, dredging (n=1) & COTS, Storm (n=1)
levels(a2$disturbance)
a3<-subset(a2, disturbance != "Disease, Dredging")
droplevels(a3$disturbance)

a4<-subset(a3, disturbance != "COTS, Storm")
droplevels(a4$disturbance)


#remove NAs from data
recov<-a4 %>% drop_na(calculated.recovery.rate) %>% droplevels()
tail(recov)
nrow(recov) #182 rows of data for recovery!
str(recov)
levels(recov$disturbance)

#exploring the structure of the data
#HISTOGRAM OF RECOVERY RATES
hist(recov$calculated.recovery.rate) #positively skewed normal

ggplot(recov, aes(x=calculated.recovery.rate))+
  geom_histogram(binwidth=5)

#center (standardise) explanatory variable(s) (mean of zero=centering, sd=1 = scaling --doing both here)
recov$hii100km2<-scale(recov$hii100km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(recov$hii100km2)

recov$hii500km2<-scale(recov$hii500km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(recov$hii500km2)

recov$recovery_time2<-scale(recov$recovery_time, center=TRUE, scale=TRUE)
hist(recov$recovery_time2)

hist(recov$distance_to_shore_m)
recov$distance_to_shore_m2<-scale(recov$distance_to_shore_m, center=TRUE, scale=TRUE)
hist(recov$distance_to_shore_m2)

recov$dist_to_riv_m2<-scale(recov$dist_to_riv_m, center=TRUE, scale=TRUE)
hist(recov$dist_to_riv_m2)

hist(recov$grav_NC)
recov$grav_NC2<-scale(recov$grav_NC, center=TRUE, scale=TRUE)
hist(recov$grav_NC2)

hist(recov$cml_scr)
recov$cml_scr2<-scale(recov$cml_scr, center=TRUE, scale=TRUE)
hist(recov$cml_scr2)

hist(recov$gravity.Grav_tot)
recov$gravity.Grav_tot2<-scale(recov$gravity.Grav_tot, center=TRUE, scale=TRUE)
hist(recov$gravity.Grav_tot2)

hist(recov$travel_time.tt_pop)
recov$travel_time.tt_pop2<-scale(recov$travel_time.tt_pop, center=TRUE, scale=TRUE)
hist(recov$travel_time.tt_pop2)

head(recov)


#How many studies are in the recov dataset?
length(unique(recov$study)) #57 studies

#write.csv(recov, 'recov_2021.csv')


####EXPLORING DATA BEFORE MODEL TESTING####

recov<-read.csv('recov_2021.csv')


#basic lm
#how many unique lat/longs are there?
nrow(recov) #182
head(recov)

lm1<-lm(calculated.recovery.rate ~ hii100km2, data=recov)
summary(lm1)

lm2<-lm(calculated.recovery.rate ~ cml_scr2, data=recov)
summary(lm2)

lm3<-lm(calculated.recovery.rate ~ grav_NC2, data=recov)
summary(lm3)


plot1<-ggplot(recov, aes(x=hii100km, y=calculated.recovery.rate))+
  geom_point()+
  geom_smooth(method="lm")
plot1

plot2<-ggplot(recov, aes(x=cml_scr2, y=calculated.recovery.rate))+
  geom_point()+
  geom_smooth(method="lm")
plot2

plot3<-ggplot(recov, aes(x=grav_NC2, y=calculated.recovery.rate))+
  geom_point()+
  geom_smooth(method="lm")
plot3


plot1+plot2+plot3

#are assumptions met?
#residuals
plot(lm1, which=1) #gray line is flat, red should be nearly flat (mimicking gray)
#looks ok
plot(lm2, which=1) # looks ok
plot(lm3, which=1) # looks ok

#qq
plot(lm1, which=2) #points should be close to the line. They diverge at the ends a little
plot(lm2, which=2) #same as above
plot(lm3, which=2) #same as above

#check for observation independence (use categorical vars here)
#if data from within each category are more similar to each other than to data from different categories then they are correlated!

#region
boxplot(calculated.recovery.rate ~ region, data=recov) #maybe correlated? Possibly not.

#disturbance
boxplot(calculated.recovery.rate ~ disturbance, data= recov) #most likely correlated though sample sizes might be small

#plot w/ colors by category to see
color1<-ggplot(recov, aes(x=hii100km2, y=calculated.recovery.rate, color=region))+
  geom_point(size=2)+
  theme_bw()
color1
#regions vary by recovery rate and hii, so observations within are not independent

color2<-ggplot(recov, aes(x=hii100km2, y=calculated.recovery.rate, color=disturbance))+
  geom_point(size=2)+
  theme_bw()
color2
#disturbances vary by recovery rate AND hii100km, so observations within these are NOT INDEPENDENT 

#add fixed effects into the model
lm2.1<-lm(calculated.recovery.rate ~ hii100km2 + region + disturbance, data=recov)
summary(lm2.1)



######MODEL TESTING#######
###GLM MODELS FOR RECOVERY
#RUN ONCE! JUST READ .csv IN FOR GRAPHING!
# mixed effect model run using lmer()
#model testing (REML=FALSE for this step)
#REML=restricted (or residual) maximum likelihood
#fixed effects fo after the ~, random go after that var ~ fixed + (1| random)

#When VIFs are very high, collinearity is an issue. Remove VIFs > ~ 2 and retest. 

#######
#Building the "full" model. Can remove  variables one by one from here to find optimal model 
#######

#with all vars in
r1<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r1) 
performance::check_model(r1) 
check_outliers(r1)
summary(r1)
#distance to river has high VIF (18.52). So does region. So let's remove dist to river and try again

r2<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r2) #region still high, so something colinear w/ that... 
performance::check_model(r2) #this look pretty good
check_outliers(r2) #1 outlier detected 
summary(r2)

#cut gravity
r3<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r3) 
performance::check_model(r3) 
check_outliers(r3)  
summary(r3)

#cut distance from shore (and not gravity)
r4<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4) 
performance::check_model(r4) 
check_outliers(r4)  #no outliers
summary(r4)

#cut both distance from shore and gravity
r4.1<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1) #max VIF is region at 4.92
performance::check_model(r4.1) 
check_outliers(r4.1)  
summary(r4.1)

#cut travel time
r4.2<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.2) #max VIF region at 5.35
performance::check_model(r4.2) 
check_outliers(r4.2)  #none
summary(r4.2)

#cut travel time and gravity, ends up being just the old model +cml_scr
r4.3<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.3) #VERY low VIF (nothing about 2.06--we like this)
performance::check_model(r4.3) 
check_outliers(r4.3)  #none
summary(r4.3)

#original final model
r5<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r5) #all VIF below 2.38 :)
performance::check_model(r5) 
check_outliers(r5)  #1 outlier
summary(r5)

perf<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5, rank=TRUE)
perf

#testing some simpler models for comparisons sake

#only region, dist, hii, and cml_scr
r6<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+cml_scr2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r6) #all VIF are low
performance::check_model(r6) 
check_outliers(r6)  #3 outliers!
summary(r6)

#only region, disturbance, and hii
r7<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r7) #all VIF are low
performance::check_model(r7) 
check_outliers(r7)  #5 outliers!
summary(r7)

perf<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r6,r7, rank=TRUE)
perf

#old models (from model building)

r2.0<-lmer(calculated.recovery.rate ~ region*disturbance*hii100km2*MPA_during_recov + (1|study), data = recov, REML=FALSE)
summary(r2.0)
#if modeled slope is < modeled error, then the effect cannot be distinguished from zero!
#variance from random effects / total variance #Tells us how much var left over AFTER fixed effects is explained by random effects
#rank deficient

r2.1<-lmer(calculated.recovery.rate ~ region*disturbance*hii100km2*MPA_during_recov*recovery_time2 + (1|study), data = recov, REML=FALSE)# rank deficient
r2.2<-lmer(calculated.recovery.rate ~ region*disturbance*hii100km2*MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_model(r2.2)
check_collinearity(r2.2) #VIFs are very high, indicating collinearity is an issue. Remove VIFs > ~ 2 and retest. 
#can use the vif function in the cars package to see the VIFS as well
vif(r2.2)
#REMOVE THE INTERACTION TERMS 1 at a time and compare VIFs
r2.2.1<-lmer(calculated.recovery.rate ~ region*disturbance*hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.1) #still high
performance::check_model(r2.2.1)

r2.2.2<-lmer(calculated.recovery.rate ~ region*disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.2) #still high
performance::check_model(r2.2.2) 

r2.2.3<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.3) #This finally looks good!
performance::check_model(r2.2.3) #THIS LOOKS GOOD. Use this model moving forward.
summary(r2.2.3)

#add in both distance to nearest river and distance to shore
r2.2.4<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)
performance::check_model(r2.2.4) #multicollinearity is high
check_collinearity(r2.2.4) #region and dist to river are highly correlated. Remove dist to river

#just distance to shore added
r2.2.5<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)
performance::check_model(r2.2.5) #multicollinearity is fine
check_collinearity(r2.2.5)

#remove recovery time
r2.2.6<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)


#remove distance to shore
r2.2.7<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+ (1|study), data = recov, REML=FALSE)

#remove hii
r2.2.8<-lmer(calculated.recovery.rate ~ region+disturbance+MPA_during_recov+ (1|study), data = recov, REML=FALSE)

#remove MPA status
r2.2.9<-lmer(calculated.recovery.rate ~ region+disturbance+ (1|study), data = recov, REML=FALSE)

#remove region
r2.2.10<-lmer(calculated.recovery.rate ~ disturbance+ (1|study), data = recov, REML=FALSE)


##COMPARE model performance 
perf1<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r2.0,r2.1,r2.2,r2.2.1,r2.2.2,r2.2.3,r2.2.4,r2.2.5,r2.2.6,r2.2.7,r2.2.8,r2.2.9,r2.2.10, rank= TRUE)
perf1

perf2<-performance::compare_performance(r1,r2,r3,r4,r5,r4.1,r4.2,r4.3,rank=TRUE)
perf2

#compare_performance explained
#AIC, BIC are both criterion for model selection-> smaller is better
#R2_marginal is R2 consdiering the fixed factors only
#R2 condtiional takes both fixed and random factors into account
#ICC is Intraclass correlation coefficient - "the proportion of the variance explained by the grouping structure in the population" (Hox 2010)
#RMSE= root mean squared error (smaller is generally better)
#BF is bayes factor (not sure it is relevant)


#model r4.1 is best AIC and BIC. IT does not contain gravity or distance to shore
r4.1<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1) #highest VIF=4.92, suggests some multicolinearity
performance::check_model(r4.1) 
check_outliers(r4.1)  #non3
summary(r4.1)


#replace travel time w/ dist to shore 
r4.1.2<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+distance_to_shore_m2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1.2) #highest VIF = 2.06, so this is better than the above
performance::check_model(r4.1.2) 
check_outliers(r4.1.2)  #none
summary(r4.1.2)

perf3<-performance::compare_performance(r4.1,r4.1.2,rank=TRUE)
perf3

##4.1 has lowest AIC, so let's move forward with that one. 
#Differs from previous final model because cml_scr is added and dist to shore replaced w/ travel time

summary(r4.1)
0.883/(0.883+3.152) #21.88%

summary(r4.1.2)
3.978/(3.978+3.171) #55.64%

#add REML in
finr4.1<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=TRUE)
check_collinearity(finr4.1)#VIF max =3.85
performance::check_model(finr4.1) 
performance::check_outliers(finr4.1, method=c("cook","pareto"))  #none
summary(finr4.1)

2.148/(2.148+3.334) #39.18% of variance
#HII estimate = 0.80077 +/-0.41651
#cml_scr estimate = -0.02998 +/- 0.38993
#travel time estimate = 0.347 +/- 0.33291


compare_performance(r4.1,finr4.1, rank=TRUE)
#finr4.1 has lower AIC :)

###CHECK FOR SPATIAL AUTOCORRELATION###

#method 1: use the check_autocorrelation() feature of the performance package
check_autocorrelation(finr4.1)# OK: Residuals appear to be independent and not autocorrelated (p = 0.118).

#method 2:
#based on: https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
ggplot(recov, aes(x=long, y=lat, size=calculated.recovery.rate))+
  geom_point()+
  theme_bw()
#maybe size looks very close to the same throughout but maybe a pattern? let's investigate.

#We will check this using our "non-spatial" model (finr4.1) above

#first we need the residuals from the model
#option 1-> the resid() function generates just a list of residiuals
recov$resid<-resid(finr4.1)
resids<-as.data.frame(resid(finr4.1))
resids
#option 2: augment from the broom package makes a df of all model stuff
lmerresid<-broom::augment(finr4.1)
head(lmerresid)
nrow(lmerresid)#109
nrow(recov)#182

#problem: neither contains lat/long so a join or bind is needed. 
#luckily column X.2 on recov and column .rownames on the augmented data are the same! So we can bind by that is we name the columns the same thing
recov1<-rename(recov, .rownames = X.2)
head(recov1)

spatialdf<-merge(recov1, lmerresid, by='.rownames')
head(spatialdf)
#now we have resids and data together so we can carry on

#using tri2nb in the spdep package to find nearest neighbors of points
#https://rdrr.io/rforge/spdep/man/tri2nb.html


#make new column for longlat 
spatialdf$longlat<-paste(spatialdf$long, spatialdf$lat, sep="_")

#make a df with no duplicate lat/long
spdf2<-spatialdf %>%
  distinct(longlat, .keep_all=TRUE)
head(spdf2)

#make recov2 into a spatial object
WGScoor<-spdf2

coordinates(WGScoor)=~long+lat
proj4string(WGScoor)<-CRS("+proj=longlat +datum=WGS84")

#raster::shapefile(WGScoor, "recovshape.shp")

#WGScoor is a spatialpointsdataframe and was saved as a shape file

coords<-coordinates(WGScoor)

#find nearest neighbors
tri.nb<-tri2nb(coords, row.names=rownames(WGScoor))
tri.nb #THIS WORKED AND MADE nearest neighbors!
summary(tri.nb)

nb2listw(tri.nb) #this also worked!

head(spdf2)

vect=spdf2$.resid #vector of model residuals
vect

vect1=spdf2$calculated.recovery.rate.y #vector of response var
vect1

#MORANS test for spatial autocorrelation!
moran.test(vect, nb2listw(tri.nb))
# Moran I statistic standard deviate = -1.8338, p-value = 0.9667
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.122120186      -0.010989011       0.003672439 

moran.test(vect1, nb2listw(tri.nb))
# Moran I statistic standard deviate = 0.21191, p-value = 0.4161
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.001747535      -0.010989011       0.003612590 

#######################################################################################################
#######################################################################################################
#######################################################################################################

###RESISTANCE####

#read in data
a=read.csv('coral_recov_master_2021.csv') 
a1<-a %>% mutate_if(is.character, as.factor)
nrow(a1)
head(a1)
str(a1)

#BUILDING RESISTANCE DATASET#
#remove Arabian Gulf sample (n=1)
levels(a1$region)
a2<-subset(a1, region != "Arabian Gulf")
droplevels(a2$region)

#remove disturbance= disease, dredging (n=1) & COTS, Storm (n=1)
levels(a2$disturbance)
a3<-subset(a2, disturbance != "Disease, Dredging")
droplevels(a3$disturbance)

a4<-subset(a3, disturbance != "COTS, Storm")
droplevels(a4$disturbance)

#remove distrubance = Bleaching, Disease (n=1 for resistance only)
a5<-subset(a4, disturbance != "Bleaching, Disease")
droplevels(a5$disturbance)

#remove NAs from data
resist<-a5 %>% drop_na(resistance) %>% droplevels()
tail(resist)
nrow(resist) #184 rows of data
str(resist)
levels(resist$region)

#exploring the structure of the data
#HISTOGRAM OF RECOVERY RATES
hist(resist$resistance) #slight negative skew but approx normal

ggplot(resist, aes(x=resistance))+
  geom_histogram(binwidth=5)

#center (standardise) explanatory variable(s) (mean of zero=centering, sd=1 = scaling --doing both here)
resist$hii100km2<-scale(resist$hii100km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(resist$hii100km2)

resist$hii500km2<-scale(resist$hii500km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(resist$hii500km2)

resist$resistance_time_2<-scale(resist$resistance.time, center=TRUE, scale=TRUE)
hist(resist$resistance_time2)

hist(resist$distance_to_shore_m)
resist$distance_to_shore_m2<-scale(resist$distance_to_shore_m, center=TRUE, scale=TRUE)
hist(resist$distance_to_shore_m2)

recov$dist_to_riv_m2<-scale(recov$dist_to_riv_m, center=TRUE, scale=TRUE)
hist(recov$dist_to_riv_m2)

hist(resist$grav_NC)
resist$grav_NC2<-scale(resist$grav_NC, center=TRUE, scale=TRUE)
hist(resist$grav_NC2)

hist(resist$cml_scr)
resist$cml_scr2<-scale(resist$cml_scr, center=TRUE, scale=TRUE)
hist(resist$cml_scr2)

hist(resist$gravity.Grav_tot)
resist$gravity.Grav_tot2<-scale(resist$gravity.Grav_tot, center=TRUE, scale=TRUE)
hist(resist$gravity.Grav_tot2)

hist(resist$travel_time.tt_pop)
resist$travel_time.tt_pop2<-scale(resist$travel_time.tt_pop, center=TRUE, scale=TRUE)
hist(resist$travel_time.tt_pop2)

resist$dist_to_riv_m2<-scale(resist$dist_to_riv_m, center=TRUE, scale=TRUE)
hist(resist$dist_to_rivm2)


head(resist)

#How many studies are in the recov dataset?
length(unique(resist$study)) #59 studies

#write.csv(resist, 'resist_2021.csv')

######MODEL TESTING#######
###GLM MODELS FOR RECOVERY
#RUN ONCE! JUST READ .csv IN FOR GRAPHING!
# mixed effect model run using lmer()
#model testing (REML=FALSE for this step)
#REML=restricted (or residual) maximum likelihood
#fixed effects fo after the ~, random go after that var ~ fixed + (1| random)

#When VIFs are very high, collinearity is an issue. Remove VIFs > ~ 2 and retest. 

#######
#Building the "full" model. Can remove  variables one by one from here to find optimal model 
#######
resist<-read.csv('resist_2021.csv')

#with all vars in
s1<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+dist_to_riv_m2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s1) #dist to riv is HIGH correlation. Will need to remove
performance::check_model(s1) 
check_outliers(s1)#none
summary(s1)
#distance to river has high VIF (12.28). So does region. So let's remove dist to river and try again

s2<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s2) #Looks like we are good on this
performance::check_model(s2) #this look pretty good
check_outliers(s2) #none 
summary(s2)

#cut gravity since it built into cml_scr
s3<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s3) #good
performance::check_model(s3) 
check_outliers(s3)  #none
summary(s3)

#cut distance from shore (and not gravity)-- since travel time is used for same thing
s4<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4) #good
performance::check_model(s4) 
check_outliers(s4)  #no outliers
summary(s4)

#cut both distance from shore and gravity
s4.1<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1) #good
performance::check_model(s4.1) 
check_outliers(s4.1)  #none
summary(s4.1)

#cut travel time (see if this differs from using distance from shore)
s4.2<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.2) #good
performance::check_model(s4.2) 
check_outliers(s4.2)  #none
summary(s4.2)

#cut travel time and gravity, ends up being just the old model +cml_scr
s4.3<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+cml_scr2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.3) #good
performance::check_model(s4.3) 
check_outliers(s4.3)  #none
summary(s4.3)

#original final model
s5<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s5) #good)
performance::check_model(s5) 
check_outliers(s5)  #1 outlier
summary(s5)

#testing some simpler models for comparisons sake
#only region, dist, hii, and cml_scr
s6<-lmer(resistance ~ region+disturbance+hii100km2+cml_scr2 + (1|study), data = resist, REML=FALSE)
check_collinearity(s6) #all VIF are low
performance::check_model(s6) 
check_outliers(s6)  #3 outliers!
summary(s6)

#only region, disturbance, and hii
s7<-lmer(resistance ~ region+disturbance+hii100km2 + (1|study), data = resist, REML=FALSE)
check_collinearity(s7) #all VIF are low
performance::check_model(s7) 
check_outliers(s7)  #5 outliers!
summary(s7)

rperf<-compare_performance(s1,s2,s3,s4,s4.1,s4.2,s4.3,s5,s6,s7, rank=TRUE)
rperf

#old models (from model building)
s2.0<-lmer(resistance ~ region*disturbance*hii100km2*MPA_during_resist + (1|study), data = resist, REML=FALSE)
summary(s2.0)
#if modeled slope is < modeled error, then the effect cannot be distinguished from zero!
#variance from random effects / total variance #Tells us how much var left over AFTER fixed effects is explained by random effects
#rank deficient

s2.1<-lmer(resistance ~ region*disturbance*hii100km2*MPA_during_resist*resistance_time_2 + (1|study), data = resist, REML=FALSE)# rank deficient
s2.2<-lmer(resistance ~ region*disturbance*hii100km2*MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_model(s2.2)
check_collinearity(s2.2) #VIFs are very high, indicating collinearity is an issue. Remove VIFs > ~ 2 and retest. 
#can use the vif function in the cars package to see the VIFS as well
vif(s2.2)
#REMOVE THE INTERACTION TERMS 1 at a time and compare VIFs
s2.2.1<-lmer(resistance ~ region*disturbance*hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.1) #still high
performance::check_model(r2.2.1)

s2.2.2<-lmer(resistance ~ region*disturbance+hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.2) #still high
performance::check_model(r2.2.2) 

s2.2.3<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.3) #This finally looks good!
performance::check_model(r2.2.3) #THIS LOOKS GOOD. Use this model moving forward.
summary(r2.2.3)

#add in both distance to nearest river and distance to shore
s2.2.4<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+dist_to_riv_m2+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)
performance::check_model(r2.2.4) #multicollinearity is high
check_collinearity(r2.2.4) #region and dist to river are highly correlated. Remove dist to river

#just distance to shore added
s2.2.5<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)
performance::check_model(r2.2.5) #multicollinearity is fine
check_collinearity(r2.2.5)

#remove resistance time
s2.2.6<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)

#remove distance to shore
s2.2.7<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+ (1|study), data = resist, REML=FALSE)

#remove hii
s2.2.8<-lmer(resistance ~ region+disturbance+MPA_during_resist+ (1|study), data = resist, REML=FALSE)

#remove MPA status
s2.2.9<-lmer(resistance ~ region+disturbance+ (1|study), data = resist, REML=FALSE)

#remove region
s2.2.10<-lmer(resistance ~ disturbance+ (1|study), data = resist, REML=FALSE)


##COMPARE model performance 
perfs1<-performance::compare_performance(s1,s2,s3,s4,s4.1,s4.2,s4.3,s5,s2.0,s2.1,s2.2,s2.2.1,s2.2.2,s2.2.3,s2.2.4,s2.2.5,s2.2.6,s2.2.7,s2.2.8,s2.2.9,s2.2.10, rank= TRUE)
perfs1 #Model S3 (w/ distance to shore and travel time performs best but S2, S4.1, S4, S4.2, and S1 all VERY similar--essentially not different)
#WE WILL SELECT s4.1 as it matches what we did for recovery! s4.1 uses travel time in place of distance to shore as a proxy for remoteness and includes cml_scr from the WCS pre-print in addition to HII



#compare_performance explained
#AIC, BIC are both criterion for model selection-> smaller is better
#R2_marginal is R2 consdiering the fixed factors only
#R2 condtiional takes both fixed and random factors into account
#ICC is Intraclass correlation coefficient - "the proportion of the variance explained by the grouping structure in the population" (Hox 2010)
#RMSE= root mean squared error (smaller is generally better)
#BF is bayes factor (not sure it is relevant)


#does not contain gravity or distance to shore
s4.1<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1) #good
performance::check_model(s4.1) 
check_outliers(s4.1)  #none
summary(s4.1)

#replace travel time w/ dist to shore 
s4.1.2<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+distance_to_shore_m2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1.2) #good
performance::check_model(s4.1.2) 
check_outliers(s4.1.2)  #none
summary(s4.1.2)

perf3<-performance::compare_performance(s4.1,s4.1.2,rank=TRUE)
perf3

##S4.1 has lowest AIC, so let's move forward with that one. 
#Differs from previous final model because cml_scr is added and dist to shore replaced w/ travel time

summary(s4.1)
38.72/(38.72+75) #34.05%

summary(s4.1.2)
39.05/(39.05+84.87) #31.51%

#add REML in
fins4.1<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=TRUE)
check_collinearity(fins4.1)#good
performance::check_model(fins4.1) 
check_outliers(fins4.1)  #none
summary(fins4.1)

69.98/(69.98+77.00) #47.61% of variance
#HII estimate = -.57795 +/-1.90827
#cml_scr estimate = 1.70562 +/- 1.66645
#travel time estimate = -0.72707 +/- 1.27091

####SPATIAL AUTOCORRELATION FOR RESISTANCE
#Method 1: use the check_autocorrelation() function in the performance package
check_autocorrelation(fins4.1) #OK: Residuals appear to be independent and not autocorrelated (p = 0.610).

#method 2: 
ggplot(resist, aes(x=long, y=lat, size=resistance))+
  geom_point()+
  theme_bw()
#Doesn't look like much but we should investigate residuals

#We will check this using our "non-spatial" model (fins4.1) above

#first we need the residuals from the model
#option 1-> the resid() function generates just a list of residiuals
#resist$resid<-resid(fins4.1)
resids<-as.data.frame(resid(fins4.1))
resids
#option 2: augment from the broom package makes a df of all model stuff
lmerresidres<-broom::augment(fins4.1)
head(lmerresidres)
nrow(lmerresidres)#134
nrow(resist)#184

#problem: neither contains lat/long so a join or bind is needed. 
#luckily column X.2 on recov and column .rownames on the augmented data are the same! So we can bind by that is we name the columns the same thing
head(resist)
resist1<-rename(resist, .rownames = X.1)
head(resist1)

spatialdfres<-merge(resist1, lmerresidres, by='.rownames')
head(spatialdfres)
#now we have resids and data together so we can carry on

#using tri2nb in the spdep package to find nearest neighbors of points
#https://rdrr.io/rforge/spdep/man/tri2nb.html

#make new column for longlat 
spatialdfres$longlat<-paste(spatialdfres$long, spatialdfres$lat, sep="_")


#make a df with no duplicate lat/long
spdf2res<-spatialdfres %>%
  distinct(longlat, .keep_all=TRUE)
head(spdf2res)

#make recov2 into a spatial object
WGScoorres<-spdf2res
coordinates(WGScoorres)=~long+lat
proj4string(WGScoorres)<-CRS("+proj=longlat +datum=WGS84")

raster::shapefile(WGScoorres, "resistshape.shp")

#WGScoor is a spatialpointsdataframe and was saved as a shape file
coords<-coordinates(WGScoorres)

#find nearest neighbors
tri.nbres<-tri2nb(coords, row.names=rownames(WGScoorres))
tri.nbres #THIS WORKED AND MADE nearest neighbors!
summary(tri.nbres)

nb2listw(tri.nbres) #this also worked!

head(spdf2res)

vectr=spdf2res$.resid #vector of model residuals
vectr

vectr1=spdf2res$resistance.y #vector of response var
vectr1

#MORANS test for spatial autocorrelation!
moran.test(vectr, nb2listw(tri.nbres))
# Moran I statistic standard deviate = 0.60469, p-value = 0.2727
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.024899265      -0.009523810       0.003240714 

moran.test(vectr1, nb2listw(tri.nbres))
# Moran I statistic standard deviate = 0.47645, p-value = 0.3169
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.017786786      -0.009523810       0.003285649 
#############################################################
#############################################################

### IGR ######################################################

# Read in data
recov<-read.csv('recov_2021.csv')
head(recov)
### Instantaneous Growth Rate (IGR) from Ortiz et al (2018)
#calculated as: r = LN ((recovery coral cover + 5) / (post dist coral cover + 5))/recovery time

#histogram of IGR
ggplot(recov, aes(x=IGR))+
  geom_histogram(binwidth=1) #actually appears approx normal

#using same model structure we used for other vars:

######MODEL TESTING#######
###GLM MODELS FOR RECOVERY
#RUN ONCE! JUST READ .csv IN FOR GRAPHING!
# mixed effect model run using lmer()
#model testing (REML=FALSE for this step)
#REML=restricted (or residual) maximum likelihood
#fixed effects fo after the ~, random go after that var ~ fixed + (1| random)

#When VIFs are very high, collinearity is an issue. Remove VIFs > ~ 2 and retest. 

#######
#Building the "full" model. Can remove  variables one by one from here to find optimal model 
#######

#with all vars in
r1<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r1) 
performance::check_model(r1) 
check_outliers(r1) #1 outlier
summary(r1)
#dist to riv and region have HIGH VIF

r2<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r2) #looks good actually
performance::check_model(r2) #this look pretty good
check_outliers(r2) #1 outlier detected 
summary(r2)

#cut gravity
r3<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r3) #looks good
performance::check_model(r3) 
check_outliers(r3)  # 1 outlier
summary(r3)

#cut distance from shore (and not gravity)
r4<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4) #looks good
performance::check_model(r4) 
check_outliers(r4)  #no outliers
summary(r4)

#cut both distance from shore and gravity
r4.1<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1) #looks good
performance::check_model(r4.1) 
check_outliers(r4.1)  #none
summary(r4.1)

#cut travel time
r4.2<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.2) #looks good
performance::check_model(r4.2) 
check_outliers(r4.2)  #none
summary(r4.2)

#cut travel time and gravity, ends up being just the old model +cml_scr
r4.3<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.3) #VERY low VIF
performance::check_model(r4.3) 
check_outliers(r4.3)  #1 outlier
summary(r4.3)

#original final model
r5<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r5) #all VIF below 2.38 :)
performance::check_model(r5) 
check_outliers(r5)  #1 outlier
summary(r5)

compare_performance(r1,r2)

perf<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5, rank=TRUE)
perf
#model r4 wins, though lowest AIC and BIC are model r4.1

#testing some simpler models for comparisons sake

#only region, dist, hii, and cml_scr
r6<-lmer(IGR ~ region+disturbance+hii100km2+cml_scr2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r6) #all VIF are low
performance::check_model(r6) 
check_outliers(r6)  #1 outlier!
summary(r6)

#only region, disturbance, and hii
r7<-lmer(IGR ~ region+disturbance+hii100km2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r7) #all VIF are low
performance::check_model(r7) 
check_outliers(r7)  #1 outlier!
summary(r7)

perf<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r6,r7, rank=TRUE)
perf #still r4.1 and r4

#old models (from model building)

r2.0<-lmer(IGR ~ region*disturbance*hii100km2*MPA_during_recov + (1|study), data = recov, REML=FALSE)
summary(r2.0)
#if modeled slope is < modeled error, then the effect cannot be distinguished from zero!
#variance from random effects / total variance #Tells us how much var left over AFTER fixed effects is explained by random effects
#rank deficient

r2.1<-lmer(IGR ~ region*disturbance*hii100km2*MPA_during_recov*recovery_time2 + (1|study), data = recov, REML=FALSE)# rank deficient
r2.2<-lmer(IGR ~ region*disturbance*hii100km2*MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_model(r2.2) #singular fit
check_collinearity(r2.2) #VIFs are very high, indicating collinearity is an issue. Remove VIFs > ~ 2 and retest. 
#can use the vif function in the cars package to see the VIFS as well
vif(r2.2)
#REMOVE THE INTERACTION TERMS 1 at a time and compare VIFs
r2.2.1<-lmer(IGR ~ region*disturbance*hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.1) #still high
performance::check_model(r2.2.1) #rank deficient

r2.2.2<-lmer(IGR ~ region*disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.2) #still high
performance::check_model(r2.2.2) #rank deficient

r2.2.3<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = recov, REML=FALSE)
check_collinearity(r2.2.3) #This finally looks good!
performance::check_model(r2.2.3) #THIS LOOKS GOOD. Use this model moving forward.
summary(r2.2.3)

#add in both distance to nearest river and distance to shore
r2.2.4<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)
performance::check_model(r2.2.4) #multicollinearity is high
check_collinearity(r2.2.4) #region and dist to river are highly correlated. Remove dist to river

#just distance to shore added
r2.2.5<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)
performance::check_model(r2.2.5) #multicollinearity is fine
check_collinearity(r2.2.5)

#remove recovery time
r2.2.6<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+distance_to_shore_m2 + (1|study), data = recov, REML=FALSE)


#remove distance to shore
r2.2.7<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+ (1|study), data = recov, REML=FALSE)

#remove hii
r2.2.8<-lmer(IGR ~ region+disturbance+MPA_during_recov+ (1|study), data = recov, REML=FALSE)

#remove MPA status
r2.2.9<-lmer(IGR ~ region+disturbance+ (1|study), data = recov, REML=FALSE)

#remove region
r2.2.10<-lmer(IGR ~ disturbance+ (1|study), data = recov, REML=FALSE)


##COMPARE model performance 
perf1<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r2.0,r2.1,r2.2,r2.2.1,r2.2.2,r2.2.3,r2.2.4,r2.2.5,r2.2.6,r2.2.7,r2.2.8,r2.2.9,r2.2.10, rank= TRUE)
perf1

perf2<-performance::compare_performance(r1,r2,r3,r4,r5,r4.1,r4.2,r4.3,rank=TRUE)
perf2 #r4 and 4.1 still the best! Let's use 4.1 as it is the best for recov and resist

#compare_performance explained
#AIC, BIC are both criterion for model selection-> smaller is better
#R2_marginal is R2 consdiering the fixed factors only
#R2 condtiional takes both fixed and random factors into account
#ICC is Intraclass correlation coefficient - "the proportion of the variance explained by the grouping structure in the population" (Hox 2010)
#RMSE= root mean squared error (smaller is generally better)
#BF is bayes factor (not sure it is relevant)


#model r4.1 is best AIC and BIC. IT does not contain gravity or distance to shore
r4.1<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1) #all good 
performance::check_model(r4.1) 
check_outliers(r4.1)  #NONE
summary(r4.1)


#replace travel time w/ dist to shore 
r4.1.2<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+distance_to_shore_m2+pre.dist.cc + (1|study), data = recov, REML=FALSE)
check_collinearity(r4.1.2) #good
performance::check_model(r4.1.2) 
check_outliers(r4.1.2)  #1 outlier
summary(r4.1.2)

perf3<-performance::compare_performance(r4.1,r4.1.2,rank=TRUE)
perf3

##4.1 has lowest AIC, so let's move forward with that one. 
#Differs from previous final model because cml_scr is added and dist to shore replaced w/ travel time

summary(r4.1)
0.1348/(0.1348+0.1580) #46.04%

#add REML in
finr4.1igr<-lmer(IGR ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=TRUE)
check_collinearity(finr4.1igr)#ALL good
performance::check_model(finr4.1igr) 
performance::check_outliers(finr4.1igr, method=c("cook","pareto"))  #none
summary(finr4.1igr)

0.2619/(0.2619+0.1659) #61.22% of variance

compare_performance(r4.1,finr4.1igr, rank=TRUE)
#finr4.1 has lower AIC :)

##########################################################################################################

# Relative Recovery and Resistance ########################################################################

####EXPLORING DATA BEFORE MODEL TESTING####

relrecov<-read.csv('rel_recov_2021.csv')


#basic lm
#how many unique lat/longs are there?
nrow(relrecov) #151 (182 in the non-relative recov df--had to cut when we didn't have a pre-dist value to standardize to)
head(relrecov)

lm1<-lm(rel_rec_rate ~ hii100km2, data=relrecov)
summary(lm1)

lm2<-lm(rel_rec_rate ~ cml_scr2, data=relrecov)
summary(lm2)

lm3<-lm(rel_rec_rate ~ grav_NC2, data=relrecov)
summary(lm3)


plot1<-ggplot(relrecov, aes(x=hii100km, y=rel_rec_rate))+
  geom_point()+
  geom_smooth(method="lm")
plot1

plot2<-ggplot(relrecov, aes(x=cml_scr2, y=rel_rec_rate))+
  geom_point()+
  geom_smooth(method="lm")
plot2

plot3<-ggplot(relrecov, aes(x=grav_NC2, y=rel_rec_rate))+
  geom_point()+
  geom_smooth(method="lm")
plot3


plot1+plot2+plot3

#are assumptions met?
#residuals
plot(lm1, which=1) #gray line is flat, red should be nearly flat (mimicking gray)
#looks ok
plot(lm2, which=1) # looks ok
plot(lm3, which=1) # looks ok

#qq
plot(lm1, which=2) #points should be close to the line. They diverge at the ends a little
plot(lm2, which=2) #same as above
plot(lm3, which=2) #same as above

#check for observation independence (use categorical vars here)
#if data from within each category are more similar to each other than to data from different categories then they are correlated!

#region
boxplot(rel_rec_rate ~ region, data=relrecov) #maybe correlated? Possibly not.

#disturbance
boxplot(rel_rec_rate ~ disturbance, data= relrecov) #most likely correlated though sample sizes might be small

#plot w/ colors by category to see
color1<-ggplot(relrecov, aes(x=hii100km2, y=rel_rec_rate, color=region))+
  geom_point(size=2)+
  theme_bw()
color1
#regions vary by relrecovery rate and hii, so observations within are not independent

color2<-ggplot(relrecov, aes(x=hii100km2, y=rel_rec_rate, color=disturbance))+
  geom_point(size=2)+
  theme_bw()
color2
#disturbances vary by relrecovery rate AND hii100km, so observations within these are NOT INDEPENDENT 

#add fixed effects into the model
lm2.1<-lm(rel_rec_rate ~ hii100km2 + region + disturbance, data=relrecov)
summary(lm2.1)



######MODEL TESTING#######
###GLM MODELS FOR relrecovERY
#RUN ONCE! JUST READ .csv IN FOR GRAPHING!
# mixed effect model run using lmer()
#model testing (REML=FALSE for this step)
#REML=restricted (or residual) maximum likelihood
#fixed effects fo after the ~, random go after that var ~ fixed + (1| random)

#When VIFs are very high, collinearity is an issue. Remove VIFs > ~ 2 and retest. 

#######
#Building the "full" model. Can remove  variables one by one from here to find optimal model 
#######

#with all vars in
r1<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r1) 
performance::check_model(r1) 
check_outliers(r1) #1 outlier
summary(r1)
#Singular fit -- so need to simplify. VIFs for region and dist to riv are very high. CUT dist to riv

r2<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r2) #region still high, so something colinear w/ that... 
performance::check_model(r2) #this look pretty good
check_outliers(r2) #1 outlier detected 
summary(r2)

#cut gravity
r3<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r3) #still moderate multicollinearity (VIF for region = 7.05, disturbance=5.34)
performance::check_model(r3) 
check_outliers(r3)#1 outlier   
summary(r3)

#cut distance from shore (and not gravity)
r4<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4) #still moderate VIF for region and dist
performance::check_model(r4) 
check_outliers(r4)  #1 outlier
summary(r4)

#cut both distance from shore and gravity
r4.1<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4.1) #still moderate VIFS (region =6.79, dist=5.16)
performance::check_model(r4.1) 
check_outliers(r4.1)  #no outliers
summary(r4.1)

#cut travel time
r4.2<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4.2) #Still moderate Vifs for region and dist
performance::check_model(r4.2) 
check_outliers(r4.2)  #none
summary(r4.2)

#cut travel time and gravity, ends up being just the old model +cml_scr
r4.3<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+cml_scr2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4.3) #Vifs are less bad BUT SINGULAR FIT. So cannot use
performance::check_model(r4.3) 
check_outliers(r4.3)  #none
summary(r4.3)

#original final model
r5<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r5) #VIFS are low but singular fit!
performance::check_model(r5) 
check_outliers(r5)  #none
summary(r5)

perf<-performance::compare_performance(r4.1,r4.2,r4.3,r5)
perf

#testing some simpler models for comparisons sake

#only region, dist, hii, and cml_scr
r6<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+cml_scr2 + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r6) #all VIF are low but SINGULAR FIT
performance::check_model(r6) 
check_outliers(r6)  #none
summary(r6)

#only region, disturbance, and hii
r7<-lmer(rel_rec_rate ~ region+disturbance+hii100km2 + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r7) #all VIF are low but SINGULAR FIT
performance::check_model(r7) 
check_outliers(r7)  #1 outlier
summary(r7)

perf<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r6,r7, rank=TRUE)
perf

#old models (from model building)

r2.0<-lmer(rel_rec_rate ~ region*disturbance*hii100km2*MPA_during_recov + (1|study), data = relrecov, REML=FALSE)
summary(r2.0)
#if modeled slope is < modeled error, then the effect cannot be distinguished from zero!
#variance from random effects / total variance #Tells us how much var left over AFTER fixed effects is explained by random effects
#rank deficient

r2.1<-lmer(rel_rec_rate ~ region*disturbance*hii100km2*MPA_during_recov*recovery_time2 + (1|study), data = relrecov, REML=FALSE)# rank deficient
r2.2<-lmer(rel_rec_rate ~ region*disturbance*hii100km2*MPA_during_recov+recovery_time2 + (1|study), data = relrecov, REML=FALSE)
check_model(r2.2)
check_collinearity(r2.2) #VIFs are very high, indicating collinearity is an issue. Remove VIFs > ~ 2 and retest. 
#can use the vif function in the cars package to see the VIFS as well
vif(r2.2)
#REMOVE THE INTERACTION TERMS 1 at a time and compare VIFs
r2.2.1<-lmer(rel_rec_rate ~ region*disturbance*hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r2.2.1) #still high
performance::check_model(r2.2.1)

r2.2.2<-lmer(rel_rec_rate ~ region*disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r2.2.2) #still high
performance::check_model(r2.2.2) 

r2.2.3<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2 + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r2.2.3) #This finally looks good!
performance::check_model(r2.2.3) #THIS LOOKS GOOD. Use this model moving forward.
summary(r2.2.3)

#add in both distance to nearest river and distance to shore
r2.2.4<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+dist_to_riv_m2+distance_to_shore_m2 + (1|study), data = relrecov, REML=FALSE)
performance::check_model(r2.2.4) #multicollinearity is high
check_collinearity(r2.2.4) #region and dist to river are highly correlated. Remove dist to river

#just distance to shore added
r2.2.5<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+distance_to_shore_m2 + (1|study), data = relrecov, REML=FALSE)
performance::check_model(r2.2.5) #multicollinearity is fine
check_collinearity(r2.2.5)

#remove relrecovery time
r2.2.6<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+distance_to_shore_m2 + (1|study), data = relrecov, REML=FALSE)


#remove distance to shore
r2.2.7<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+ (1|study), data = relrecov, REML=FALSE)

#remove hii
r2.2.8<-lmer(rel_rec_rate ~ region+disturbance+MPA_during_recov+ (1|study), data = relrecov, REML=FALSE)

#remove MPA status
r2.2.9<-lmer(rel_rec_rate ~ region+disturbance+ (1|study), data = relrecov, REML=FALSE)

#remove region
r2.2.10<-lmer(rel_rec_rate ~ disturbance+ (1|study), data = relrecov, REML=FALSE)


##COMPARE model performance 
perf1<-performance::compare_performance(r1,r2,r3,r4,r4.1,r4.2,r4.3,r5,r2.0,r2.1,r2.2,r2.2.1,r2.2.2,r2.2.3,r2.2.4,r2.2.5,r2.2.6,r2.2.7,r2.2.8,r2.2.9,r2.2.10, rank= TRUE)
perf1

perf2<-performance::compare_performance(r1,r2,r3,r4,r5,r4.1,r4.2,r4.3,rank=TRUE)
perf2

#compare_performance explained
#AIC, BIC are both criterion for model selection-> smaller is better
#R2_marginal is R2 consdiering the fixed factors only
#R2 condtiional takes both fixed and random factors into account
#ICC is Intraclass correlation coefficient - "the proportion of the variance explained by the grouping structure in the population" (Hox 2010)
#RMSE= root mean squared error (smaller is generally better)
#BF is bayes factor (not sure it is relevant)


#model r4.1 is best AIC and BIC. IT does not contain gravity or distance to shore
r4.1<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4.1) #highest VIF=6.79, suggest some multicollinearity!
performance::check_model(r4.1) #these look ok though except for multicollinearity
check_outliers(r4.1)  #none
summary(r4.1) #AIC = 776.1, BIC=824.6


#replace travel time w/ dist to shore 
r4.1.2<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+distance_to_shore_m2+pre.dist.cc + (1|study), data = relrecov, REML=FALSE)
check_collinearity(r4.1.2) #VIFs a little lower but SINGULAR FIT. 
performance::check_model(r4.1.2) 
check_outliers(r4.1.2)  #none
summary(r4.1.2)

perf3<-performance::compare_performance(r4.1,r4.1.2,rank=TRUE)
perf3

##4.1 has lowest AIC and is not a singular fit. In spite of collinearity concerns, lets use it use it. 
#Differs from previous final model because cml_scr is added and dist to shore replaced w/ travel time

summary(r4.1)
3.932/(3.932+48.775) #7.46%



#add REML in
finr4.1<-lmer(rel_rec_rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = relrecov, REML=TRUE)
check_collinearity(finr4.1)#VIF max =5.29 (region). We can live with that I think
performance::check_model(finr4.1) #looks good for (note 1 VIF > 5--region)
performance::check_outliers(finr4.1, method=c("cook","pareto"))  #none
summary(finr4.1)

11.60/(11.60+53.68) #17.76961% of variance
#HII estimate = 2.456 +/-1.34970
#cml_scr estimate = 0.22960 +/- 1.40278
#travel time estimate = 1.12155 +/- 1.10504


compare_performance(r4.1,finr4.1, rank=TRUE)
#finr4.1 has lower AIC :)

###CHECK FOR SPATIAL AUTOCORRELATION###

#method 1: use the check_autocorrelation() feature of the performance package
check_autocorrelation(finr4.1)# OK: Residuals appear to be independent and not autocorrelated (p = 0.812).

#method 2:
#based on: https://datascienceplus.com/spatial-regression-in-r-part-1-spamm-vs-glmmtmb/
ggplot(relrecov, aes(x=long, y=lat, size=rel_rec_rate))+
  geom_point()+
  theme_bw()
#maybe size looks very close to the same throughout but maybe a pattern? let's investigate.

#We will check this using our "non-spatial" model (finr4.1) above

#first we need the residuals from the model
#option 1-> the resid() function generates just a list of residiuals
relrecov$resid<-resid(finr4.1)
resids<-as.data.frame(resid(finr4.1))
resids

#option 2: augment from the broom package makes a df of all model stuff
#due to updates and changes we now need the broom.mixed package
library(broom.mixed)

lmerresid<-broom.mixed::augment(finr4.1)
head(lmerresid)
nrow(lmerresid)#109
nrow(relrecov)#151

#problem: neither contains lat/long so a join or bind is needed. 
#luckily column X.2 on relrecov and column .rownames on the augmented data are the same! So we can bind by that is we name the columns the same thing
relrecov1<-rename(relrecov, .rownames = X.2)
head(relrecov1)

spatialdf<-merge(relrecov1, lmerresid, by='.rownames')
head(spatialdf)
#now we have resids and data together so we can carry on

#using tri2nb in the spdep package to find nearest neighbors of points
#https://rdrr.io/rforge/spdep/man/tri2nb.html


#make new column for longlat 
spatialdf$longlat<-paste(spatialdf$long, spatialdf$lat, sep="_")

#make a df with no duplicate lat/long
spdf2<-spatialdf %>%
  distinct(longlat, .keep_all=TRUE)
head(spdf2)

#make relrecov2 into a spatial object
WGScoor<-spdf2

coordinates(WGScoor)=~long+lat
proj4string(WGScoor)<-CRS("+proj=longlat +datum=WGS84")

#raster::shapefile(WGScoor, "relrecovshape.shp")

#WGScoor is a spatialpointsdataframe and was saved as a shape file

coords<-coordinates(WGScoor)

#find nearest neighbors
tri.nb<-tri2nb(coords, row.names=rownames(WGScoor))
tri.nb #THIS WORKED AND MADE nearest neighbors!
summary(tri.nb)

nb2listw(tri.nb) #this also worked!

head(spdf2)

vect=spdf2$.resid #vector of model residuals
vect

vect1=spdf2$rel_rec_rate.y #vector of response var
vect1

#MORANS test for spatial autocorrelation!
moran.test(vect, nb2listw(tri.nb))
# Moran I statistic standard deviate = -0.39564, p-value = 0.6538
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# -0.037737131      -0.012500000       0.004068959

moran.test(vect1, nb2listw(tri.nb))
# Moran I statistic standard deviate = 1.2743, p-value = 0.1013
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.067806433      -0.012500000       0.003971795 

#######################################################################################################
#######################################################################################################
#######################################################################################################

###RESISTANCE####

#read in data
a=read.csv('coral_relrecov_master_2021.csv') 
a1<-a %>% mutate_if(is.character, as.factor)
nrow(a1)
head(a1)
str(a1)

#BUILDING RESISTANCE DATASET#
#remove Arabian Gulf sample (n=1)
levels(a1$region)
a2<-subset(a1, region != "Arabian Gulf")
droplevels(a2$region)

#remove disturbance= disease, dredging (n=1) & COTS, Storm (n=1)
levels(a2$disturbance)
a3<-subset(a2, disturbance != "Disease, Dredging")
droplevels(a3$disturbance)

a4<-subset(a3, disturbance != "COTS, Storm")
droplevels(a4$disturbance)

#remove distrubance = Bleaching, Disease (n=1 for resistance only)
a5<-subset(a4, disturbance != "Bleaching, Disease")
droplevels(a5$disturbance)

#remove NAs from data
resist<-a5 %>% drop_na(resistance) %>% droplevels()
tail(resist)
nrow(resist) #184 rows of data
str(resist)
levels(resist$region)

#exploring the structure of the data
#HISTOGRAM OF relrecovERY RATES
hist(resist$resistance) #slight negative skew but approx normal

ggplot(resist, aes(x=resistance))+
  geom_histogram(binwidth=5)

#center (standardise) explanatory variable(s) (mean of zero=centering, sd=1 = scaling --doing both here)
resist$hii100km2<-scale(resist$hii100km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(resist$hii100km2)

resist$hii500km2<-scale(resist$hii500km, center=TRUE, scale=TRUE) #only works for numeric variables
hist(resist$hii500km2)

resist$resistance_time_2<-scale(resist$resistance.time, center=TRUE, scale=TRUE)
hist(resist$resistance_time2)

hist(resist$distance_to_shore_m)
resist$distance_to_shore_m2<-scale(resist$distance_to_shore_m, center=TRUE, scale=TRUE)
hist(resist$distance_to_shore_m2)

relrecov$dist_to_riv_m2<-scale(relrecov$dist_to_riv_m, center=TRUE, scale=TRUE)
hist(relrecov$dist_to_riv_m2)

hist(resist$grav_NC)
resist$grav_NC2<-scale(resist$grav_NC, center=TRUE, scale=TRUE)
hist(resist$grav_NC2)

hist(resist$cml_scr)
resist$cml_scr2<-scale(resist$cml_scr, center=TRUE, scale=TRUE)
hist(resist$cml_scr2)

hist(resist$gravity.Grav_tot)
resist$gravity.Grav_tot2<-scale(resist$gravity.Grav_tot, center=TRUE, scale=TRUE)
hist(resist$gravity.Grav_tot2)

hist(resist$travel_time.tt_pop)
resist$travel_time.tt_pop2<-scale(resist$travel_time.tt_pop, center=TRUE, scale=TRUE)
hist(resist$travel_time.tt_pop2)

resist$dist_to_riv_m2<-scale(resist$dist_to_riv_m, center=TRUE, scale=TRUE)
hist(resist$dist_to_rivm2)


head(resist)

#How many studies are in the relrecov dataset?
length(unique(resist$study)) #59 studies

#write.csv(resist, 'rel_resist_2021.csv')

######MODEL TESTING#######
###GLM MODELS FOR relrecovERY
#RUN ONCE! JUST READ .csv IN FOR GRAPHING!
# mixed effect model run using lmer()
#model testing (REML=FALSE for this step)
#REML=restricted (or residual) maximum likelihood
#fixed effects fo after the ~, random go after that var ~ fixed + (1| random)

#When VIFs are very high, collinearity is an issue. Remove VIFs > ~ 2 and retest. 

#######
#Building the "full" model. Can remove  variables one by one from here to find optimal model 
#######
resist<-read.csv('rel_resist_2021.csv')

#with all vars in
s1<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+dist_to_riv_m2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s1) #dist to riv is HIGH correlation. Will need to remove
performance::check_model(s1) 
check_outliers(s1)#none
summary(s1)
#distance to river has high VIF (13.25). So does region. So let's remove dist to river and try again

s2<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s2) #Looks like we are good on this
performance::check_model(s2) #this look pretty good
check_outliers(s2) #none 
summary(s2)

#cut gravity since it built into cml_scr
s3<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s3) #good
performance::check_model(s3) 
check_outliers(s3)  #1 outlier
summary(s3)

#cut distance from shore (and not gravity)-- since travel time is used for same thing
s4<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+gravity.Grav_tot2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4) #good
performance::check_model(s4) 
check_outliers(s4)  #no outliers
summary(s4)

#cut both distance from shore and gravity
s4.1<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1) #good
performance::check_model(s4.1) 
check_outliers(s4.1)  #none
summary(s4.1) #AIC = 1275.1, BIC=1324.4

#cut travel time (see if this differs from using distance from shore)
s4.2<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+gravity.Grav_tot2+cml_scr2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.2) #good
performance::check_model(s4.2) 
check_outliers(s4.2)  #1 outlier
summary(s4.2) #AIC=1277.1, BIC=1329.3

#cut travel time and gravity, ends up being just the old model +cml_scr
s4.3<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+cml_scr2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.3) #good
performance::check_model(s4.3) 
check_outliers(s4.3)  #none
summary(s4.3)

#original final model
s5<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s5) #good)
performance::check_model(s5) 
check_outliers(s5)  #none
summary(s5)

#testing some simpler models for comparisons sake
#only region, dist, hii, and cml_scr
s6<-lmer(rel_res ~ region+disturbance+hii100km2+cml_scr2 + (1|study), data = resist, REML=FALSE)
check_collinearity(s6) #all VIF are low
performance::check_model(s6) 
check_outliers(s6)  #3 outliers!
summary(s6)

#only region, disturbance, and hii
s7<-lmer(rel_res ~ region+disturbance+hii100km2 + (1|study), data = resist, REML=FALSE)
check_collinearity(s7) #all VIF are low
performance::check_model(s7) 
check_outliers(s7)  #5 outliers!
summary(s7)

rperf<-compare_performance(s1,s2,s3,s4,s4.1,s4.2,s4.3,s5,s6,s7, rank=TRUE)
rperf

#old models (from model building)
s2.0<-lmer(rel_res ~ region*disturbance*hii100km2*MPA_during_resist + (1|study), data = resist, REML=FALSE)
summary(s2.0)
#if modeled slope is < modeled error, then the effect cannot be distinguished from zero!
#variance from random effects / total variance #Tells us how much var left over AFTER fixed effects is explained by random effects
#rank deficient

s2.1<-lmer(rel_res ~ region*disturbance*hii100km2*MPA_during_resist*resistance_time_2 + (1|study), data = resist, REML=FALSE)# rank deficient
s2.2<-lmer(rel_res ~ region*disturbance*hii100km2*MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_model(s2.2)
check_collinearity(s2.2) #VIFs are very high, indicating collinearity is an issue. Remove VIFs > ~ 2 and retest. 
#can use the vif function in the cars package to see the VIFS as well
vif(s2.2)
#REMOVE THE INTERACTION TERMS 1 at a time and compare VIFs
s2.2.1<-lmer(rel_res ~ region*disturbance*hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.1) #still high
performance::check_model(r2.2.1)

s2.2.2<-lmer(rel_res ~ region*disturbance+hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.2) #still high
performance::check_model(r2.2.2) 

s2.2.3<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2 + (1|study), data = resist, REML=FALSE)
check_collinearity(r2.2.3) #This finally looks good!
performance::check_model(r2.2.3) #THIS LOOKS GOOD. Use this model moving forward.
summary(r2.2.3)

#add in both distance to nearest river and distance to shore
s2.2.4<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+dist_to_riv_m2+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)
performance::check_model(r2.2.4) #multicollinearity is high
check_collinearity(r2.2.4) #region and dist to river are highly correlated. Remove dist to river

#just distance to shore added
s2.2.5<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)
performance::check_model(r2.2.5) #multicollinearity is fine
check_collinearity(r2.2.5)

#remove resistance time
s2.2.6<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+distance_to_shore_m2 + (1|study), data = resist, REML=FALSE)

#remove distance to shore
s2.2.7<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+ (1|study), data = resist, REML=FALSE)

#remove hii
s2.2.8<-lmer(rel_res ~ region+disturbance+MPA_during_resist+ (1|study), data = resist, REML=FALSE)

#remove MPA status
s2.2.9<-lmer(rel_res ~ region+disturbance+ (1|study), data = resist, REML=FALSE)

#remove region
s2.2.10<-lmer(rel_res ~ disturbance+ (1|study), data = resist, REML=FALSE)


##COMPARE model performance 
perfs1<-performance::compare_performance(s1,s2,s3,s4,s4.1,s4.2,s4.3,s5,s2.0,s2.1,s2.2,s2.2.1,s2.2.2,s2.2.3,s2.2.4,s2.2.5,s2.2.6,s2.2.7,s2.2.8,s2.2.9,s2.2.10, rank= TRUE)
perfs1 #Model S3 (w/ distance to shore and travel time performs best but S2, S4.1, S4, S4.2, and S1 all VERY similar--essentially not different)
#WE WILL SELECT s4.1 as it matches what we did for relrecovery! s4.1 uses travel time in place of distance to shore as a proxy for remoteness and includes cml_scr from the WCS pre-print in addition to HII



#compare_performance explained
#AIC, BIC are both criterion for model selection-> smaller is better
#R2_marginal is R2 consdiering the fixed factors only
#R2 condtiional takes both fixed and random factors into account
#ICC is Intraclass correlation coefficient - "the proportion of the variance explained by the grouping structure in the population" (Hox 2010)
#RMSE= root mean squared error (smaller is generally better)
#BF is bayes factor (not sure it is relevant)


#does not contain gravity or distance to shore
s4.1<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1) #good
performance::check_model(s4.1) 
check_outliers(s4.1)  #none
summary(s4.1) #AIC=1275.1, BIC=1324.4

#replace travel time w/ dist to shore 
s4.1.2<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+distance_to_shore_m2+pre.dist.cc + (1|study), data = resist, REML=FALSE)
check_collinearity(s4.1.2) #good
performance::check_model(s4.1.2) 
check_outliers(s4.1.2)  #none
summary(s4.1.2) #AIC = 1722.3, BIC=1776.3

perf3<-performance::compare_performance(s4.1,s4.1.2,rank=TRUE)
perf3

##S4.1 has lowest AIC, so let's move forward with that one. 
#Differs from previous final model because cml_scr is added and dist to shore replaced w/ travel time

summary(s4.1)
191.9/(191.9+498.1) #27.81%

#add REML in
fins4.1<-lmer(rel_res ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=TRUE)
check_collinearity(fins4.1)#good
performance::check_model(fins4.1) 
check_outliers(fins4.1)  #none
summary(fins4.1)

400.3/(400.3+502.9) #44.32% of variance
#HII estimate = -5.1338 +/-4.7084
#cml_scr estimate = 8.9376 +/- 4.1876
#travel time estimate = -0.1736 +/-3.1860

####SPATIAL AUTOCORRELATION FOR RESISTANCE
#Method 1: use the check_autocorrelation() function in the performance package
check_autocorrelation(fins4.1) #OK: Residuals appear to be independent and not autocorrelated (p = 0.344).

#method 2: 
ggplot(resist, aes(x=long, y=lat, size=resistance))+
  geom_point()+
  theme_bw()
#Doesn't look like much but we should investigate residuals

#We will check this using our "non-spatial" model (fins4.1) above

#first we need the residuals from the model

#option 1-> the resid() function generates just a list of residiuals
#resist$resid<-resid(fins4.1)
resids<-as.data.frame(resid(fins4.1))
resids

#option 2: augment from the broom package makes a df of all model stuff
lmerresidres<-broom.mixed::augment(fins4.1)
head(lmerresidres)
nrow(lmerresidres)#134
nrow(resist)#184

#problem: neither contains lat/long so a join or bind is needed. 
#luckily column X.2 on relrecov and column .rownames on the augmented data are the same! So we can bind by that is we name the columns the same thing
head(resist)

resist1<-rename(resist, .rownames = X.1)
head(resist1)

spatialdfres<-merge(resist1, lmerresidres, by='.rownames')
head(spatialdfres)
#now we have resids and data together so we can carry on

#using tri2nb in the spdep package to find nearest neighbors of points
#https://rdrr.io/rforge/spdep/man/tri2nb.html

#make new column for longlat 
spatialdfres$longlat<-paste(spatialdfres$long, spatialdfres$lat, sep="_")


#make a df with no duplicate lat/long
spdf2res<-spatialdfres %>%
  distinct(longlat, .keep_all=TRUE)
head(spdf2res)

#make relrecov2 into a spatial object
WGScoorres<-spdf2res
coordinates(WGScoorres)=~long+lat
proj4string(WGScoorres)<-CRS("+proj=longlat +datum=WGS84")

raster::shapefile(WGScoorres, "resistshape.shp")

#WGScoor is a spatialpointsdataframe and was saved as a shape file
coords<-coordinates(WGScoorres)

#find nearest neighbors
tri.nbres<-tri2nb(coords, row.names=rownames(WGScoorres))
tri.nbres #THIS WORKED AND MADE nearest neighbors!
summary(tri.nbres)

nb2listw(tri.nbres) #this also worked!

head(spdf2res)

vectr=spdf2res$.resid #vector of model residuals
vectr

vectr1=spdf2res$resistance.y #vector of response var
vectr1

#MORANS test for spatial autocorrelation!
moran.test(vectr, nb2listw(tri.nbres))
# Moran I statistic standard deviate = 0.47746, p-value = 0.3165
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.017704076      -0.009523810       0.003252007 

moran.test(vectr1, nb2listw(tri.nbres))
#FAILS- is not a numeric vector






