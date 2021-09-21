### Baumann et al (2021)-GCB - Figure and graphic code###
#By: Justin Baumann

#########################################################

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
library(ggridges)
library(formattable)
library(RColorBrewer)

### Table 1 ##########################################
#NEW TABLE 1

tab1<-read.csv('table1.csv')

customRed = "#ff7f7f"
customBlue = "#006699" 
customRed2 = "#cc0033"



tab1<- tab1%>%
  arrange(desc(Recovery.Rate))
head(tab1)

formattable(tab1,
            align =c('l','c','c','c','c','c','c','c','c'),
            list(`Location` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")),
                 'Recovery.Rate'= color_tile('firebrick','steelblue3'),
                 'IGR' = color_tile('firebrick','steelblue3'),
                 'Resistance' = color_tile('firebrick','steelblue3'),
                 'Human.Influence.Index' = color_bar(customRed),
                 'WCS.Cumulative.Score' = color_bar(customRed)
            ))


### Figure 1 #########################################

# Fig 1 is a conceptual diagram. Made in inkscape. 



### FIGURE 2: Histograms ##############################

##Build the Recovery Histogram
recov$region<-factor(recov$region, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

#overall recovery histo
recovr<-ggplot(recov, aes(x=calculated.recovery.rate))+
  geom_density(fill='gray')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "Recovery Rate", y='Frequency')+
  theme(axis.text = element_text(size=10))+
  theme(axis.title = element_text(size=12))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  theme(panel.border=element_blank())

recovr

#by region

#calculate outliers for recovery rate and generate a new column (outlier)
recouts<-recov%>%
  group_by(region) %>%
  mutate(outlier = ifelse(is_outlier(calculated.recovery.rate), location, as.numeric(NA)))
recouts
#write.csv(recouts, 'recouts.csv')



#calculate IQRs for each region (Q1 and Q3)
recouts1<-recov%>%
  group_by(region)%>%
  summarise_at(vars(calculated.recovery.rate),
               list(IQR=IQR, Q1=~quantile(.,probs=0.25), median=median, Q3=~quantile(.,probs=0.75)))
recouts1
#add outlier lines (1.5*Q1 or Q3)
recouts1$outlierval=1.5*recouts1$IQR
recouts1$highout=recouts1$Q3+recouts1$outlierval
recouts1$lowout=recouts1$Q1-recouts1$outlierval
head(recouts1)

head(recouts)
#merge iqr info to recouts

recouts2<-merge(recouts, recouts1, by='region')
head(recouts2)

#histo for recovery by region (can add outliers with commented out code)
rechist<-ggplot(recouts, aes(x=calculated.recovery.rate, y=region, fill=region))+
  geom_density_ridges(jittered_points=TRUE, position=position_points_jitter(width=0.1, height=0), point_shape='|',alpha=0.8, quantile_lines=TRUE)+
  geom_segment(data=recouts1, aes(x=lowout, xend=lowout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  geom_segment(data=recouts1, aes(x=highout, xend=highout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "Recovery Rate", y='Region')+
  scale_color_jama()+
  scale_fill_jama()
rechist

#+
#add place names to the outliers!
# geom_text(aes(color=region),na.rm=TRUE, size=2.5, vjust=2, position= position_jitter(width=0, height=0.2))



##Build the Resistance Histrogram
resist$region<-factor(resist$region, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

#overall resistance
resovr<-ggplot(resist, aes(x=resistance))+
  geom_density(fill='gray')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "Resistance", y="Frequency")+
  theme(axis.text = element_text(size=10))+
  theme(axis.title = element_text(size=12))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  theme(panel.border = element_blank())
resovr

###RES outliers
#calculate outliers for recovery rate and generate a new column (outlier)
resouts<-resist%>%
  group_by(region) %>%
  mutate(outlier = ifelse(is_outlier(resistance), location, as.numeric(NA)))
resouts
#write.csv(resouts, 'resouts.csv')


#calculate IQRs for each region (Q1 and Q3)
resouts1<-resist%>%
  group_by(region)%>%
  summarise_at(vars(resistance),
               list(IQR=IQR, Q1=~quantile(.,probs=0.25), median=median, Q3=~quantile(.,probs=0.75)))
resouts1
#add outlier lines (1.5*Q1 or Q3)
resouts1$outlierval=1.5*resouts1$IQR
resouts1$highout=resouts1$Q3+resouts1$outlierval
resouts1$lowout=resouts1$Q1-resouts1$outlierval
head(resouts1)

head(resouts)
#merge iqr info to recouts

resouts2<-merge(resouts, resouts1, by='region')
head(resouts2)

##by region
reshist<-ggplot(resist, aes(x=resistance, y=region, fill=region))+
  geom_density_ridges(jittered_points=TRUE, position=position_points_jitter(width=0.1, height=0), point_shape='|',alpha=0.8, quantile_lines=TRUE)+
  geom_segment(data=resouts1, aes(x=lowout, xend=lowout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  geom_segment(data=resouts1, aes(x=highout, xend=highout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "Resistance", y='Region')+
  scale_color_jama()+
  scale_fill_jama()

reshist

##Build global maps of Recov and Resist

###plotting them together in a panel plot

#histograms
recovr
resovr

#Inset of overall histos
layout<- c(
  patchwork::area(t=2, l=1, b=9, r=8),
  patchwork::area(t = 1, l = 6, b = 2, r=8),
  patchwork::area(t= 2, l=9, b=9, r=17),
  patchwork::area(t = 1, l = 15, b = 2, r =17))

#make final Fig 2
Fig2<-rechist + recovr + reshist + resovr+
  plot_layout(design=layout)
Fig2

### Figure 3: Coef Plots ##############################

summary(finr4.1)#rec
summary(fins4.1)#res

finr4.1mp<-parameters::model_parameters(finr4.1)
fmp=as.data.frame(finr4.1mp)
fmp

fmplot<-ggplot(fmp, aes(x=Coefficient, y= Parameter))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbarh(aes(xmin=CI_low, xmax= CI_high, y=Parameter),height=0.5)+
  theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Recovery Rate')
fmplot


fmpres<-parameters::model_parameters(fins4.1)
fmpresmp<-as.data.frame(fmpres)
fmpresmp

fmplotres<-ggplot(fmpresmp, aes(x=Coefficient, y= Parameter))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbarh(aes(xmin=CI_low, xmax= CI_high, y=Parameter),height=0.5)+
  theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('Resistance')

fmplotres

fmplot/fmplotres + plot_annotation(tag_levels = 'A')

### Figure 4: HII and WCS score x Recovery and Resistance ####################

#To make this figure you must first run the final recovery and resistance linear models
#Then use ggeffects to extract the data we want

##RECOVERY HII
#need the recov dataframe to run this
head(recov)
finalmodel2<-lmer(calculated.recovery.rate ~ region+disturbance+hii100km2+MPA_during_recov+recovery_time2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = recov, REML=TRUE)
#note this is called finr4.1 in the model selection code

#run the linear model
# Extract the prediction data frame
pred.hii <- ggpredict(finalmodel2, terms = c("hii100km2"))  # this gives overall predictions for the model
pred.hii

# Plot the predictions 
ggpredhii<-ggplot(pred.hii) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.3) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = hii100km2, y = calculated.recovery.rate),alpha=0.5) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Human Influence Index", y = "Recovery Rate")+
  #labs(color= "MPA Status", shape= "MPA Status")+
  #scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_aaas(labels=c('Outside','Inside'))+
  scale_fill_aaas()

ggpredhii

##RESISTANCE HII
#need the resist dataframe to run this
#run the linear model
resfinalmodel2<-lmer(resistance ~ region+disturbance+hii100km2+MPA_during_resist+resistance_time_2+cml_scr2+travel_time.tt_pop2+pre.dist.cc + (1|study), data = resist, REML=TRUE)
#note this is called fins4.1 in the model selection code

# Extract the prediction data frame
pred.hiires <- ggpredict(resfinalmodel2, terms = c("hii100km2"))  # this gives overall predictions for the model

# Plot the predictions 
ggpredhiires<-ggplot(pred.hiires) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = resist,                      # adding the raw data (scaled values)
             aes(x = hii100km2, y = resistance),alpha=0.3) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Human Influence Index", y = "Resistance")+
  #labs(color= "MPA Status", shape= "MPA Status")+
  scale_color_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_aaas(labels=c('Outside','Inside'))+
  scale_fill_aaas()

ggpredhiires


##just HII (recov + resist)
ggpredhii
ggpredhiires

Fig4hii<-(ggpredhii / ggpredhiires)+
  plot_layout(widths=c(1,1))+
  plot_annotation(tag_level='A')+
  plot_layout(guides='collect')
Fig4hii


#### RECOVERY WCS CML SCORE
#run the linear model
# Extract the prediction data frame
head(recov)
pred.wcs <- ggpredict(finalmodel2, terms = c("cml_scr2"))  # this gives overall predictions for the model
pred.wcs

ggpredwcs<-ggplot(pred.wcs) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.3) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = cml_scr2, y = calculated.recovery.rate),alpha=0.5) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled WCS Cumulative Stress Score", y = "Recovery Rate")+
  #labs(color= "MPA Status", shape= "MPA Status")+
  #scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_aaas(labels=c('Outside','Inside'))+
  scale_fill_aaas()

ggpredwcs

##RESISTANCE WCS SCORE
# Extract the prediction data frame
pred.wcsres <- ggpredict(resfinalmodel2, terms = c("cml_scr2"))  # this gives overall predictions for the model

# Plot the predictions 
ggpredwcsres<-ggplot(pred.wcsres) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = resist,                      # adding the raw data (scaled values)
             aes(x = cml_scr2, y = resistance),alpha=0.3) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled WCS Cumulative Stress Score", y = "Resistance")+
  #labs(color= "MPA Status", shape= "MPA Status")+
  scale_color_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_aaas(labels=c('Outside','Inside'))+
  scale_fill_aaas()

ggpredwcsres

###JUST WCS SCORE
ggpredwcs
ggpredwcsres

Fig4wcs<-(ggpredwcs / ggpredwcsres)+
  plot_layout(widths=c(1,1))+
  plot_annotation(tag_level='A')+
  plot_layout(guides='collect')
Fig4wcs

#Patchwork Fig (final fig)
ggpredhii
ggpredhiires
ggpredwcs
ggpredwcsres

Fig4<-(ggpredhii | ggpredhiires) / (ggpredwcs | ggpredwcsres)+  
  plot_layout(widths=c(1,1,1,1))+
  plot_annotation(tag_level='A')+
  plot_layout(guides='collect')
Fig4
######################################################################

### Figure S1: Maps ##################################################


#BUILDING MAPS
pal <- wes_palette("Zissou1", 100, type = "continuous")

pal3 <- wes_palette("Zissou1", 5, type = "continuous")



#extract world hi res map data for plotting
# world<-map_data("worldHires")
# str(world)

world1 <- ne_countries(scale = "medium", returnclass = "sf")
str(world1)

#recovery sites by hii
str(recov)
recov$hii100km=as.numeric(recov$hii100km)

#global
rechii<-
  #ggplot()+
  ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-165, 160), ylim = c(-31, 31), expand = FALSE)+
  #geom_polygon(data=world1, aes(x=lon, y=lat, group= group), fill='lightgrey', color='black')+ #makes world map
  #geom_polygon(data=world1, aes(x=long, y=lat, group= group), fill='lightgrey', color='black')+ #makes world map
  geom_point(data=recov, aes(x=long, y=lat, color=hii100km, fill=hii100km),size=4, stroke=2, position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_color_gsea(name="Human Influence Index\n             100Km")+
  scale_colour_gradientn(colours = pal,name="Human Influence Index\n             100Km")+ 
  scale_fill_gradientn(colours = pal,name="Human Influence Index\n             100Km")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  #xlim(8, 35) +  
  #ylim(-40, 35) +
  #coord_fixed(1.3)+
  #scale_colour_gradient(low='forestgreen', high='firebrick1')+ #diverging (from middle) gradient
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))
rechii
#res-> 12.0x5.4

#Making paneled figures for REC

#RECOV MAP
summary(recov$calculated.recovery.rate)

#Note: I will use hii as a size option here. Not log(hii) which makes it harder to see differences. 

recmap<-
  #ggplot()+
  ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-165, 160), ylim = c(-31, 31), expand = FALSE)+
  #geom_polygon(data=world, aes(x=long, y=lat, group= group), fill='lightgrey', color='black')+ #makes world map
  geom_point(data=recov, aes(x=long, y=lat, color=calculated.recovery.rate, size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_color_gsea(name="Recovery Rate")+
  scale_colour_gradientn(colours = rev(pal),name="Recovery Rate")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  # scale_shape_discrete(name  ="MPA Status",
  #                      breaks=c("Inside", "Outside"),
  #                      labels=c("Inside", "Outside"))+
  #xlim(8, 35) +  
  #ylim(-40, 35) +
  #coord_fixed(1.3)+
  #scale_colour_gradient(low='forestgreen', high='firebrick1')+ #diverging (from middle) gradient
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))
recmap

#Carib
library(RColorBrewer)
newpal<-brewer.pal(n=11, name='RdYlBu')

recc<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-95, -58), ylim = c(8, 30), expand = FALSE)+
  geom_point(data=recov, aes(x=long, y=lat, color=sin(calculated.recovery.rate), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_color_gradientn(colors = rainbow(20))+
  #scale_color_gradient2()+
  scale_colour_gradientn(colours = rev(pal), name="Sin(Recovery Rate)")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))+
  ggtitle('Caribbean Sea')
recc

cos(recov$calculated.recovery.rate)


#Indian
reci<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(30,123), ylim = c(-33,15), expand = FALSE)+
  geom_point(data=recov, aes(x=long, y=lat, color=sin(calculated.recovery.rate), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_colour_gradientn(colours = newpal, name="Recovery Rate")+ 
  scale_colour_gradientn(colours = rev(pal),name="Recovery Rate")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12),legend.position = 'none')+
  ggtitle('Indian Ocean')

reci

#W Pac
recw<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(120,155), ylim = c(-25, 30), expand = FALSE)+
  geom_point(data=recov, aes(x=long, y=lat, color=sin(calculated.recovery.rate), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_colour_gradientn(colours = newpal, name="Recovery Rate")+ 
  scale_colour_gradientn(colours = rev(pal),name="Recovery Rate")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12),legend.position = 'none')+
  ggtitle('W. Pacific')

recw

#E Pac
rece<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179,-79), ylim = c(-20, 25), expand = FALSE)+
  geom_point(data=recov, aes(x=long, y=lat, color=sin(calculated.recovery.rate), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  #scale_colour_gradientn(colours = newpal, name="Recovery Rate")+ 
  scale_colour_gradientn(colours = rev(pal),name="Recovery Rate")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12),legend.position = 'none')+
  ggtitle('E. Pacific')

rece


###Making paneled figures for RES
#RESIST MAP
summary(resist$resistance)
resist$hii100km

resmap<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-165, 160), ylim = c(-31, 31), expand = FALSE)+
  geom_point(data=resist, aes(x=long, y=lat, color=sin(resistance), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  scale_colour_gradientn(colours = rev(pal),name="Resistance")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))
resmap

#Carib
resc<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-95, -58), ylim = c(8, 30), expand = FALSE)+
  geom_point(data=resist, aes(x=long, y=lat, color=sin(resistance), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  scale_colour_gradientn(colours = rev(pal),name="Sin(Resistance)")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))+
  ggtitle('Caribbean Sea')
resc

#Indian
resi<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(30,123), ylim = c(-33,15), expand = FALSE)+
  geom_point(data=resist, aes(x=long, y=lat, color=sin(resistance), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  scale_colour_gradientn(colours = rev(pal),name="Resistance")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12), legend.position = 'none')+
  ggtitle('Indian Ocean')

resi

#W Pac
resw<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(120,155), ylim = c(-25, 30), expand = FALSE)+
  geom_point(data=resist, aes(x=long, y=lat, color=sin(resistance), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  scale_colour_gradientn(colours = rev(pal),name="Resistance")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12), legend.position = 'none')+
  ggtitle('W. Pacific')
resw

#E Pac
rese<-ggplot(data=world1,)+
  geom_sf(color=0) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(xlim = c(-179,-79), ylim = c(-20, 25), expand = FALSE)+
  geom_point(data=resist, aes(x=long, y=lat, color=sin(resistance), size=hii100km), stroke=2,position = position_jitter(width=1, height=1),alpha=0.6)+
  scale_colour_gradientn(colours = rev(pal),name="Sin(Resistance)")+ 
  scale_shape_manual(values=c(19,1), name  ="MPA Status",
                     breaks=c("Inside", "Outside"),
                     labels=c("Inside", "Outside"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12),legend.position = 'none')+
  ggtitle('E. Pacific')

rese  



######################

#paneling the figures together w/ patchwork

FigS1A<-((recc/reci/rece) | (recw))+
  plot_layout(guides='collect')+
  plot_layout(widths=c(1,1),heights=c(1,1))
FigS1A

FigS1B<-((resc/resi/rese) | (resw))+
  plot_layout(guides='collect')+
  plot_layout(widths=c(1,1),heights=c(1,1))
FigS1B

#patchwork of A and B together was difficult to view. A and B were instead combined using Inkscape

##########################################################################################################

### Figure S2: multipanel model results ##################################################################

#NOTE: need models to be run to get the predictions we need to run these. 

# Extract the prediction data frame
pred.hii <- ggpredict(finalmodel2, terms = c("hii100km2"))  # this gives overall predictions for the model
pred.hii
summary(pred.hii)

head(recov)

#plot predictions
ggpred<-ggplot(pred.rec) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = recovery_time2, y = calculated.recovery.rate, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope  
  theme_bw()+
  labs(x = "Scaled Recovery Time", y = "Recovery")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()
ggpred


# Plot the predictions 
ggpredhii<-ggplot(pred.hii) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = hii100km2, y = calculated.recovery.rate, colour = region)) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Human Influence Index", y = "Recovery Rate")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpredhii

###travel time to nearest pop
head(recov)
pred.dist <- ggpredict(finalmodel2, terms = c("travel_time.tt_pop2"))  # this gives overall predictions for the model
pred.dist
ggpreddist<-ggplot(pred.dist) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = travel_time.tt_pop2, y = calculated.recovery.rate, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Travel time", y = "Recovery Rate")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpreddist

###pre.dist.cc
pred.precc <- ggpredict(finalmodel2, terms = c("pre.dist.cc"))  # this gives overall predictions for the model
pred.precc

ggprecc<-ggplot(pred.precc) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = pre.dist.cc, y = calculated.recovery.rate, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Pre-Disturbance Coral Cover (%)", y = "Recovery Rate")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside', 'Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggprecc




###Predictions of factor effects

###REGION
pred.region <- ggpredict(finalmodel2, terms = c("region"))  # this gives overall predictions for the model
summary(pred.region)
pred.region
pred.region$x<-factor(pred.region$x, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

ggpredreg<-ggplot(pred.region) + 
  geom_point(data=recov, aes(x=region, y= calculated.recovery.rate), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, colour = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high, color=x),width=0.5)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  theme_bw()+
  labs(x = "Region", y = "Recovery Rate")+
  #facet_grid(~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank())+
  scale_color_jama()+
  scale_fill_jama()

ggpredreg


####DISTURBANCE
pred.disturbance<- ggpredict(finalmodel2, terms=c("disturbance"))
pred.disturbance$x<-factor(pred.disturbance$x, levels=c("Bleaching", "Disease", "Storm","COTS", "Bleaching, Disease", "Bleaching, Storm", "Bleaching, COTS"))
summary(pred.disturbance)
pred.disturbance
pred.disturbance$conf.low
pred.disturbance=as.data.frame(pred.disturbance)
pred.disturbance

ggpreddisturb<-ggplot(pred.disturbance) + 
  geom_point(data=recov, aes(x=disturbance, y= calculated.recovery.rate), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.5)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  geom_text_repel(aes(x=x, y=conf.low, label=x),size=3,force=1,vjust=0, direction='y', nudge_y=-1, segment.size=0.0, segment.colour = 'White')+
  theme_bw()+
  labs(x = "Disturbance", y = "Recovery Rate")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ylim(-12,5)

ggpreddisturb

##final recplots including coral cover
recplots<-ggpred + ggpreddist+ ggprecc+ ggpredreg  +ggpreddisturb+
  plot_layout(ncol=3)+
  plot_annotation(tag_levels='A')+
  plot_layout(guides='collect')
recplots 



#### PLOTS FOR RESISTANCE
resist$region<-factor(resist$region, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))


###RESISTANCE TIME
# Extract the prediction data frame
head(resist)
pred.res <- ggpredict(resfinalmodel2, terms = c("resistance_time_2"))  # this gives overall predictions for the model
summary(pred.res)
pred.res


# Plot the predictions 
ggpred1<-ggplot(pred.res) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = resist,                      # adding the raw data (scaled values)
             aes(x = resistance_time_2, y = resistance, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope  
  theme_bw()+
  labs(x = "Scaled Resistance Time", y = "Resistance")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()
ggpred1



####Human Influence Index
# Extract the prediction data frame
pred.hiires <- ggpredict(resfinalmodel2, terms = c("hii100km2"))  # this gives overall predictions for the model
pred.hiires

# Plot the predictions 
ggpredhiires<-ggplot(pred.hiires) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = resist,                      # adding the raw data (scaled values)
             aes(x = hii100km2, y = resistance, colour = region)) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Human Influence Index", y = "Resistance")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpredhiires

###distance to shore
head(resist)
pred.distres <- ggpredict(resfinalmodel2, terms = c("travel_time.tt_pop2"))  # this gives overall predictions for the model


ggpreddistres<-ggplot(pred.distres) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = resist,                      # adding the raw data (scaled values)
             aes(x = distance_to_shore_m2, y = resistance, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Travel Time", y = "Resistance")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpreddistres

###pre.dist.cc
pred.preccres <- ggpredict(resfinalmodel2, terms = c("pre.dist.cc"))  # this gives overall predictions for the model
pred.preccres

ggpreccres<-ggplot(pred.preccres) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = pre.dist.cc, y = resistance, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Pre-Disturbance Coral Cover (%)", y = "Resistance")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpreccres

summary(resfinalmodel2)



###Predictions of factor effects

###REGION
pred.regionres <- ggpredict(resfinalmodel2, terms = c("region"))  # this gives overall predictions for the model
summary(pred.regionres)
pred.regionres
pred.regionres$x<-factor(pred.regionres$x, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

ggpredregres<-ggplot(pred.regionres) + 
  geom_point(data=resist, aes(x=region, y= resistance), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, colour = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high, color=x),width=0.5)+
  theme_bw()+
  labs(x = "Region", y = "Resistance (% cover loss due to disturbance)")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank())+
  scale_color_jama()+
  scale_fill_jama()

ggpredregres


####DISTURBANCE
pred.disturbanceres<- ggpredict(resfinalmodel2, terms=c("disturbance"))
levels(pred.disturbanceres$x)
pred.disturbanceres$x<-factor(pred.disturbanceres$x, levels=c("Bleaching", "Disease", "Storm","COTS", "Bleaching, Disease", "Bleaching, Storm", "Bleaching, COTS"))
pred.disturbanceres

ggpreddisturbres<-ggplot(pred.disturbanceres) + 
  geom_point(data=resist, aes(x=disturbance, y= resistance), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.5)+
  geom_text_repel(aes(x=x, y=conf.low, label=x),size=3,force=1,vjust=0, direction='y', nudge_y=-2, segment.size=0.0, segment.colour = 'White')+
  theme_bw()+
  labs(x = "Disturbance", y = "Resistance")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())
#ylim(-30, 20)

ggpreddisturbres

resplots<-ggpred1 + ggpreddistres + ggpreccres+ ggpredregres +ggpreddisturbres+
  plot_layout(ncol=3)+
  plot_annotation(tag_levels='A')+
  plot_layout(guides='collect')
resplots  

##FINAL FIG S1
FigS2<-(recplots/resplots)+
  plot_layout(widths=c(1,1), heights = c(1,1))+
  plot_annotation(tag_level='A')+
  plot_layout(guides='collect')
FigS2

##################################################################################################

### Fig S3: MPA status ###########################################################################


pred.mpa<- ggpredict(finalmodel2, terms=c("MPA_during_recov"))
pred.mpa
summary(pred.mpa)
summary(finalmodel)

ggpredmpa<-ggplot(pred.mpa) + 
  geom_point(data=recov, aes(x=MPA_during_recov, y= calculated.recovery.rate), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, shape = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.3)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  theme_bw()+
  labs(x = "MPA Status")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(y='Recovery Rate')+
  scale_x_discrete(labels=c('Outside','Inside'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))
#theme(axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggpredmpa


pred.mpares<- ggpredict(resfinalmodel2, terms=c("MPA_during_resist"))
pred.mpares
summary(resfinalmodel2)

ggpredmpares<-ggplot(pred.mpares) + 
  geom_point(data=resist, aes(x=MPA_during_recov, y= resistance), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, shape = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.3)+
  theme_bw()+
  labs(x = "MPA Status")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(y='Resistance')+
  scale_x_discrete(labels=c('Outside','Inside'))+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))

ggpredmpares



pred.mpaigr<- ggpredict(finr4.1igr, terms=c("MPA_during_recov"))

ggpredmpaigr<-ggplot(pred.mpaigr) + 
  geom_point(data=recov, aes(x=MPA_during_recov, y= IGR), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, shape = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.3)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  theme_bw()+
  labs(x = "MPA Status")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(y='IGR')+
  scale_x_discrete(labels=c('Outside','Inside'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text  = element_text(size=12),axis.title = element_text(size=14),
        legend.title = element_text(size=14), legend.text = element_text(size= 12))
#theme(axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggpredmpaigr

# all 3 together

FigS3<-(ggpredmpa / ggpredmpaigr / ggpredmpares)+
  plot_layout(widths=c(1,1,1))+
  plot_annotation(tag_level='A')
FigS3

#################################################################################################

### Fig S4: Relative Recovery and Resistance plot
summary(finr4.1)#rec
summary(fins4.1)#res

finr4.1mp<-parameters::model_parameters(finr4.1)
fmp=as.data.frame(finr4.1mp)
fmp

fmplot<-ggplot(fmp, aes(x=Coefficient, y= Parameter))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbarh(aes(xmin=CI_low, xmax= CI_high, y=Parameter),height=0.5)+
  theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
fmplot


fmpres<-parameters::model_parameters(fins4.1)
fmpresmp<-as.data.frame(fmpres)
fmpresmp

fmplotres<-ggplot(fmpresmp, aes(x=Coefficient, y= Parameter))+
  geom_vline(xintercept=0, linetype='dashed')+
  geom_point()+
  geom_errorbarh(aes(xmin=CI_low, xmax= CI_high, y=Parameter),height=0.5)+
  theme_bw()+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fmplotres

FigS4<-fmplot/fmplotres + plot_annotation(tag_levels = 'A')

FigS4

########################################################################################

### Fig S5: IGR Histo ##################################################################

##Build the IGR Histogram
recov$region<-factor(recov$region, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

#overall IGR histo
igrovr<-ggplot(recov, aes(x=IGR))+
  geom_density(fill='gray')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "IGR", y='Frequency')+
  theme(axis.text = element_text(size=10))+
  theme(axis.title = element_text(size=12))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  theme(panel.border=element_blank())

igrovr

#by region

#calculate outliers for IGR and generate a new column (outlier)
igrouts<-recov%>%
  group_by(region) %>%
  mutate(outlier = ifelse(is_outlier(IGR), location, as.numeric(NA)))
igrouts
#write.csv(recouts, 'recouts.csv')



#calculate IQRs for each region (Q1 and Q3)
igrouts1<-recov%>%
  group_by(region)%>%
  summarise_at(vars(IGR),
               list(IQR=IQR, Q1=~quantile(.,probs=0.25), median=median, Q3=~quantile(.,probs=0.75)))
igrouts1
#add outlier lines (1.5*Q1 or Q3)
igrouts1$outlierval=1.5*igrouts1$IQR
igrouts1$highout=igrouts1$Q3+igrouts1$outlierval
igrouts1$lowout=igrouts1$Q1-igrouts1$outlierval
head(igrouts1)

head(igrouts)
#merge iqr info to recouts

igrouts2<-merge(igrouts, igrouts1, by='region')
head(igrouts2)

#histo for IGR by region (can add outliers with commented out code)
IGRhist<-ggplot(igrouts, aes(x=IGR, y=region, fill=region))+
  geom_density_ridges(jittered_points=TRUE, position=position_points_jitter(width=0.1, height=0), point_shape='|',alpha=0.8, quantile_lines=TRUE)+
  geom_segment(data=igrouts1, aes(x=lowout, xend=lowout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  geom_segment(data=igrouts1, aes(x=highout, xend=highout, y=as.numeric(region), yend=as.numeric(region)+.9),color='red')+
  theme_ridges()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  labs(x = "IGR", y='Region')+
  scale_color_jama()+
  scale_fill_jama()
IGRhist

layout1<- c(
  patchwork::area(t=2, l=1, b=9, r=8),
  patchwork::area(t = 1, l = 6, b = 2, r=8))

FigS5<-IGRhist + igrovr+
  plot_layout(design=layout1)
FigS5

#################################################################################################

### Fig S6: IGR multipanel ######################################################################


# Extract the prediction data frame
pred.recigr <- ggpredict(finr4.1igr, terms = c("recovery_time2"))  # this gives overall predictions for the model
pred.recigr
summary(pred.recigr)

head(recov)

#plot predictions
ggpredigr<-ggplot(pred.recigr) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = recovery_time2, y = IGR, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope  
  theme_bw()+
  labs(x = "Scaled Recovery Time", y = "IGR")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()
ggpredigr


#HII

pred.hiiigr <- ggpredict(finr4.1igr, terms = c("hii100km2"))  # this gives overall predictions for the model
pred.hiiigr
summary(pred.hiiigr)

ggpredhii<-ggplot(pred.hiiigr) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = hii100km2, y = IGR, colour = region)) +
  geom_line(aes(x = x, y = predicted),size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Human Influence Index", y = "IGR")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpredhiiigr

###travel time to nearest pop
head(recov)
pred.distigr <- ggpredict(finr4.1igr, terms = c("travel_time.tt_pop2"))  # this gives overall predictions for the model
pred.distigr
ggpreddistigr<-ggplot(pred.distigr) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = travel_time.tt_pop2, y = IGR, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Scaled Travel time", y = "IGR")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside','Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpreddistigr

###pre.dist.cc
pred.preccigr <- ggpredict(finr4.1igr, terms = c("pre.dist.cc"))  # this gives overall predictions for the model
pred.preccigr

ggpreccigr<-ggplot(pred.preccigr) + 
  geom_ribbon(aes(x = x, ymin = conf.low, ymax = conf.high), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = recov,                      # adding the raw data (scaled values)
             aes(x = pre.dist.cc, y = IGR, colour = region)) +
  geom_line(aes(x = x, y = predicted), size=1, color='black') +          # slope
  theme_bw()+
  labs(x = "Pre-Disturbance Coral Cover (%)", y = "IGR")+
  labs(color= "Region", shape= "MPA Status")+
  scale_shape_discrete(labels=c('Outside', 'Inside'))+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_jama()+
  scale_fill_jama()

ggpreccigr




###Predictions of factor effects

###REGION
pred.regionigr <- ggpredict(finr4.1igr, terms = c("region"))  # this gives overall predictions for the model
summary(pred.regionigr)
pred.regionigr
pred.regionigr$x<-factor(pred.regionigr$x, levels= c("Caribbean","Indian Ocean", "W. Pacific", "E. Pacific"))

ggpredregigr<-ggplot(pred.regionigr) + 
  geom_point(data=recov, aes(x=region, y= IGR), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted, colour = x),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high, color=x),width=0.5)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  theme_bw()+
  labs(x = "Region", y = "IGR")+
  #facet_grid(~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank())+
  scale_color_jama()+
  scale_fill_jama()

ggpredregigr


####DISTURBANCE
pred.disturbanceigr<- ggpredict(finr4.1igr, terms=c("disturbance"))
pred.disturbanceigr$x<-factor(pred.disturbanceigr$x, levels=c("Bleaching", "Disease", "Storm","COTS", "Bleaching, Disease", "Bleaching, Storm", "Bleaching, COTS"))
summary(pred.disturbanceigr)
pred.disturbanceigr
pred.disturbanceigr$conf.low
pred.disturbanceigr=as.data.frame(pred.disturbanceigr)
pred.disturbanceigr

ggpreddisturbigr<-ggplot(pred.disturbanceigr) + 
  geom_point(data=recov, aes(x=disturbance, y= IGR), color='grey', alpha=0.3)+
  geom_point(aes(x = x, y = predicted),size=3) +
  geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high),width=0.5)+
  #geom_errorbar(aes(x=x, ymin=predicted-std.error, ymax=predicted+std.error),width=0.3)+
  geom_text_repel(aes(x=x, y=conf.low, label=x),size=3,force=1,vjust=0, direction='y', nudge_y=-1, segment.size=0.0, segment.colour = 'White')+
  theme_bw()+
  labs(x = "Disturbance", y = "IGR")+
  #facet_grid(region~MPA_status)+
  #theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8))+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(axis.title.y=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ylim(-12,5)

ggpreddisturbigr

##final IGRplots including coral cover
igrplots<-ggpredigr + ggpreddistigr+ ggpreccigr+ ggpredregigr  +ggpreddisturbigr+
  plot_layout(ncol=3)+
  plot_annotation(tag_levels='A')+
  plot_layout(guides='collect')
igrplots 


### adding mpa status,  ggpredhiiigr, and ggpredwcsigr to the mix (with IGR ONLY)
igrplots2<-ggpredigr + ggpreddistigr+ ggpreccigr+ ggpredregigr  +ggpreddisturbigr+ ggpredhiiigr + ggpredwcsigr+ggpredmpaigr+
  plot_layout(ncol=3)+
  plot_annotation(tag_levels='A')+
  plot_layout(guides='collect')
igrplots2 

FigS6 <-igrplots2

FigS6

















