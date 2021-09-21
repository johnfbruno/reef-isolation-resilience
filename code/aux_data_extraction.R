# aux data extraction


###Extract data from other local stress / human influence metrics
##What is in HII? --> pop density, land use (build up areas, nighttime lights, land use, land cover), human acess (coastalines, roads, railroads, nav. rivers).
#Total parameters (8 or 9)

##WCS local reef pressures (from Andrello et al pre-print): Cumulative pressure index that includes fishing (market gravity from Cinner et al 2018--this IS the same as above metric...), coastal development (gridded pop of the world from CIESIN 2018--same as HII only updated --total human pop within 5km of each coral reef location), industrial development (World port location from Google data--total port # within 5km of each reef location), tourism (economic value of each reef loc from Spalding et al 2017), sedimentation (estaimate sed delivery using a sediment delivery ratio from the InVEST SDR model, Also used the sed deilivery ratio values from Reefs at Risk (Burke et al 2011)--used R to caluc sediment plumes from each river mouth and averaged these values for each reef location), nitrogen (Gathered crop cover data using methods from Halpern 2008, gathered N fertilizer use values for each country from FAO data (for each reef cathcment), Estimated plant N uptake and estimted N devlivery to the river mouth assuming 20% and 80% N fertilzer uptake by plants. Modeled the N plume following similar methods as SDR above, summed plume values from river mouths to provide a cumulative N exposure value. Averaged the plume values for each reef location). Cumulative score: weighted average of the percentiles of the 6 layers
#Total parameters= 7. Shared w. HII = 2.5 or so (pop den is in gravity, dev is proxied differently but in both) 
#These are data summarized for each 5kmx5km reef pixel. The code below will take my coordinates, cross reference them to a global dbf of 5kmx5km for which values were calculated and find the nearest pixel centroid. My coords are then assigned to that pixel!
#**NOTE: Gravity includes a population density component (pop within 500km radius / travel time squared) and so does coastal dev (Sum of 5km coastl pop grid cells that were within 5km of each reef grid cell)

library(here)
library(tidyverse)
library(readr)
library(rgdal)



data=read.csv('coral_recov_master.csv') #long=10, lat=9
head(data)

allreefs <- rgdal::readOGR("allreefs.shp")
summary(allreefs)

#there is a custom function in the original extract_values.R code from the paper. Most of which is below. I stepped through it step by step so I could understand how it worked. 
#I also made some changes from the orignal to match my needs!
max.radius=5000

# Read CRS for allreefs
prj4 <- sp::proj4string(allreefs)
cat("CRS for allreefs is\n",prj4,"\n")

# Convert points to spatialPoints. Assuming input data are in WGS84 (EPSG: 4326)
points <- sp::SpatialPoints(data[,c(10,9)], #can alter these columns for long, lat on any data :)
                            proj4string=sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))

# Change CRS for points
points_prj4 <- sp::spTransform(points,sp::CRS(prj4))

# Extract allreefs centroids
allreefs_centroids <- rgeos::gCentroid(allreefs, byid=T)

# Find nearest neighboring allreefs for each point
nn <- RANN::nn2(sp::coordinates(allreefs_centroids),
                sp::coordinates(points_prj4),
                k=1,
                searchtype = "radius",
                radius = max.radius)

# Initializing output dataframe
# LEGEND:
# "score", Climate: composite score
# "scorecn", Climate: connectivity
# "scorecy", Climate: Cyclone Risk
# "scorepfc", Climate: Thermal future
# "scoreth", Climate: Thermal history
# "scoretr", Climate: Recent stress
# "grav_NC", Fishing: Market Pressure
# "sediment", Pollution: Sedimentation
# "nutrient", Pollution: Nutrients
# "pop_count", Coastal Development: Human Population
# "num_ports", Industrial Development: Ports
# "reef_value", Tourism: Reef Value
out.data <- matrix(NA, nrow=nrow(nn$nn.idx), ncol = 8)
colnames(out.data) <-c("grav_NC",
                       "pop_cnt",
                       "nm_prts",
                       "reef_vl",
                       "sedimnt",
                       "nutrint",
                       "cml_scr",
                       "tp_thrt")
out.data <- as.data.frame(out.data)
# Loop on points to fill values into output dataframe
for (i in 1 : nrow(nn$nn.idx)) {
  if (nn$nn.idx[i,1] == 0) next
  out.data[i, ] <- as.list(allreefs@data[nn$nn.idx[i,1],
                                         c("grav_NC",
                                           "pop_cnt",
                                           "nm_prts",
                                           "reef_vl",
                                           "sedimnt",
                                           "nutrint",
                                           "cml_scr",
                                           "tp_thrt")])
}
out.data <- cbind(data,out.data)


write.csv(out.data, "coral_recov_master_2021.csv")

#######################################################

##Global human gravity (Cinner et al, 2018)-- combines population density and approx. travel time to location. Gives data within 500km radius
#Total parameters = 2. Shared w/ Hii= 1 (pop density), ratio= 50% shared
#***NOTE: this index is intergrated into the WCS index below. I output it as grav_nc into my dataframe. So I can compare it and the WCS data!
#**Note that this gravity metric is some kind of intergrated product that considers a 500 km radius

data=read.csv('coral_recov_master_2021.csv')
head(data) #long,lat =11,10***CHECK THIS EACH TIME

grav <- rgdal::readOGR("Total Gravity of Coral Reefs 1.0.shp")
summary(grav)

#change the coordinate reference system to WGS84 (standard)
# Define spatial projections to standard WGS84 datum
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

# Project for Human Infl Index data
proj4string(grav) <- crs.geo

#get coords of data (long, lat) 
coord<-data[,c(11,10)]
head(coord)

# Define coordinates columns
coordinates(coord)<- c("long","lat")
head(coord)

# Define spatial projections for sites
#this projection must match the map for calculations to be accurate!
proj4string(coord) <- crs.geo
summary(coord)

data$gravity <- extract(grav, coord, fun = sum)
head(data)

write.csv(data, "coral_recov_master_2021.csv")
######################################################################

##Marie et al (2016) Travel time to reefs data
data=read.csv('coral_recov_master_2021.csv')

head(data) #check for random added columns at beginning... 
data<-data[,-c(1:3)]

#long = column 9
#lat = column 8

trav <- rgdal::readOGR("Global Accessibility of Coral Reefs 2.0.shp")
summary(trav)

#backup df for if sites are too close to shore (within 10km) 
trav2<-rgdal::readOGR("FullPolygons.shp")
sumary(trav2)

#change the coordinate reference system to WGS84 (standard)
# Define spatial projections to standard WGS84 datum
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

# Project for Human Infl Index data
proj4string(trav) <- crs.geo
proj4string(trav2) <- crs.geo


#get coords of data (long, lat) 
head(data)
coord<-data[,c(7,6)]
head(coord)

# Define coordinates columns
coordinates(coord)<- c("long","lat")
head(coord)

# Define spatial projections for sites
#this projection must match the map for calculations to be accurate!
proj4string(coord) <- crs.geo
summary(coord)

#travel time from final global accessibility shapefile
data$travel_time_test <- extract(trav, coord, fun = sum)

#travel time from full polygons shapefile
data$travel_time_backup2 <- extract(trav2, coord, fun = mean)

head(data)
ncol(data)

#inspect data
dat2<-data[,c(1,2,3,4,5,6,7,65,66,67,68,69,70,71,72,73)]
head(dat2)
ncol(dat2)

dat2<-data[,c(1,2,3,4,5,70,71,72,73)]
summary(dat2)

write.csv(dat2, "coral_recov_test.csv")


write.csv(data, "coral_recov_master_2021.csv")
################################################