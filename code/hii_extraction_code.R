### HII extraction code ####

#Step 1: Read in data
a=read.csv('Coral_Resilience_database_2.0.csv') #database of coordinates used in this study. You can use any coordinates!
head(a)


#STEP 2: Human influence index
##Human Influence Index Code ##
#this code modified from Abel Valdivia's code: https://github.com/johnfbruno/Extract_GHII/blob/master/Extract%20GHII.R
#data from SEDAC: http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic
#once downloaded, data are stored in ESRI file formarts (ArcGIS)
#Cannot open the pyramid file (.rrd) in r, but you do not need to. 
#1.) Extract files to a folder
#2.) open the largest .adf file (w001001.adf) in r using the raster package

library(sp)
library(rgdal)
library(raster)
library(ggmap)
library(maps)
library(maptools)
library(mapdata)
library(ggplot2)


#navigate to directory containing .adf file and load using raster
hii_global<-raster(paste0("w001001.adf"))
hii_global #confirm that it loads, note the reference system for the coorindates is an unusual one

#change the coordinate reference system to WGS84 (standard)
# Define spatial projections to standard WGS84 datum
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

# Project for Human Infl Index data
proj4string(hii_global) <- crs.geo

# Plot data only for the Mesoamerican Barrier Region as a test
plot(hii_global, ylim=c(5,30), xlim=c(-91,-80),  main="Human Influence Index")

#get coords of data (long, lat) 
coord<-a[,c(6,5)] #columns 6 and 5 are long and lat from this data file. 
head(coord)

# Define coordinates columns
coordinates(coord)<- c("long","lat")
head(coord)

# Define spatial projections for sites
#this projection must match the map for calculations to be accurate!
proj4string(coord) <- crs.geo
summary(coord)

#mapping data to check it
#mapping data to check it
map("worldHires", xlim=c(-90.0,-80.0),ylim=c(6,30), col="gray90", fill=TRUE) 
#add coords to confirm they seem right
plot(coord, add=T,  cex=0.5, col="blue", pch=16)


###Sum of HHI within a buffer zone of the sites###
# Extract Human Influence Index for the site coordinates within 100 km buffer
a$hii100km <- extract(hii_global, coord, buffer = 100000, fun = sum)

# Extract Human Influence Index for the site coordinates within 75 km buffer
a$hii75km <- extract(hii_global, coord, buffer = 75000, fun = sum)

# Extract Human Influence Index for the workd within 50 km buffer
a$hii50km <- extract(hii_global, coord, buffer = 50000, fun = sum)

# Extract Human Influence Index for the workd within 25 km buffer
a$hii25km <- extract(hii_global, coord, buffer = 25000, fun = sum)

# Extract Human Influence Index for the workd within 10 km buffer
a$hii10km <- extract(hii_global, coord, buffer = 10000, fun = sum)

head(a)

write.csv(a, file='HHI_coral_recov.csv')