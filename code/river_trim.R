###Trimming global rivers to global river mouths###

library(ggplot2)
library(ggmap)
#detach("package:ggmap", unload=TRUE)
library(maps)
library(mapdata)
library(maptools)
library(scales) #for transparency
library(mapproj) #for projected maps
library(rgdal)
library(GISTools)  
library(devtools)
library(ggsn) #for north arrows and scale
library(ggthemes)
library(Rmisc)
library(sp)
library(raster)
library(reshape)
library(plyr)
library(sf)
library(gdistance) #used to calc geographic distance
library(rgeos) #used to calc geographic distance
library(ggmap)
library(maptools)
library(PBSmapping)



rivers<-readOGR("ne_10m_rivers_lake_centerlines.shp", layer="ne_10m_rivers_lake_centerlines")

#remove lakes from the shapefile data
rivers1<-subset(rivers, featurecla=='River')

#read in global coastline data
coasts<-readOGR("combined_coastline.shp", layer="combined_coastline")

#Trim rivers to only river mouths (where coasts and rivers intersect)
test1<-raster::intersect(rivers1, coasts) #error when run w/ rivers1 (spatiallines)

#write the resulting file out
writeOGR(obj=test1, dsn=".", layer="test1", driver="ESRI Shapefile")

