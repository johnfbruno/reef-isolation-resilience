###Distance to nearest river###
#source: https://gis.stackexchange.com/questions/225102/calculate-distance-between-points-and-nearest-polygon-in-r

#devtools::install_github("valentinitnelav/geobuffer")
library(geobuffer)
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
library(tidyverse)

#Note that source data file has changed (coral_recov_master_2021.csv will be needed). This may change the column numbers that need to be kept. 


b=read.csv('coral_recov_master.csv')
head(b)
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot() +   mapWorld
mp
Lon<-b$long
Lat<-b$lat
mp <- mp+ geom_point(aes(x=Lon, y=Lat),color="blue", size=2) 
mp




head(b)
#c<-b[ -c(1)] #remove site column
#c
#fals<-c[-1,]
#silk<-c[-2,]

crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #MAke things default to WGS84

#READ IN RIVERS SHAPEFILE
rivers<-readOGR("ne_10m_rivers_lake_centerlines.shp", layer="ne_10m_rivers_lake_centerlines")
plot(rivers)
summary(rivers)
str(rivers)
class(rivers)
str(rivers@data)
rivers2<-rivers[c(1:5),]

#remove lakes from the shapefile data
rivers1<-subset(rivers, featurecla=='River')


#make into spatial points (from spatial lines)
riverpt<-as(rivers1, "SpatialPointsDataFrame")
plot(riverpt)
str(riverpt@data)

#convert to spatial polygons
riverply<-SpatialLines2PolySet(rivers1) 
plot(riverply)


###MAKE RIVERS SHAPEFILE INTO RIVER MOUTHS SHAPEFILE!!!!####
#Try to trim river lines to river mouths using intersect in raster package and global coastline layer.
#read in global coastline data
coasts<-readOGR("combined_coastline.shp", layer="combined_coastline")
summary(coasts)

#plot coatlines and rivers (as lines or points)
plot(coasts)
#plot(rivers1, add=TRUE)
plot(riverpt, add=TRUE)


###
#Trim rivers to only river mouths (where coasts and rivers intersect)
#test1<-raster::intersect(rivers1, coasts) #This doesn't work for 2 spatial lines objects, try rgeos:gIntersection instead
test1<-rgeos::gIntersection(rivers1, coasts) #works! Makes spatialpoints which are river mouths!
plot(test1)

test2<-rgeos::gIntersection(rivers1, coasts, byid=TRUE, id)
row.names(rivers1)
row.names(coasts)
row.names(test2) #each row name is rowname of rivers1 then coasts. Can convert to data.frame and add data back in by ID later if needed. 

names(geometry(rivers1))
plot(test2)
head(test2)
str(test2)

#my points
head(b)
#c=b[,c(6,5)]
#head(c)

# #test data
# c1=b[c(1:10),]
# head(c1)

#get coordinates for points
c2<-b[,c(8,7)]
head(c2)

#coordinates(c1)<-~long+lat
#c2<-coordinates(c1)
#head(c2)

c3<-SpatialPointsDataFrame(coords=c2, data=c2, proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

#make test2 (spatialpoints) into a dataframe or extract coordinates
head(test2@coords)
t2coord=as.data.frame(test2@coords)
head(t2coord)

t2coords<-rename(t2coord, c(x="long", y="lat"))
head(t2coords)

head(test2)
head(c2)


#Distance to nearest river
#attempts to calculate distance between points seem to result in warning in geosphere package
#warning: In cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2], p3[, 1], p3[, 2], r) :
#these return matrixes? Not super useful here. 
#warning message: number of rows of result is not a multiple of vector length (arg 1)
# dist.1<-distHaversine(c2, t2coords)
# dist.2<-distGeo(c2, test2)
# dist.3<-distCosine(c2,test2)
# dist.4<-dist2gc(c2,test2, c(0,0))
# dist.m <- geosphere::distm(c3, test2)
# dist.m

#calculate distance to nearest river (in meters)
head(t2coords)

#make an empty list of rows of C2 to append distances to.
res<-numeric(nrow(c2))
#get the min distance
for (i in 1:length(res)) res[i]<-min(distHaversine(c2[i,,drop=FALSE],t2coord)) #gives distance of nearest river for each row of C2
res
#make into a data frame (column)
res=as.data.frame(res)

res1<-numeric(nrow(c2))
for (i in 1:length(res1)) res1[i]<-which.min(distHaversine(c2[i,,drop=FALSE],t2coord)) #gives index number (row number) of nearest river
res1
res1=as.data.frame(res1)

#append res and res 1 together!
res2=cbind(res, res1)
res2 #res= distance (m) to nearest river. Res1= row number (from t2coord) of nearest river. Rownumber of res2= row number of C2

#Rename columns
res2<-plyr::rename(res2, c("res"="dist_to_riv_m", "res1"="riv_row"))
head(res2)
#append c2 row names as point IDs (note, C2 is just coords of C1, so row names are the same!)
res2$point_ID=row.names(c2)

#append river data by row name
  #first must give t2coord an index column
t2coord$riv_row<-seq.int(nrow(t2coord))
head(t2coord)

  #next add ID column (to determine which river we are talking about) to t2coord
t2coord$ID=row.names(t2coord)
head(t2coord)

  #now we can join the t2coord (rivers) and res2 (nearest river to each point) together by riv_row
res3<-merge(res2, t2coord, by='riv_row')
head(res3)
res3

  #remove the x and y (long and lat of river points)
res4<-res3[,-c(4,5)]

#convert distance to km
res4$dist_to_riv_km=res4$dist_to_riv_m/1000
head(res4)

#Split ID column (river ID column= X, river ID, coast ID)
res5<-separate(data=res4, col=ID, into=c("river","coast"), sep="\\.")
res5

#remove the X from the beginning of the river ID
res5$river=sub('.', '', res5$river)
head(res5)

#add ID column to rivers1
rivers1$river=row.names(rivers1)
head(rivers1)

#trim rivers 1 to only river name and river ID
rivers2<-rivers1[,c(4,5,35)]
head(rivers2)

#Merge river nanes and distance data
fin_res<-merge(res5, rivers2, by='river')
head(fin_res)

#Append river data to master data
# head(c2)
# c2$point_ID=row.names(c2)

head(b)
b$point_ID=row.names(b)

# c3<-merge(c2, fin_res, by='point_ID')
# head(c3)

c4<-merge(b, fin_res, by='point_ID')
head(c4)
max(c4$long)
min(c4$long)

write.csv(c4, file='coral_recov_master.csv')


###Rivers within 100km 
plot(test1) #this is river mouths
head(b) #full data (includes lat and long)

#convert test1 as a raster layer
test1.2<-raster(test1)
test1.2
head(b)

head(c2)
head(c3) #spatial points dataframe of coordinates from b

#add 100 km buffer to each points in b
pts_buf_100km<-geobuffer_pts(xy=c2, dist_m=100*10^3)
pts_buf_100km #spatialpolygons
plot(pts_buf_100km)
head(pts_buf_100km)

test1 #spatialpoints

#may need to convert to spatialpoints and spatialpolygons dataframes!

crs(test1)=crs(pts_buf_100km)

res<-over(test1, pts_buf_100km)
head(res)
summary(res)
str(res)

b$riv100km<-over(test1, pts_buf_100km)
b$riv100km


#proj4string(c3)
coordinates(c2)<-c('long','lat')
proj4string(c2)<-CRS("+proj=longlat +datum=WGS84")
pc<-spTransform(c2, CRS("+proj=longlat +datum=WGS84"))

distInMeters<-1000
pc100km<-gBuffer(pc, width=100*distInMeters, byid=TRUE)



#trying to count number of points (river mouths) within a buffer (THIS DOES NOT WORK)
#b$riv100km <- raster::extract(test1.2, coord, buffer = 100000, fun = count)
#b$riv100km



# ##############Visualizing river distances#################
# ###PLUS failed attempts###################################
# 
# #visual confirmation that these distances are correct!
# pts.sp <- sp::SpatialPoints(coords      = c2[,c("long","lat")], # order matters
#                             proj4string = rivers1@proj4string)
# 
# #plot(rivers1)
# plot(test1)
# plot(pts.sp, col="red", add=TRUE)
# 
# # plot arrows to indicate the direction of the great-circle-distance
# for (i in 1:nrow(pts.wit.dist)) {
#   arrows(x0 = pts.wit.dist[i,1], 
#          y0 = pts.wit.dist[i,2], 
#          x1 = pts.wit.dist[i,4], 
#          y1 = pts.wit.dist[i,5],
#          length = 0.1,
#          col = "green")
# }
# 
# 
# #find names of rivers that are nearest (note: some will be NA). NOTE 2: Some are not river mouths. 
# 
# pts.wit.dist
# pts.wit.dist$lat1=pts.wit.dist$lat
# pts.wit.dist=pts.wit.dist[c(1,2,3,4,7,6)]
# 
# pts.wit.dist$ID
# 
# test<-rivers1[382,]
# head(test)
# 
# head(rivers1)
# tail(rivers1)
# 
# #add an ID column to rivers 1 that is just numbers 1-x (starting at 1)
# rivers1$ID <- seq.int(nrow(rivers1))
# rivs=as.data.frame(rivers1)
# 
# riv<-rivs[,c(4,5,35)]
# head(riv)
# riv[382,]
# rivs$name #note: there are some rivers with no name!
# 
# #append river name and river number onto dist.mat
# dist.mat1<-merge(pts.wit.dist, riv, by='ID')
# dist.mat1
# #pts.wit.dist
# 
# dist.mat1$distkm=dist.mat1$distance / 1000
# dist.mat1

