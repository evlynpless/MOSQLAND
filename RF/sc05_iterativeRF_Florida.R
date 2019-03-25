#Building iterative RF model with Florida as test case

library("sp")
library("spatstat")
library("maptools")
library("raster")
library("randomForest")
library("gdistance")
library("SDraw")
library("tidyverse")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system


###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

#Add a loop that will make FST_list_Florida_reorder.csv from FL_points_list.csv??

G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FST_list_Florida_reorder.csv", sep=",", header=T) 

#create dataframes of begin and end coordinates from a file:
begin.table <- G.table[,c(4,3)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- G.table[,c(6,5)]
end.coord <- end.table
coordinates(end.coord) <- c("long2", "lat2")

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data


###############################################
#Create raster stack 
###############################################

arid = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/Florida_clips/AI_annual_FloridaClip.tif")
proj4string(arid) <- crs.geo

access = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips/accessibility_to_cities_2015_v1.0_FloridaClip.tif")
proj4string(access) <- crs.geo

prec = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/Florida_clips/bio12_mean_FloridaClip.tif")
proj4string(prec) <- crs.geo

mean.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/Florida_clips/bio1_mean_FloridaClip.tif")
proj4string(mean.temp) <- crs.geo

human.density = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/Florida_clips/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip.tif")
proj4string(human.density) <- crs.geo

crop = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_7_FloridaClip.tif")
proj4string(crop) <- crs.geo

urban = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_9_FloridaClip.tif")
proj4string(urban) <- crs.geo

friction = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/Florida_Clips/friction_surface_2015_v1.0_FloridaClip.tif")
proj4string(friction) <- crs.geo

min.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/bio6_mean_FloridaClip.tif")
proj4string(friction) <- crs.geo

ABSHUM = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/Florida_clips/ABS50_res_FloridaClip.tif")
proj4string(friction) <- crs.geo

env=stack(arid, access, prec, mean.temp, human.density, crop, urban, friction, min.temp, ABSHUM)

names(env) [1] <- "arid"
names(env) [2] <- "access"
names(env) [3] <- "prec"
names(env) [4] <- "mean.temp"
names(env) [5] <- "human.density"
names(env) [6] <- "crop"
names(env) [7] <- "urban"
names(env) [8] <- "friction"
names(env) [9] <- "min.temp"
names(env) [10] <- "ABSHUM"


########################################
#Calculate mean of straight lines and making initial RF model
#######################################
StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

StraightMeanDF <- as.data.frame(StraightMean)

StraightMeanDF$FST_arl <- G.table$FST_arl

#option of trying DPS
#StraightMeanDF$DPS <- G.table$DPS
  
Straight_RF = randomForest(FST_arl ~   arid + access  +   prec  +   mean.temp  +   human.density  +   crop   +    urban  +   friction + min.temp + ABSHUM, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

StraightPred <- predict(env, Straight_RF)

StraightPred[is.na(StraightPred[])] <- 0.4 #can delete after Florida Keys highway is added

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL_points_list.csv", sep=",", header=T)
G.coordinates1 <- G.table[,c(3,2)]
G.points <- SpatialPoints(G.table[,c(3,2)])  # ... converted into a spatial object
proj4string(G.points) <- crs.geo  
plot(G.points)

it <- 1
for (it in 1:2) {
  
  trFlorida <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix
  trFloridaC <- geoCorrection(trFlorida, type="c") 

  AtoT <- shortestPath(trFloridaC, G.points[1,], G.points[1,], output="SpatialLines")
  for (x in 1:13) {  
    for (y in (x+1):14) { 
     Ato <- shortestPath(trFloridaC, G.points[x,], G.points[y,], output="SpatialLines")
      AtoT <- AtoT + Ato  
    }
  }
  AtoT = (AtoT[-1])

  LcpLoop <- raster::extract(env, AtoT, fun=mean, na.rm=TRUE)

  LcpLoopDF <- as.data.frame(LcpLoop)

  LcpLoopDF$FST_arl <- G.table$FST_arl

  LCP_RF = randomForest(FST_arl ~  arid + access  +   prec  +   mean.temp  +   human.density  +   crop   +    urban  +   friction + min.temp + ABSHUM, importance=TRUE, na.action=na.omit, data=LcpLoopDF)

  pred = predict(env, LCP_RF)
  
  pred[is.na(pred[])] <- 0.4
  
  pred.cond <- 1/pred 

  it <- it+1
  
  print("round done")
}

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/image.RData")

test = summary(LCP_RF)

write.csv(test, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/test.csv")