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

G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FST_list_Florida_reorder_noKeys.csv", sep=",", header=T) 

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

print("spatial lines done")
###############################################
#Create raster stack 
###############################################

aridI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/Florida_clips/AI_annual_FloridaClip.tif")
arid = aridI*1
proj4string(arid) <- crs.geo

accessI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips/accessibility_to_cities_2015_v1.0_FloridaClip.tif")
access <- accessI*1
proj4string(access) <- crs.geo

precI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/Florida_clips/bio12_mean_FloridaClip.tif")
prec <- precI*1
proj4string(prec) <- crs.geo

mean.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/Florida_clips/bio1_mean_FloridaClip.tif")
mean.temp <- mean.tempI*1
proj4string(mean.temp) <- crs.geo

human.densityI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/Florida_clips/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip.tif")
human.density <- human.densityI*1
proj4string(human.density) <- crs.geo

frictionI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/Florida_Clips/friction_surface_2015_v1.0_FloridaClip.tif")
friction <- frictionI*1
proj4string(friction) <- crs.geo

min.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/Florida_clips/bio6_mean_FloridaClip.tif")
min.temp <- min.tempI*1
proj4string(min.temp) <- crs.geo

NeedleleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_1_FloridaClip.tif")
Needleleaf = NeedleleafI*1
proj4string(Needleleaf) <- crs.geo

EvBroadleafI  = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_2_FloridaClip.tif")
EvBroadleaf = EvBroadleafI*1
proj4string(EvBroadleaf) <- crs.geo

DecBroadleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_3_FloridaClip.tif")
DecBroadleaf = DecBroadleafI*1
proj4string(DecBroadleaf) <- crs.geo

MiscTreesI =  raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_4_FloridaClip.tif")
MiscTrees = MiscTreesI*1
proj4string(MiscTrees) <- crs.geo

ShrubsI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_5_FloridaClip.tif")
Shrubs = ShrubsI*1
proj4string(Shrubs) <- crs.geo

HerbI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_6_FloridaClip.tif")
Herb = HerbI*1
proj4string(Herb) <- crs.geo

CropI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_7_FloridaClip.tif")
Crop <- CropI*1
proj4string(Crop) <- crs.geo

FloodI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_8_FloridaClip.tif")
Flood = FloodI*1
proj4string(Flood) <- crs.geo

UrbanI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_9_FloridaClip.tif")
Urban <- UrbanI*1
proj4string(Urban) <- crs.geo

SnowI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_10_FloridaClip.tif")
Snow = SnowI*1
proj4string(Snow) <- crs.geo

BarrenI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_11_FloridaClip.tif")
Barren = BarrenI*1
proj4string(Barren) <- crs.geo

WaterI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips/consensus_full_class_12_FloridaClip.tif")
Water = WaterI*1
proj4string(Water) <- crs.geo

SlopeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/Florida_clips/slope_1KMmedian_MERIT_FloridaClip_positive.tif")
Slope = SlopeI*1
proj4string(Slope) <- crs.geo

AltitudeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/Florida_clips/altitude_1KMmedian_MERIT_FloridaClip.tif")
Altitude = AltitudeI*1
proj4string(Altitude) <- crs.geo

PETI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET/Florida_clips/pet_mean_FloridaClip.tif")
PET = PETI*1
proj4string(PET) <- crs.geo

DailyTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio2/Florida_clips/bio2_mean_FloridaClip.tif")
DailyTempRange = DailyTempRangeI*1
proj4string(DailyTempRange) <- crs.geo

max.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5/Florida_clips/bio5_mean_FloridaClip.tif")
max.temp = max.tempI*1
proj4string(max.temp) <- crs.geo

AnnualTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio7/Florida_clips/bio7_mean_FloridaClip.tif")
AnnualTempRange = AnnualTempRangeI*1
proj4string(AnnualTempRange) <- crs.geo

prec.wetI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio13/Florida_clips/bio13_mean_FloridaClip.tif")
prec.wet = prec.wetI*1
proj4string(prec.wet) <- crs.geo

prec.dryI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio14/Florida_clips/bio14_mean_FloridaClip.tif")
prec.dry = prec.dryI*1
proj4string(prec.dry) <- crs.geo

GPPI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly_mean/Florida_clips/GPP_mean_Florida_clip.tif")
GPP = GPPI*1
proj4string(GPP) <- crs.geo

env=stack(arid, access, prec, mean.temp, human.density, friction, min.temp, Needleleaf, EvBroadleaf, DecBroadleaf, MiscTrees, Shrubs, Herb, Crop, Flood, Urban, Snow, Barren, Water, Slope, Altitude, PET, DailyTempRange, max.temp, AnnualTempRange, prec.wet, prec.dry, GPP)

names(env) [1] <- "arid"
names(env) [2] <- "access"
names(env) [3] <- "prec"
names(env) [4] <- "mean.temp"
names(env) [5] <- "human.density"
names(env) [6] <- "friction"
names(env) [7] <- "min.temp"
names(env) [8] <- "Needleleaf" 
names(env) [9] <- "EvBroadleaf"
names(env) [10] <- "DecBroadleaf"
names(env) [11] <- "MiscTrees"
names(env) [12] <- "Shrubs"
names(env) [13] <- "Herb"
names(env) [14] <- "Crop"
names(env) [15] <- "Flood"
names(env) [16] <- "Urban"
names(env) [17] <- "Snow"
names(env) [18] <- "Barren"
names(env) [19] <- "Water"
names(env) [20] <- "Slope"
names(env) [21] <- "Altitude"
names(env) [22] <- "PET"
names(env) [23] <- "DailyTempRange"
names(env) [24] <- "max.temp"
names(env) [25] <- "AnnualTempRange"
names(env) [26] <- "prec.wet"
names(env) [27] <- "prec.dry"
names(env) [28] <- "GPP"


print("raster stack done")

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/sc05E_sc06E_rasterstack_image.RData")


########################################
#Calculate mean of straight lines and making initial RF model
#######################################
StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

StraightMeanDF <- as.data.frame(StraightMean)

StraightMeanDF$FST_arl <- G.table$FST_arl

#option of trying DPS
#StraightMeanDF$DPS <- G.table$DPS
  
Straight_RF = randomForest(FST_arl ~   arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

#Straight_RF_text = Straight_RF

#write.csv(Straight_RF_text, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/Straight_RF.csv")


StraightPred <- predict(env, Straight_RF)

print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL_points_list_noKeys.csv", sep=",", header=T)
P.coordinates1 <- P.table[,c(3,2)]
P.points <- SpatialPoints(P.table[,c(3,2)])  # ... converted into a spatial object
proj4string(P.points) <- crs.geo  
#plot(P.points)

print("starting loops")

it <- 1
for (it in 1:10) {
  
  trFlorida <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix
  trFloridaC <- geoCorrection(trFlorida, type="c") 

  AtoT <- shortestPath(trFloridaC, P.points[1,], P.points[1,], output="SpatialLines")
  for (x in 1:10) {  
    for (y in (x+1):11) { 
     Ato <- shortestPath(trFloridaC, P.points[x,], P.points[y,], output="SpatialLines")
     AtoT <- AtoT + Ato
     assign(paste0("Ato_px",x,"_py",y,"_it",it) , Ato)   

    }
  }
  AtoT = (AtoT[-1])

  LcpLoop <- raster::extract(env, AtoT, fun=mean, na.rm=TRUE)

  LcpLoopDF <- as.data.frame(LcpLoop)

  LcpLoopDF$FST_arl <- G.table$FST_arl

  LCP_RF = randomForest(FST_arl ~  arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf +  MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP, importance=TRUE, na.action=na.omit, data=LcpLoopDF)

  assign(paste0("LCP_RF", it), LCP_RF )

  pred = predict(env, LCP_RF)

  assign(paste0("pred", it), pred)
  
  pred.cond <- 1/pred 
  
  print("round done")

}

#make table of RF variance explained

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/sc05E_sc06E.RData")

test = summary(LCP_RF)

write.csv(test, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/test.csv")