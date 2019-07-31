#Building iterative RF model with Florida, no Keys
#70:30 Training:testing data

library("sp")
library("spatstat")
library("maptools")
library("raster")
library("randomForest")
library("gdistance")
library("SDraw")
library("tidyverse")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

rmr=function(x){
## function to truly delete raster and temporary files associated with them
if(class(x)=="RasterLayer"&grepl("^/tmp",x@file@name)&fromDisk(x)==T){
file.remove(x@file@name,sub("grd","gri",x@file@name))
rm(x)
}
}

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

#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/FL_rasterstack_image.RData")

#load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/FL_rasterstack_image.RData")

###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/FST_list_Florida_reorder_noKeys.csv", sep=",", header=T)

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

########################################
#Calculate mean of straight lines and making initial RF model
#######################################
StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

StraightMeanDF <- as.data.frame(StraightMean)

StraightMeanDF$FST_arl <- G.table$FST_arl

#option of trying DPS
#StraightMeanDF$DPS <- G.table$DPS
 
NumPairs = nrow(StraightMeanDF)
Training = NumPairs * 0.7
TrainingInt = round(Training)
TrainingPairs = sample(1:NumPairs, TrainingInt, replace = FALSE)
StraightMeanDF.train = StraightMeanDF[TrainingPairs,]
StraightMeanDF.valid = StraightMeanDF[-TrainingPairs,]

set.seed(NULL)
 
Straight_RF = randomForest(FST_arl ~   arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP, importance=TRUE, na.action=na.omit, data=StraightMeanDF.train)

#Define empty vectors
RSQ_vec = c()
RSQeq_vec = c()
MSE_vec	    = c()
MSEeq_vec = c()
MSE2_vec  = c()
Cor1_vec  = c()
Cor2_vec  = c()

#Performance measures
RSQ = tail(Straight_RF$rsq ,1 ) 
RSQeq = 1 - (tail(Straight_RF$mse , 1) / var(StraightMeanDF.train$FST_arl))   
MSE = tail(Straight_RF$mse ,1 )  
MSEeq = mean ((predict(Straight_RF, StraightMeanDF.train) - StraightMeanDF.train$FST_arl)^2) 
MSE2 = mean ((predict(Straight_RF, StraightMeanDF.valid) - StraightMeanDF.valid$FST_arl)^2)
Cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_arl)
Cor2 = cor(predict(Straight_RF, StraightMeanDF.valid), StraightMeanDF.valid$FST_arl)


#Add straight line parameters to the vectors
RSQ_vec	        = c(RSQ)
RSQeq_vec = c(RSQeq)
MSE_vec   = c(MSE)
MSEeq_vec = c(MSEeq)
MSE2_vec  = c(MSE2)
Cor1_vec  = c(Cor1)
Cor2_vec  = c(Cor2)


fit = lm(Straight_RF$predicted ~ StraightMeanDF.train$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run5_StraightRF_TrainingScatter.pdf", 5, 5)
plot(StraightMeanDF.train$FST_arl, Straight_RF$predicted,  xlab ="Observed FST (training)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(Straight_RF, StraightMeanDF.valid) ~ StraightMeanDF.valid$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run5_StraightRF_ValidScatter.pdf", 5, 5)
plot(StraightMeanDF.valid$FST_arl, predict(Straight_RF, StraightMeanDF.valid),  xlab ="Observed FST (validation)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

StraightPred <- predict(env, Straight_RF)

	     	print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/FL_points_list_noKeys.csv", sep=",", header=T)
P.coordinates1 <- P.table[,c(3,2)]
P.points <- SpatialPoints(P.table[,c(3,2)])  # ... converted into a spatial object
proj4string(P.points) <- crs.geo  
#plot(P.points)

print("starting loops")

it <- 1
for (it in 1:10) {

  trFlorida <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix
  trFloridaC <- geoCorrection(trFlorida, type="c") 

  	                 print("transition matrix done")
  	                 rm(trFlorida)
  	                 gc()

  AtoT <- shortestPath(trFloridaC, P.points[1,], P.points[1,], output="SpatialLines")
  for (x in 1:10) {  
    for (y in (x+1):11) { 
     Ato <- shortestPath(trFloridaC, P.points[x,], P.points[y,], output="SpatialLines")
     AtoT <- AtoT + Ato
     #assign(paste0("Ato_px",x,"_py",y,"_it",it) , Ato)   

    }
  }
  AtoT = (AtoT[-1])

  LcpLoop <- raster::extract(env, AtoT, fun=mean, na.rm=TRUE)

  LcpLoopDF <- as.data.frame(LcpLoop)

  LcpLoopDF$FST_arl <- G.table$FST_arl


  #Break data 70/30 here

  LcpLoopDF.train = LcpLoopDF[TrainingPairs,]
  LcpLoopDF.valid = LcpLoopDF[-TrainingPairs,]


  LCP_RF = randomForest(FST_arl ~  arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf +  MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP, importance=TRUE, na.action=na.omit, data=LcpLoopDF.train)

  assign(paste0("LCP_RF", it), LCP_RF )

  			  print(paste0("finishing RF for iteration #", it))
  			  gc()
			  rm(trFloridaC)
			  gc()

RSQ = tail(LCP_RF$rsq ,1 )
RSQeq = 1 - (tail(LCP_RF$mse , 1) / var(LcpLoopDF.train$FST_arl))
MSE = tail(LCP_RF$mse ,1 )
MSEeq = mean ((predict(LCP_RF, LcpLoopDF.train) - LcpLoopDF.train$FST_arl)^2)
MSE2 = mean ((predict(LCP_RF, LcpLoopDF.valid) - LcpLoopDF.valid$FST_arl)^2)
Cor1 = cor(predict(LCP_RF, LcpLoopDF.train), LcpLoopDF.train$FST_arl)
Cor2 = cor(predict(LCP_RF, LcpLoopDF.valid), LcpLoopDF.valid$FST_arl)

RSQ_vec   = append(RSQ_vec, RSQ)
RSQeq_vec = append(RSQeq_vec, RSQeq)
MSE_vec   = append(MSE_vec, MSE)
MSEeq_vec = append(MSEeq_vec, MSEeq)
MSE2_vec  = append(MSE2_vec, MSE2)
Cor1_vec  = append(Cor1_vec, Cor1)
Cor2_vec  = append(Cor2_vec, Cor2)  


pred = predict(env, LCP_RF)

       	 	      rm(LCP_RF)

  assign(paste0("pred", it), pred)
  
  pred.cond <- 1/pred 

  	       	      rmr(pred)
 		      gc()
		      print(paste0("end of loop for iteration #", it))

}


save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run5.RData")


d = data.frame(RSQ = RSQ_vec, RSQeq = RSQeq_vec, MSE = MSE_vec, MSEeq = MSEeq_vec, MSE2 = MSE2_vec, Cor1 = Cor1_vec, Cor2=Cor2_vec) 
write.csv(d, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run5_ValidationTable.csv", row.names =FALSE)


sdfdsf

pos_max = which.max(RSQ_vec)

#debug this in R

if(pos_max > 1) {
best_it = pos_max - 1
RF = paste0("LCP_RF", best_it)
ResistanceMap = paste0("pred", best_it)
} else {
  print("The straight lines are the best model")
  best_it = pos_max - 1
  RF = Straight_RF
  ResistanceMap = StraightPred
}

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_BestPred_Run2.pdf", 5, 5)
plot(ResistanceMap)
dev.off()

fit = lm(RF$predicted ~ LcpLoopDF.train$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run3_BestModelScatter.pdf", 5, 5)
plot(LcpLoopDF.train$FST_arl, RF$predicted,  xlab ="Observed FST (train)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(RF, LcpLoopDF.valid) ~ LcpLoopDF.valid$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run2_BestModelScatter.pdf", 5, 5)
plot(LcpLoopDF.valid$FST_arl, predict(RF, LcpLoopDF.valid),  xlab ="Observed FST (valid)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/FL/SplitData_Run2_ImpVars.pdf", 5, 5)
varImpPlot(RF)
dev.off()