#Building iterative RF model with SW2 (11 points)

run=6

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
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SW2_FST_list.csv", sep=",", header=T) 

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

aridI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/SW_clip/AI_annual_SWclip.tif")
arid = aridI*1
proj4string(arid) <- crs.geo

accessI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/SW_clip/accessibility_to_cities_2015_v1.0_res_SWclip.tif")
access <- accessI*1
proj4string(access) <- crs.geo

precI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/SW_clip/bio12_mean_SWclip.tif")
prec <- precI*1
proj4string(prec) <- crs.geo

mean.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/SW_clip/bio1_mean_SWclip.tif")
mean.temp <- mean.tempI*1
proj4string(mean.temp) <- crs.geo

human.densityI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/SW_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_SWclip.tif")
human.density <- human.densityI*1
proj4string(human.density) <- crs.geo

frictionI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/SW_clip/friction_surface_2015_v1.0_SWclip.tif")
friction <- frictionI*1
proj4string(friction) <- crs.geo

min.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/SW_clip/bio6_mean_SWclip.tif")
min.temp <- min.tempI*1
proj4string(min.temp) <- crs.geo

NeedleleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_1_SWclip.tif")
Needleleaf = NeedleleafI*1
proj4string(Needleleaf) <- crs.geo

EvBroadleafI  = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_2_SWclip.tif")
EvBroadleaf = EvBroadleafI*1
proj4string(EvBroadleaf) <- crs.geo

DecBroadleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_3_SWclip.tif")
DecBroadleaf = DecBroadleafI*1
proj4string(DecBroadleaf) <- crs.geo

MiscTreesI =  raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_4_SWclip.tif")
MiscTrees = MiscTreesI*1
proj4string(MiscTrees) <- crs.geo

ShrubsI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_5_SWclip.tif")
Shrubs = ShrubsI*1
proj4string(Shrubs) <- crs.geo

HerbI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_6_SWclip.tif")
Herb = HerbI*1
proj4string(Herb) <- crs.geo

CropI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_7_SWclip.tif")
Crop <- CropI*1
proj4string(Crop) <- crs.geo

FloodI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_8_SWclip.tif")
Flood = FloodI*1
proj4string(Flood) <- crs.geo

UrbanI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_9_SWclip.tif")
Urban <- UrbanI*1
proj4string(Urban) <- crs.geo

SnowI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_10_SWclip.tif")
Snow = SnowI*1
proj4string(Snow) <- crs.geo

BarrenI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_11_SWclip.tif")
Barren = BarrenI*1
proj4string(Barren) <- crs.geo

WaterI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/SW_clip/consensus_full_class_12_SWclip.tif")
Water = WaterI*1
proj4string(Water) <- crs.geo

SlopeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/SW_clip/slope_1KMmedian_MERIT_SWclip.tif")
Slope = SlopeI*1
proj4string(Slope) <- crs.geo

AltitudeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/SW_clip/altitude_1KMmedian_MERIT_SWclip.tif")
Altitude = AltitudeI*1
proj4string(Altitude) <- crs.geo

PETI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET/SW_clip/pet_mean_SWclip.tif")
PET = PETI*1
proj4string(PET) <- crs.geo

DailyTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio2/SW_clip/bio2_mean_SWclip.tif")
DailyTempRange = DailyTempRangeI*1
proj4string(DailyTempRange) <- crs.geo

max.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5/SW_clip/bio5_mean_SWclip.tif")
max.temp = max.tempI*1
proj4string(max.temp) <- crs.geo

AnnualTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio7/SW_clip/bio7_mean_SWclip.tif")
AnnualTempRange = AnnualTempRangeI*1
proj4string(AnnualTempRange) <- crs.geo

prec.wetI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio13/SW_clip/bio13_mean_SWclip.tif")
prec.wet = prec.wetI*1
proj4string(prec.wet) <- crs.geo

prec.dryI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio14/SW_clip/bio14_mean_SWclip.tif")
prec.dry = prec.dryI*1
proj4string(prec.dry) <- crs.geo

GPPI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly_mean/SW_clip/GPP_mean_SWclip.tif")
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

#save.image("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/sc14_rasterstack_image.RData")


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
MSE_vec	      = c()
MSEeq_vec = c()
MSE2_vec  = c()
Cor1_vec  = c()
Cor2_vec  = c()

#Performance measures
RSQ = tail(Straight_RF$rsq ,1 ) 
RSQeq = 1 - (tail(Straight_RF$mse , 1) / var(StraightMeanDF$FST_arl))   
MSE = tail(Straight_RF$mse ,1 )  
MSEeq = mean ((predict(Straight_RF, StraightMeanDF) - StraightMeanDF$FST_arl)^2) 
MSE2 = mean ((predict(Straight_RF, StraightMeanDF.valid) - StraightMeanDF.valid$FST_arl)^2)
Cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_arl)
Cor2 = cor(predict(Straight_RF, StraightMeanDF.valid), StraightMeanDF.valid$FST_arl)

#Add straight line parameters to the vectors
RSQ_vec	              = c(RSQ)
RSQeq_vec = c(RSQeq)
MSE_vec   = c(MSE)
MSEeq_vec = c(MSEeq)
MSE2_vec  = c(MSE2)
Cor1_vec  = c(Cor1)
Cor2_vec  = c(Cor2)

fit = lm(Straight_RF$predicted ~ StraightMeanDF.train$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_StraightRF_TrainingScatter.pdf"), 5, 5)
plot(StraightMeanDF.train$FST_arl, Straight_RF$predicted,  xlab ="Observed FST (training)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(Straight_RF, StraightMeanDF.valid) ~ StraightMeanDF.valid$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_StraightRF_ValidScatter.pdf"), 5, 5)
plot(StraightMeanDF.valid$FST_arl, predict(Straight_RF, StraightMeanDF.valid),  xlab ="Observed FST (validation)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()


StraightPred <- predict(env, Straight_RF)

print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SW2_points_list.csv", sep=",", header=T)
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
  RSQeq = 1 - (tail(LCP_RF$mse , 1) / var(LcpLoopDF$FST_arl))
  MSE = tail(LCP_RF$mse ,1 )
  MSEeq = mean ((predict(LCP_RF, LcpLoopDF) - LcpLoopDF$FST_arl)^2)
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

save.image(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,".RData"))

d = data.frame(RSQ = RSQ_vec, RSQeq = RSQeq_vec, MSE = MSE_vec, MSEeq = MSEeq_vec, MSE2 = MSE2_vec, Cor1 = Cor1_vec, Cor2=Cor2_vec) 
write.csv(d, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_ValidationTable.csv"), row.names =FALSE)


RF0 = Straight_RF
RF1 = LCP_RF1 
RF2 = LCP_RF2 
RF3 = LCP_RF3 
RF4 = LCP_RF4 
RF5 = LCP_RF5 
RF6 = LCP_RF6 
RF7 = LCP_RF7
RF8 = LCP_RF8
RF9 = LCP_RF9
RF10 = LCP_RF10
resist0 = StraightPred
resist1 = pred1 
resist2 = pred2 
resist3 = pred3 
resist4 = pred4 
resist5 = pred5 
resist6 = pred6 
resist7 = pred7
resist8 = pred8
resist9 = pred9
resist10 = pred10

#Best iteration based on Cor2
pos_max = which.max(Cor2_vec)

best_it = pos_max - 1
RF = paste0("RF", best_it)
ResistanceMap = paste0("resist", best_it)


pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestCor2_Pred_it",best_it,".pdf"), 5, 5)
plot(get(ResistanceMap))
dev.off()

fit = lm(get(RF)$predicted ~ LcpLoopDF.train$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestCor2_TrainingScatter_it", best_it, ".pdf"), 5, 5)
plot(LcpLoopDF.train$FST_arl,get(RF)$predicted,  xlab ="Observed FST (train)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(get(RF), LcpLoopDF.valid) ~ LcpLoopDF.valid$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestCor2_ValidScatter_it", best_it,".pdf"), 5, 5)
plot(LcpLoopDF.valid$FST_arl, predict(get(RF), LcpLoopDF.valid),  xlab ="Observed FST (valid)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestCor2_ImpVars_it",best_it,".pdf"), 5, 5)
varImpPlot(get(RF))
dev.off()


#Best iteration based on RSQ
pos_max = which.max(RSQ_vec)

best_it = pos_max - 1
RF = paste0("RF", best_it)
ResistanceMap = paste0("resist", best_it)


pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestRSQ_Pred_it",best_it,".pdf"), 5, 5)
plot(get(ResistanceMap))
dev.off()

fit = lm(get(RF)$predicted ~ LcpLoopDF.train$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestRSQ_TrainingScatter_it", best_it, ".pdf"), 5, 5)
plot(LcpLoopDF.train$FST_arl,get(RF)$predicted,  xlab ="Observed FST (train)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(get(RF), LcpLoopDF.valid) ~ LcpLoopDF.valid$FST_arl)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestRSQ_ValidScatter_it", best_it,".pdf"), 5, 5)
plot(LcpLoopDF.valid$FST_arl, predict(get(RF), LcpLoopDF.valid),  xlab ="Observed FST (valid)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/SW2/SplitData_Run",run,"_BestRSQ_ImpVars_it",best_it,".pdf"), 5, 5)
varImpPlot(get(RF))
dev.off()