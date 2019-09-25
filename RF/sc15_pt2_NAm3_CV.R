#Import foldnum for 10-fold cross validation

foldnum<-Sys.getenv(c('foldnum'))
print(foldnum)


#Import packages
library("sp")
library("spatstat")
library("maptools")
library("raster")
library("randomForest")
library("gdistance")
library("SDraw")
library("tidyverse")
library("foreach")
library("doParallel")
library("doMC")

load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_pt1.RData")


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

#need to download test

Test.table <- read.table(file=paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/testData_", foldnum, ".csv"), sep=",", header=T)

Train.table <- read.table(file=paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/trainData_", foldnum, ".csv"), sep=",", header=T)

#For train data
#create dataframes of begin and end coordinates from a file:
begin.table <- Train.table[,c(8,7)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- Train.table[,c(10,9)]
end.coord <- end.table
coordinates(end.coord) <- c("long2", "lat2")

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

spatial.p.train <- as(p, "SpatialLines")
proj4string(spatial.p.train) <- crs.geo  # define projection system of our data

print("spatial lines done")


#For test data
begin.table <- Test.table[,c(8,7)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- Test.table[,c(10,9)]
end.coord <- end.table
coordinates(end.coord) <- c("long2", "lat2")

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

spatial.p.test <- as(p, "SpatialLines")
proj4string(spatial.p.test) <- crs.geo  # define projection system of our data

print("spatial lines done")

########################################
#Calculate mean of straight lines and making initial RF model
#######################################
StraightMean.train <- raster::extract(env, spatial.p.train, fun=mean, na.rm=TRUE)

DistVar.train <- raster::extract(GeoDist, spatial.p.train, fun=sum, na.rm=TRUE)

StraightMeanDF.train <- as.data.frame(StraightMean.train)

StraightMeanDF.train$GeoDist <- DistVar.train

StraightMeanDF.train$FST_lin <- Train.table$FST_lin

#For testing
StraightMean.test <- raster::extract(env, spatial.p.test, fun=mean, na.rm=TRUE)

DistVar.test <- raster::extract(GeoDist, spatial.p.test, fun=sum, na.rm=TRUE)

StraightMeanDF.test <- as.data.frame(StraightMean.test)

StraightMeanDF.test$GeoDist <- DistVar.test

StraightMeanDF.test$FST_lin <- Test.table$FST_lin



set.seed(NULL)

#Tuning RF for straight lines
tune_x <- StraightMeanDF.train[,1:29]
tune_y <- StraightMeanDF.train[,29]
bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

Straight_RF = randomForest(FST_lin ~   arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + GeoDist, importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=StraightMeanDF.train)

gc()

save.image(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/StopAfterStraightRF_Fold",foldnum,".RData"))


sdfsdfsd


#Define empty vectors
RSQ_vec = c()
RMSE_vec = c()
RMSE2_vec = c()
MAE_vec = c()
MAE2_vec = c()
MAE3_vec = c()
Cor1_vec  = c()
Cor2_vec  = c()

#Validation parameters
RSQ = tail(Straight_RF$rsq ,1 )
RMSE = sqrt(tail(Straight_RF$mse ,1 ))
RMSE2 = sqrt(((predict(Straight_RF, StraightMeanDF.test) - StraightMeanDF.test$FST_lin)^2))
MAE = mean(abs(Straight_RF$predicted - StraightMeanDF.train$FST_lin))
MAE2 =  mean(abs(predict(Straight_RF, StraightMeanDF.train) - StraightMeanDF.train$FST_lin))
MAE3 = mean(abs(predict(Straight_RF, StraightMeanDF.test) - StraightMeanDF.test$FST_lin))
Cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_lin)
Cor2 = cor(predict(Straight_RF, StraightMeanDF.test), StraightMeanDF.test$FST_lin)

#Add straight line parameters to the vectors
RSQ_vec   = c(RSQ)
RMSE_vec   = c(RMSE)
RMSE2_vec  = c(RMSE2)
MAE_vec = c(MAE)
MAE2_vec = c(MAE2)
MAE3_vec = c(MAE3)
Cor1_vec  = c(Cor1)
Cor2_vec  = c(Cor2)

fit = lm(Straight_RF$predicted ~ StraightMeanDF.train$FST_lin)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Run",run,"_StraightRF_TrainingScatter.pdf"), 5, 5)
plot(StraightMeanDF.train$FST_lin, Straight_RF$predicted,  xlab ="Observed FST (training)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", Cor2)), cex=0.7)
dev.off()

fit = lm(predict(Straight_RF, StraightMeanDF.test) ~ StraightMeanDF.test$FST_lin)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Run",run,"_StraightRF_ValidScatter.pdf"), 5, 5)
plot(StraightMeanDF.valid$FST_lin, predict(Straight_RF, StraightMeanDF.test),  xlab ="Observed FST (testing)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", Cor2)), cex=0.7)
dev.off()


StraightPred <- predict(env, Straight_RF)


print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Check if it's working before adding LCP

save.image(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12E/LinFSTData_withGeoDistFull_Fold",foldnum,".RData"))
