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

load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_pt1_withKernel.RData")


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


########################################
#Calculate mean of straight lines and making initial RF model
#######################################
#For training
StraightMean.train <- raster::extract(env, spatial.p.train, fun=mean, na.rm=TRUE)

StraightMeanDF.train <- as.data.frame(StraightMean.train)

StraightMeanDF.train$FST_lin <- Train.table$FST_lin

#For testing
StraightMean.test <- raster::extract(env, spatial.p.test, fun=mean, na.rm=TRUE)

StraightMeanDF.test <- as.data.frame(StraightMean.test)

StraightMeanDF.test$FST_lin <- Test.table$FST_lin

set.seed(NULL)

#check these
tune_x <- StraightMeanDF.train[,names(env)]
tune_y <- StraightMeanDF.train[,c("FST_lin")]
bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

Straight_RF = randomForest(FST_lin ~   arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=StraightMeanDF.train)

gc()

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
RMSE2 = sqrt(mean((predict(Straight_RF, StraightMeanDF.test) - StraightMeanDF.test$FST_lin)^2))
MAE = mean(abs(Straight_RF$predicted - StraightMeanDF.train$FST_lin))
MAE2 =  mean(abs(predict(Straight_RF, StraightMeanDF.train) - StraightMeanDF.train$FST_lin))
MAE3 = mean(abs(predict(Straight_RF, StraightMeanDF.test) - StraightMeanDF.test$FST_lin))
#Cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_lin)
Cor1 = cor(Straight_RF$predicted, StraightMeanDF.train$FST_lin)
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
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_StraightRF_TrainingScatter.pdf"), 5, 5)
plot(StraightMeanDF.train$FST_lin, Straight_RF$predicted,  xlab ="Observed FST* (training)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Cor1,3))), cex=0.7)
dev.off()

fit = lm(predict(Straight_RF, StraightMeanDF.test) ~ StraightMeanDF.test$FST_lin)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_StraightRF_ValidScatter.pdf"), 5, 5)
plot(StraightMeanDF.test$FST_lin, predict(Straight_RF, StraightMeanDF.test),  xlab ="Observed FST (testing)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Cor2,3))), cex=0.7)
dev.off()

StraightPred <- predict(env, Straight_RF)

print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface


#Prepare points for use in least cost path loops - Training
P.points1.train <- SpatialPoints(Train.table[,c(8,7)])
P.points2.train <- SpatialPoints(Train.table[,c(10,9)])
proj4string(P.points1.train) <- crs.geo
proj4string(P.points2.train) <- crs.geo
NumPairs.train <- length(P.points1.train)


#Prepare points for use in least cost path loops - Testing
P.points1.test <- SpatialPoints(Test.table[,c(8,7)])
P.points2.test <- SpatialPoints(Test.table[,c(10,9)])
proj4string(P.points1.test) <- crs.geo
proj4string(P.points2.test) <- crs.geo
NumPairs.test		    <- length(P.points1.test)


#get parallelization set up
nw <- detectCores()
# cl <- makePSOCKcluster(nw) # is create multiple copy and it is usefull for works in multiple node
# registerDoParallel(cl)     # is create multiple copy and it is usefull for works in multiple node
registerDoMC(cores=nw)       # is create forks of the data good; for one node many cpu

print("cores registerred")


print("starting loops")

it <- 1
for (it in 1:6) {
  
  rm(trNAm1C)
  gc()
  
  trNAm1 <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix

  print("transition matrix done")

  trNAm1C <- geoCorrection(trNAm1, type="c") 

  rm(trNAm1)
  gc()


#Extract mean value from LCP for training data


Ato.all.train <- foreach(r=1:NumPairs.train, .combine='+', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
  shortestPath(trNAm1C, P.points1.train[r], P.points2.train[r]  , output="SpatialLines")
}

print("Ato.all.train complete")

assign(paste0("Ato.all.train", it), Ato.all.train )

  LcpLoop.train <- foreach(r=1:NumPairs.train, .combine='rbind', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
      raster::extract(env,  Ato.all.train[r]     , fun=mean, na.rm=TRUE)
}

print("LcpLoop.train complete")



#Extract mean value from LCP for testing data

Ato.all.test <- foreach(r=1:NumPairs.test, .combine='+', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
  shortestPath(trNAm1C, P.points1.test[r], P.points2.test[r]  , output="SpatialLines")
}

print("Ato.all.test complete")

assign(paste0("Ato.all.test", it), Ato.all.test )

  LcpLoop.test <- foreach(r=1:NumPairs.test, .combine='rbind', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
  raster::extract(env,  Ato.all.test[r]     , fun=mean, na.rm=TRUE)
}

print("LcpLoop.test complete")



	LcpLoopDF.train <- as.data.frame(LcpLoop.train)
	LcpLoopDF.train$FST_lin <- Train.table$FST_lin

        LcpLoopDF.test <- as.data.frame(LcpLoop.test)
        LcpLoopDF.test$FST_lin = Test.table$FST_lin

	tune_x <- LcpLoopDF.train[,names(env)]
	tune_y <- LcpLoopDF.train[,c("FST_lin")]
	bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
	mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

LCP_RF = randomForest(FST_lin ~  arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf +  MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, importance=TRUE, mtry=mtry_opt, na.action=na.omit, data=LcpLoopDF.train)

assign(paste0("LCP_RF", it), LCP_RF )


print(paste0("finishing RF for iteration #", it))

gc()

rm(trNAm1C)

gc()

#add validation parameters here
RSQ = tail(LCP_RF$rsq ,1 )
RMSE = sqrt(tail(LCP_RF$mse ,1 ))
RMSE2 = sqrt(mean((predict(LCP_RF, LcpLoopDF.test) - LcpLoopDF.test$FST_lin)^2))
MAE = mean(abs(LCP_RF$predicted - LcpLoopDF.train$FST_lin))
MAE2 =  mean(abs(predict(LCP_RF, LcpLoopDF.train) - LcpLoopDF.train$FST_lin))
MAE3 = mean(abs(predict(LCP_RF, LcpLoopDF.test) - LcpLoopDF.test$FST_lin))
#Cor1 = cor(predict(LCP_RF, LcpLoopDF.train), LcpLoopDF.train$FST_lin)
Cor1 = cor(LCP_RF$predicted, LcpLoopDF.train$FST_lin)
Cor2 = cor(predict(LCP_RF, LcpLoopDF.test), LcpLoopDF.test$FST_lin)


RSQ_vec   = append(RSQ_vec, RSQ)
RMSE_vec   = append(RMSE_vec, RMSE)
RMSE2_vec  = append(RMSE2_vec, RMSE2)
MAE_vec   = append(MAE_vec, MAE)
MAE2_vec  = append(MAE2_vec, MAE2)
MAE3_vec  = append(MAE3_vec, MAE3)
Cor1_vec  = append(Cor1_vec, Cor1)
Cor2_vec  = append(Cor2_vec, Cor2)


pred = predict(env, LCP_RF)

print(paste0("finishing prediction for iteration #", it))


rm(LCP_RF)

assign(paste0("pred", it), pred)
  
pred.cond <- 1/pred 
  
rmr(pred)
  
gc()

print(paste0("end of loop for iteration #", it))

}  

d = data.frame(RSQ = RSQ_vec, RMSE = RMSE_vec, RMSE2 = RMSE2_vec, MAE = MAE_vec, MAE2 = MAE2_vec, MAE3 = MAE3_vec, Cor1 = Cor1_vec,  Cor2 = Cor2_vec) 
write.csv(d, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold", foldnum, "_ValidationTable.csv"), row.names =FALSE)

RF0 = Straight_RF
RF1 = LCP_RF1 
RF2 = LCP_RF2 
RF3 = LCP_RF3 
RF4 = LCP_RF4 
RF5 = LCP_RF5 
RF6 = LCP_RF6 
#RF7 = LCP_RF7
#RF8 = LCP_RF8
#RF9 = LCP_RF9
#RF10 = LCP_RF10
resist0 = StraightPred
resist1 = pred1 
resist2 = pred2 
resist3 = pred3 
resist4 = pred4 
resist5 = pred5 
resist6 = pred6 
#resist7 = pred7
#resist8 = pred8
#resist9 = pred9
#resist10 = pred10

#Best iteration based on Cor2
pos_max = which.max(Cor2_vec)

best_it = pos_max - 1
RF = paste0("RF", best_it)
ResistanceMap = paste0("resist", best_it)

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_BestCor2_Pred_it",best_it,".pdf"), 5, 5)
plot(get(ResistanceMap))
dev.off()

fit = lm(get(RF)$predicted ~ LcpLoopDF.train$FST_lin)
#adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_BestCor2_TrainingScatter_it", best_it, ".pdf"), 5, 5)
plot(LcpLoopDF.train$FST_lin,get(RF)$predicted,  xlab ="Observed FST* (train)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Cor1,3))), cex=0.7)
dev.off()

fit = lm(predict(get(RF), LcpLoopDF.test) ~ LcpLoopDF.test$FST_lin)
#adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_BestCor2_ValidScatter_it", best_it,".pdf"), 5, 5)
plot(LcpLoopDF.test$FST_lin, predict(get(RF), LcpLoopDF.test),  xlab ="Observed FST* (test)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Cor2,3))), cex=0.7)
dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFSTData_Fold",foldnum,"_BestCor2_ImpVars_it",best_it,".pdf"), 5, 5)
varImpPlot(get(RF))
dev.off()

save.image(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_pt2_withKernel",foldnum,".RData"))
