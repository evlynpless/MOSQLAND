run=1

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

#load prepared raster stack for North America

load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/rasterstack_image.RData")

rmr=function(x){
## function to truly delete raster and temporary files associated with them
if(class(x)=="RasterLayer"&grepl("^/tmp",x@file@name)&fromDisk(x)==T){
file.remove(x@file@name,sub("grd","gri",x@file@name))
rm(x)
}
}

GeoDist <- raster(nrows = 1500, ncols = 4140, xmn = -113.5, xmx = -79, ymn = 24, ymx = 36.5, crs = crs.geo, vals = 1)
names(GeoDist) <- c('GeoDist')

rm(env)

env=stack(arid, access, prec, mean.temp, human.density, friction, min.temp, Needleleaf, EvBroadleaf, DecBroadleaf, MiscTrees, Shrubs, Herb, Crop, Flood, Urban, Snow, Barren, Water, Slope, Altitude, PET, DailyTempRange, max.temp, AnnualTempRange, prec.wet, prec.dry, GPP, GeoDist)


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
names(env) [29] <- "GeoDist"

###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/FST_list_NAmRF3.csv", sep=",", header=T)


#create dataframes of begin and end coordinates from a file:
begin.table <- G.table[,c(7,6)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- G.table[,c(9,8)]
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

DistVar <- raster::extract(GeoDist, spatial.p, fun=sum, na.rm=TRUE)

StraightMeanDF <- as.data.frame(StraightMean)

StraightMeanDF$GeoDist <- DistVar

#any way to combine using join instead?
StraightMeanDF$FST_lin <- G.table$FST_lin

#option of trying DPS
#StraightMeanDF$DPS <- G.table$DPS
  
#Creating validation and testing datasets
NumPairs = nrow(StraightMeanDF)
Training = NumPairs * 0.7
TrainingInt = round(Training)
TrainingPairs = sample(1:NumPairs, TrainingInt, replace = FALSE)
StraightMeanDF.train = StraightMeanDF[TrainingPairs,]
StraightMeanDF.valid = StraightMeanDF[-TrainingPairs,]

set.seed(NULL)

#Tuning RF for 
tune_x <- StraightMeanDF.train[,1:29]
tune_y <- StraightMeanDF.train[,29]
bestmtry <- tuneRF(x, y, stepFactor=1.5, improve=1e-5, ntree=500)
mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

Straight_RF = randomForest(FST_lin ~   arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + GeoDist, importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=StraightMeanDF.train)

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
RMSE2 = sqrt(((predict(Straight_RF, StraightMeanDF.valid) - StraightMeanDF.valid$FST_lin)^2))
MAE = mean(abs(Straight_RF$predicted - StraightMeanDF.train$FST_lin))
MAE2 =  mean(abs(predict(Straight_RF, StraightMeanDF.train) - StraightMeanDF.train$FST_lin))
MAE3 = mean(abs(predict(Straight_RF, StraightMeanDF.valid) - StraightMeanDF.valid$FST_lin))
Cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_lin)
Cor2 = cor(predict(Straight_RF, StraightMeanDF.valid), StraightMeanDF.valid$FST_lin)

#Add straight line parameters to the vectors
RSQ_vec	  = c(RSQ)
RMSE_vec   = c(RMSE)
RMSE2_vec  = c(RMSE2)
MAE_vec = c(MAE)
MAE2_vec = c(MAE2)
MAE3_vec = c(MAE3)
Cor1_vec  = c(Cor1)
Cor2_vec  = c(Cor2)

fit = lm(Straight_RF$predicted ~ StraightMeanDF.train$FST_lin)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_StraightRF_TrainingScatter.pdf"), 5, 5)
plot(StraightMeanDF.train$FST_lin, Straight_RF$predicted,  xlab ="Observed FST (training)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(Straight_RF, StraightMeanDF.valid) ~ StraightMeanDF.valid$FST_lin)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_StraightRF_ValidScatter.pdf"), 5, 5)
plot(StraightMeanDF.valid$FST_lin, predict(Straight_RF, StraightMeanDF.valid),  xlab ="Observed FST (validation)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()


StraightPred <- predict(env, Straight_RF)


print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/NAmRF3_points_list.csv", sep=",", header=T)
P.coordinates1 <- G.table[,c(3,2)]
P.points <- SpatialPoints(P.table[,c(3,2)])  # ... converted into a spatial object
proj4string(P.points) <- crs.geo  

NumPoints = length(P.points)

#get parallelization set up
nw <- detectCores()
# cl <- makePSOCKcluster(nw) # is create multiple copy and it is usefull for works in multiple node
# registerDoParallel(cl)     # is create multiple copy and it is usefull for works in multiple node
registerDoMC(cores=nw)       # is create forks of the data good; for one node many cpu

print("cores registerred")

# create list for iteration
a=c() 
for (x in 1:(NumPoints-1)) {      
  for (y in (x+1):NumPoints) {   
  a = rbind (a, c(x,y) )
}
}

FT=a[,1] != a[,2]
pointlist=a[ which(FT),]


print("starting loops")


it <- 1
for (it in 1:3) {
  
  rm(trNAm1C)
  gc()
  

  trNAm1 <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix

  print("transition matrix done")

  trNAm1C <- geoCorrection(trNAm1, type="c") 

  rm(trNAm1)
  gc()


  LcpLoop <- foreach(r=1:NumPairs, .combine='rbind', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
	Ato <- shortestPath(trNAm1C, P.points[pointlist[r,1]] ,P.points[pointlist[r,2]]  , output="SpatialLines")
        cbind (  pointlist[r,1] ,  pointlist[r,2]  , raster::extract(env,  Ato     , fun=mean, na.rm=TRUE))

}

	LcpLoopDF <- as.data.frame(LcpLoop)

	LcpLoopDF$GeoDist <- DistVar


	LcpLoopDF = merge(LcpLoopDF, G.table, by=c("V1", "V2"))

	#Break data 70/30 here

	LcpLoopDF.train = LcpLoopDF[TrainingPairs,]
	LcpLoopDF.valid = LcpLoopDF[-TrainingPairs,]

	tune_x <- LcpLoopDF.train[,1:28]
	tune_y <- LcpLoopDF.train[,29]
	bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
	mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

	LCP_RF = randomForest(FST_lin ~  arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf +  MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + GeoDist, importance=TRUE, mtry=mtry_opt, na.action=na.omit, data=LcpLoopDF.train)

assign(paste0("LCP_RF", it), LCP_RF )


print(paste0("finishing RF for iteration #", it))

gc()

rm(trNAm1C)

gc()

RSQ = tail(LCP_RF$rsq ,1 )
RMSE = sqrt(tail(LCP_RF$mse ,1 ))
RMSE2 = sqrt(((predict(LCP_RF, LcpLoopDF.valid) - LcpLoopDF.valid$FST_lin)^2))
MAE = mean(abs(LCP_RF$predicted - LcpLoopDF.train$FST_lin))
MAE2 = mean(abs(predict(LCP_RF, LcpLoopDF.train) - LcpLoopDF.train$FST_lin))
MAE3 = mean(abs(predict(LCP_RF, LcpLoopDF.valid) - LcpLoopDF.valid$FST_lin))
Cor1 = cor(predict(LCP_RF, LcpLoopDF.train), LcpLoopDF.train$FST_lin)
Cor2 = cor(predict(LCP_RF, LcpLoopDF.valid), LcpLoopDF.valid$FST_lin)

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

save.image(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_withGeoDistFull_Run",run,".RData"))


d = data.frame(RSQ = RSQ_vec, RMSE = RMSE_vec, RMSE2 = RMSE2_vec, MAE = MAE_vec, MAE2 = MAE2_vec, MAE3 = MAE3_vec, Cor1 = Cor1_vec,  Cor2 = Cor2_vec) 
write.csv(d, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run", run, "_ValidationTable.csv"), row.names =FALSE)

RF0 = Straight_RF
RF1 = LCP_RF1 
RF2 = LCP_RF2 
RF3 = LCP_RF3 
#RF4 = LCP_RF4 
#RF5 = LCP_RF5 
#RF6 = LCP_RF6 
#RF7 = LCP_RF7
#RF8 = LCP_RF8
#RF9 = LCP_RF9
#RF10 = LCP_RF10
resist0 = StraightPred
resist1 = pred1 
resist2 = pred2 
resist3 = pred3 
#resist4 = pred4 
#resist5 = pred5 
#resist6 = pred6 
#resist7 = pred7
#resist8 = pred8
#resist9 = pred9
#resist10 = pred10

#Best iteration based on Cor2
pos_max = which.max(Cor2_vec)

best_it = pos_max - 1
RF = paste0("RF", best_it)
ResistanceMap = paste0("resist", best_it)

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_BestCor2_Pred_it",best_it,".pdf"), 5, 5)
plot(get(ResistanceMap))
dev.off()

fit = lm(get(RF)$predicted ~ LcpLoopDF.train$FST_lin)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_BestCor2_TrainingScatter_it", best_it, ".pdf"), 5, 5)
plot(LcpLoopDF.train$FST_lin,get(RF)$predicted,  xlab ="Observed FST (train)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

fit = lm(predict(get(RF), LcpLoopDF.valid) ~ LcpLoopDF.valid$FST_lin)
adjr2 = round(summary(fit)$adj.r.squared, digits=3)
pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_BestCor2_ValidScatter_it", best_it,".pdf"), 5, 5)
plot(LcpLoopDF.valid$FST_lin, predict(get(RF), LcpLoopDF.valid),  xlab ="Observed FST (valid)", ylab="Predicted FST")
legend("bottomright", legend=c(paste0("Adj. R^2 = ", adjr2)), cex=0.7)
dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12D/LinFSTData_Run",run,"_BestCor2_ImpVars_it",best_it,".pdf"), 5, 5)
varImpPlot(get(RF))
dev.off()
