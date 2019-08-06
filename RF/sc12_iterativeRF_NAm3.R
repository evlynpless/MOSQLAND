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

StraightMeanDF <- as.data.frame(StraightMean)

#any way to combine using join instead?
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

gc()

#why aren't these the same value?
Straight_rsq = tail(Straight_RF$rsq ,1 ) 
#Straight_rsq_fromEquation = 1 - (tail(Straight_RF$mse , 1) / var(StraightMeanDF.train$FST_arl))     
write.table(Straight_rsq, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_RSQ_Table_", runnum, ".txt"), row.names =FALSE, col.names=FALSE)
#write.table(Straight_rsq_fromEquation, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_RSQ_fromEquation_Table.txt", row.names =FALSE, col.names=FALSE)

Straight_mse = tail(Straight_RF$mse ,1 )  
#Straight_mse_fromEquation = mean ((predict(Straight_RF, StraightMeanDF.train) - StraightMeanDF.train$FST_arl)^2) 
Straight_mse2 = mean ((predict(Straight_RF, StraightMeanDF.valid) - StraightMeanDF.valid$FST_arl)^2)
write.table(Straight_mse, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE_Table_", runnum, ".txt"), row.names =FALSE, col.names=FALSE)
#write.table(Straight_mse_fromEquation, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE_fromEquation_Table.txt", row.names =FALSE, col.names=FALSE)
write.table(Straight_mse2, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE2_Table_", runnum, ".txt"), row.names =FALSE, col.names=FALSE)

#save these as variables and save to a table
cor1 = cor(predict(Straight_RF, StraightMeanDF.train), StraightMeanDF.train$FST_arl)
cor2 = cor(predict(Straight_RF, StraightMeanDF.valid), StraightMeanDF.valid$FST_arl)
write.table(cor1, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_InternalValidation_", runnum, ".txt"), row.names =FALSE, col.names=FALSE)
write.table(cor2, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_ExternalValidation_", runnum, ".txt"), row.names =FALSE, col.names=FALSE)

#pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/Straight_scatter1.pdf", 5, 5)
#plot(StraightMeanDF.train$FST_arl, Straight_RF$predicted)
#dev.off()

#pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/Straight_scatter_test.pdf", 5, 5)
#plot(StraightMeanDF.train$FST_arl, predict(Straight_RF, StraightMeanDF.train))
#dev.off()

#pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/Straight_scatter2.pdf", 5, 5)
#plot(StraightMeanDF.valid$FST_arl, predict(Straight_RF, StraightMeanDF.valid))
#dev.off()

StraightPred <- predict(env, Straight_RF)

#write this as tif and remove it

print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

#Prepare points for use in least cost path loops
P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/NAmRF3_points_list.csv", sep=",", header=T)
P.coordinates1 <- P.table[,c(3,2)]
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

#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_BeforeLoops.RData")


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


	LcpLoopDF = merge(LcpLoopDF, G.table, by=c("V1", "V2"))

	#Break data 70/30 here

	LcpLoopDF.train = LcpLoopDF[TrainingPairs,]
	LcpLoopDF.valid = LcpLoopDF[-TrainingPairs,]

	LCP_RF = randomForest(FST_arl ~  arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf +  MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP, importance=TRUE, na.action=na.omit, data=LcpLoopDF.train)

assign(paste0("LCP_RF", it), LCP_RF )


print(paste0("finishing RF for iteration #", it))

gc()

rm(trNAm1C)

gc()


LCP_rsq = tail(LCP_RF$rsq ,1 )
#LCP_rsq_fromEquation = 1 - (tail(LCP_RF$mse , 1) / var(LcpLoopDF.train$FST_arl))
write.table(LCP_rsq, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_RSQ_Table_", runnum, ".txt"), append=TRUE, row.names =FALSE, col.names=FALSE)
#write.table(LCP_rsq_fromEquation, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_RSQ_fromEquation_Table.txt", append=TRUE, row.names =FALSE, col.names=FALSE)


LCP_mse = tail(LCP_RF$mse ,1 )
#LCP_mse_fromEquation = mean ((predict(LCP_RF, LcpLoopDF.train) - LcpLoopDF.train$FST_arl)^2)
LCP_mse2 = mean ((predict(LCP_RF, LcpLoopDF.valid) - LcpLoopDF.valid$FST_arl)^2)

write.table(LCP_mse, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE_Table_", runnum, ".txt"), append=TRUE, row.names =FALSE, col.names=FALSE)
write.table(LCP_mse2, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE2_Table_" runnum, ".txt", append=TRUE, row.names =FALSE, col.names=FALSE)
#write.table(LCP_mse_fromEquation, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_MSE_fromEquation_Table_", runnum, ".txt"), append=TRUE, row.names =FALSE, col.names=FALSE)


cor1 = cor(predict(LCP_RF, LcpLoopDF.train), LcpLoopDF.train$FST_arl)
cor2 = cor(predict(LCP_RF, LcpLoopDF.valid), LcpLoopDF.valid$FST_arl)

write.table(cor1, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_InternalValidation_", runnum, ".txt"), append=TRUE, row.names =FALSE, col.names=FALSE)
write.table(cor2, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_ExternalValidation_", runnum, ".txt"), append=TRUE, row.names =FALSE, col.names=FALSE)

#pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/LCP_scatter1_", it, ".pdf"), 5, 5)
#plot(LcpLoopDF.train$FST_arl, LCP_RF$predicted)
#dev.off()

#pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/LCP_scatter2_", it, ".pdf"), 5, 5)
#plot(LcpLoopDF.valid$FST_arl, predict(LCP_RF, LcpLoopDF.valid))
#dev.off()

pred = predict(env, LCP_RF)

print(paste0("finishing prediction for iteration #", it))


rm(LCP_RF)

assign(paste0("pred", it), pred)
  
  pred.cond <- 1/pred 
  
  rmr(pred)
  gc()

  print(paste0("end of loop for iteration #", it))


}                 

save.image(file = paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_", runnum, ".RData"))
