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

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

#load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_pt1_withKernel.RData")

Full.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/FST_list_NAmRF3.csv", sep=",", header=T)

#upload the resistance surface

mean_surfaceI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/CSE_best_preds/CSEData_mean.tif")
mean_surface = mean_surfaceI*1
proj4string(mean_surface) <- crs.geo

std_surfaceI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/CSE_best_preds/CSEData_stdev.tif")
std_surface = std_surfaceI*1
proj4string(std_surface) <- crs.geo

#convert it to transition matrix
  
trMean <- transition(mean_surface, transitionFunction=mean, directions=8) #make transitional matrix
trMeanC <- geoCorrection(trMean, type="c")


P.points1 <- SpatialPoints(Full.table[,c(7,6)])
P.points2 <- SpatialPoints(Full.table[,c(9,8)])
proj4string(P.points1) <- crs.geo
proj4string(P.points2) <- crs.geo
NumPairs <- length(P.points1)


#get parallelization set up
nw <- detectCores()
# cl <- makePSOCKcluster(nw) # is create multiple copy and it is usefull for works in multiple node
# registerDoParallel(cl)     # is create multiple copy and it is usefull for works in multiple node
registerDoMC(cores=nw)       # is create forks of the data good; for one node many cpu

print("cores registerred")


#To save the LCP lines

Ato.all <- foreach(r=1:NumPairs, .combine='+', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
  shortestPath(trMean, P.points1[r], P.points2[r]  , output="SpatialLines")
}

print("Ato.all complete")


LcpLoop <- foreach(r=1:NumPairs, .combine='rbind', .packages=c('raster', 'gdistance')  ,   .inorder=TRUE   ) %dopar% {
  raster::extract(mean_surface,  Ato.all[r]     , fun=mean, na.rm=TRUE)
}

print("LcpLoop loop complete")


save.image("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_validation_CSE.RData")

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/CSE_best_preds/CSEData_mean.pdf", 5, 5)
plot(mean_surface)
dev.off()

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/CSE_best_preds/CSEData_stdev.pdf", 5, 5)
plot(std_surface)
dev.off()