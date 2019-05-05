#Building iterative RF model with Florida as test case
# source("/gpfs/loomis/home.grace/fas/powell/esp38/scripts/MOSQLAND/RF/sc08D_iterativeRF_NAm1_MiddleStart.R")

gc()
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

nw <- detectCores()
cl <- makePSOCKcluster(nw)
registerDoParallel(cl)

#crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

print("library done")

###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

#T1 = Sys.time()
#T1

#G.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/FST_list_NAmRF1_noKeys.csv", sep=",", header=T) 

#create dataframes of begin and end coordinates from a file:
#begin.table <- G.table[,c(4,3)]
#begin.coord <- begin.table
#coordinates(begin.coord) <- c("long1", "lat1")

#end.table <- G.table[,c(6,5)]
#end.coord <- end.table
#coordinates(end.coord) <- c("long2", "lat2")

#p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

#spatial.p <- as(p, "SpatialLines")
#proj4string(spatial.p) <- crs.geo  # define projection system of our data

#print("G table done")
###############################################
#Create raster stack 
###############################################

#aridI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2.tif")
#arid = aridI*1
#proj4string(arid) <- crs.geo

#accessI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip2.tif")
#access = accessI*1
#proj4string(access) <- crs.geo

#precI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip2.tif")
#prec = precI*1
#proj4string(prec) <- crs.geo

#mean.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/NAm_clip/bio1_NAmClip2.tif")
#mean.temp = mean.tempI*1
#proj4string(mean.temp) <- crs.geo

#human.densityI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif")
#human.density = human.densityI*1
#proj4string(human.density) <- crs.geo

#cropI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip2.tif")
#crop = cropI*1
#proj4string(crop) <- crs.geo

#urbanI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip2.tif")
#urban = urbanI*1
#proj4string(urban) <- crs.geo

#frictionI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip2.tif")
#friction = frictionI*1
#proj4string(friction) <- crs.geo

#min.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip2.tif")
#min.temp = min.tempI*1
#proj4string(min.temp) <- crs.geo

#ABSHUMI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/NAm_clip/ABS50_res_cubicSP_NAmClip2.tif")
#ABSHUM = ABSHUMI*1
#proj4string(ABSHUM) <- crs.geo

#env=stack(arid, access, prec, mean.temp, human.density, crop, urban, friction, min.temp, ABSHUM)

#names(env) [1] <- "arid"
#names(env) [2] <- "access"
#names(env) [3] <- "prec"
#names(env) [4] <- "mean.temp"
#names(env) [5] <- "human.density"
#names(env) [6] <- "crop"
#names(env) [7] <- "urban"
#names(env) [8] <- "friction"
#names(env) [9] <- "min.temp"
#names(env) [10] <- "ABSHUM"

#print("raster stack done")

#T2 = Sys.time()
#T2

#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_NAmRasterStack.RData")

#load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_NAmRasterStack.RData")

########################################
#Calculate mean of straight lines and making initial RF model
#######################################


nw <- detectCores()
cl <- makePSOCKcluster(nw) # is create multiple copy and it is usefull for works in multiple node
registerDoParallel(cl)     # is create multiple copy and it is usefull for works in multiple node
registerDoMC(cores=nw)       # is create forks of the data good; for one node many cpu


#StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

#print("First extraction done")

#T3 = Sys.time()
#T3

#StraightMeanDF <- as.data.frame(StraightMean)

#StraightMeanDF$FST_arl <- G.table$FST_arl

  
#Straight_RF = randomForest(FST_arl ~   arid +   access + prec  +   mean.temp  +   human.density  +   crop   +    urban  +   friction + min.temp + ABSHUM, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

#print("Straight_RF done")

#T4 = Sys.time()
#T4

#StraightPred <- predict(env, Straight_RF)

#T5 = Sys.time()
#T5

#pred.cond <- 1/StraightPred

#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08A_StraightRFDone.RData")


#Prepare points for use in least cost path loops
#P.table <- read.table(file="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/NAmRF1_points_list_noKeys.csv", sep=",", header=T)
#P.coordinates1 <- P.table[,c(3,2)]
#P.points <- SpatialPoints(P.table[,c(3,2)])  # ... converted into a spatial object
#proj4string(P.points) <- crs.geo  
#plot(P.points)

#NumSites <- length(P.points)

#print("starting iterations")

#T6 = Sys.time()
#T6



#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_BeforeLoop.RData")

#load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_BeforeLoop.RData")

#print("image loaded")

#trNAm1 <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix

#trNAm1C <- geoCorrection(trNAm1, type="c")

#save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_BeforeLoop2.RData")

load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_BeforeLoop2.RData")


P<-foreach(i=0:3) %dopar%  i

# create list for iteration 
a=c()
for (x in 1:20) { for (y in 1:20) {a = rbind (a,   print(c(x,y)) )  } } 
FT=a[,1] != a[,2]
pointlist=a[ which(FT),]  # point list is a vector that identify the index in the P.points


it <- 1
for (it in 1:1) {
  
  #trNAm1 <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix

  #print("transition matrix done")

  #trNAm1C <- geoCorrection(trNAm1, type="c") 

LcpLoop <- foreach(r=1:8 , .combine='rbind', .packages=c('raster', 'gdistance')  , .verbose=TRUE  ,   .inorder=TRUE   ) %dopar% {
	Ato <- shortestPath(trNAm1C, P.points[pointlist[r,1]] ,P.points[pointlist[r,2]]  , output="SpatialLines")
        cbind (  pointlist[r,1] ,  pointlist[r,2]  , raster::extract(env,  Ato     , fun=mean, na.rm=TRUE))

}                 
}
               

#write.csv(LcpLoop, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc08D_LcpLoop2.csv")

zdfgzgsdfg

  LcpLoop = (LcpLoop[-1])

  T7 = Sys.time()
  T7

  LcpLoopDF <- as.data.frame(LcpLoop)

  LcpLoopDF$FST_arl <- G.table$FST_arl

  LCP_RF = randomForest(FST_arl ~  arid +   access + prec  +   mean.temp  +   human.density  +   crop   +    urban  +   friction + min.temp + ABSHUM, importance=TRUE, na.action=na.omit, data=LcpLoopDF)

  pred = predict(env, LCP_RF)

  assign(paste0("pred", it), pred)

  T8 = Sys.time()
  T8

  pred.cond <- 1/pred
  
  print("round done")
# }

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_1/sc07D_sc08D_MiddleStart.RData")


#T2-T1
#T3-T2
#T4-T3
#T5-T4
#T6-T5
#T7-T6
#T8-T7
