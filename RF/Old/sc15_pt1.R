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

load("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/rasterstack_image_withKernel.RData")

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

#Randomly shuffle the data
yourData<-G.table[sample(nrow(G.table)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  write.csv(testData, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/testData_", i, ".csv"))
  trainData <- yourData[-testIndexes, ]
  assign(paste0("trainData_", i), trainData)
  write.csv(trainData, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/trainData_", i, ".csv"))
}


save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/sc15_pt1.RData")
