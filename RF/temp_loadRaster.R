library("sp")
library("spatstat")
library("maptools")
library("raster")


setwd("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15")


load("sc15_pt2_withKernel1.RData")

print("first RData loaded")


writeRaster(pred1,filename="/Fold1_pred1.tif",options=c("COMPRESS=DEFLATE "),formats=GTiff,overwrite=TRUE)

