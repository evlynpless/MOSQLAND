library("sp")
library("spatstat")
library("maptools")
library("raster")
library("SDraw")
library("tidyverse")
library("foreach")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

mean_surfaceI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/ResidFST_best_preds/ResidData_mean.tif")
mean_surface = mean_surfaceI*1
proj4string(mean_surface) <- crs.geo

adjusted_mean_surface = mean_surface - 0.064

#make pdf

