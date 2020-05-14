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

std_surfaceI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/ResidFST_best_preds/ResidData_stdev.tif")
std_surface = std_surfaceI*1
proj4string(std_surface) <- crs.geo

save.image("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/temp_makePDF_ResidData.RData")

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/ResidFST_best_preds/ResidData_mean.pdf", 5, 5)
plot(mean_surface)
dev.off()

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/ResidFST_best_preds/ResidData_std.pdf", 5, 5)
plot(std_surface)
dev.off()