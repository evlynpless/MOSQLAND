#!/bin/bash                                        
#SBATCH -p day                                     
#SBATCH --mem=100g                                
#SBATCH -t 24:00:00                                
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/rasterstack.sh.%J.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/rasterstack.sh.%J.err
#SBATCH --mail-type=ALL                            
#SBATCH --mail-user=evlyn.pless@yale.edu           
#SBATCH --job-name=rasterstack.sh                                                                                        
# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF_EAfr/rasterstack.sh                                                     
ulimit -c 0 

module load R/3.5.3-foss-2018a-X11-20180131
module load GDAL/3.1.0-foss-2018a-Python-3.6.4
module load Rpkgs/RGDAL/1.2-5 #not working
module load Rpkgs/DOPARALLEL/1.0.3 #not working


R --vanilla -no-readline -q  << 'EOF'

#Purpose: creating North America raster stack for use in iterative RF model

#Download packages
library("sp")
#library("spatstat")
#library("maptools")
library("raster")
library("doParallel")
library("rgdal")
#library("gdistance")
#library("SDraw")
#library("tidyverse")

#Add coordinate system
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

###############################################
#Create raster stack
###############################################

#Upload each raster and define its coordinate system
#Multiply by 1 is a trick to save memory later

aridI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/AI_annual_EAfrClip.tif")
arid = aridI*1
proj4string(arid) <- crs.geo

accessI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/accessibility_to_cities_2015_v1.0_res_EAfrClip.tif")
access <- accessI*1
proj4string(access) <- crs.geo

precI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio12_mean_EAfrClip.tif")
prec <- precI*1
proj4string(prec) <- crs.geo

mean.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio1_mean_EAfrClip.tif")
mean.temp <- mean.tempI*1
proj4string(mean.temp) <- crs.geo

human.densityI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_EAfrClip.tif")
human.density <- human.densityI*1
proj4string(human.density) <- crs.geo

min.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio6_mean_EAfrClip.tif")
min.temp <- min.tempI*1
proj4string(min.temp) <- crs.geo

SlopeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/slope_1KMmedian_MERIT_EAfrClip.tif")
Slope = SlopeI*1
proj4string(Slope) <- crs.geo

AltitudeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/altitude_1KMmedian_MERIT_EAfrClip.tif")
Altitude = AltitudeI*1
proj4string(Altitude) <- crs.geo

PETI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/pet_mean_EAfrClip.tif")
PET = PETI*1
proj4string(PET) <- crs.geo

DailyTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio2_mean_EAfrClip.tif")
DailyTempRange = DailyTempRangeI*1
proj4string(DailyTempRange) <- crs.geo

max.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio5_mean_EAfrClip.tif")
max.temp = max.tempI*1
proj4string(max.temp) <- crs.geo

AnnualTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio7_mean_EAfrClip.tif")
AnnualTempRange = AnnualTempRangeI*1
proj4string(AnnualTempRange) <- crs.geo

prec.wetI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio13_mean_EAfrClip.tif")
prec.wet = prec.wetI*1
proj4string(prec.wet) <- crs.geo

prec.dryI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/bio14_mean_EAfrClip.tif")
prec.dry = prec.dryI*1
proj4string(prec.dry) <- crs.geo

GPPI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/GPP_mean_EAfrClip.tif")
GPP = GPPI*1
proj4string(GPP) <- crs.geo

kernel100I = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_clips/kernel_res_cubicSP.tif")
kernel100 = kernel100I*1
proj4string(kernel100) <- crs.geo

#Create raster stack named env
env=stack(arid, access, prec, mean.temp, human.density, min.temp, Slope, Altitude, PET, DailyTempRange, max.temp, AnnualTempRange, prec.wet, prec.dry, GPP, kernel100) 

#Rename each raster in env so they can be referenced later
names(env) [1] <- "arid"
names(env) [2] <- "access"
names(env) [3] <- "prec"
names(env) [4] <- "mean.temp"
names(env) [5] <- "human.density"
names(env) [6] <- "min.temp"
names(env) [7] <- "Slope"
names(env) [8] <- "Altitude"
names(env) [9] <- "PET"
names(env) [10] <- "DailyTempRange"
names(env) [11] <- "max.temp"
names(env) [12] <- "AnnualTempRange"
names(env) [13] <- "prec.wet"
names(env) [14] <- "prec.dry"
names(env) [15] <- "GPP"
names(env) [16] <- "kernel100"

print("raster stack done")

value.raster = as.data.frame(getValues(arid))
names(value.raster) = "arid"
value.raster$access = getValues(access)
value.raster$prec = getValues(prec)
value.raster$mean.temp = getValues(mean.temp)
value.raster$human.density = getValues(human.density)
value.raster$min.temp = getValues(min.temp)
value.raster$Slope = getValues(Slope)
value.raster$Altitude = getValues(Altitude)
value.raster$PET = getValues(PET)
value.raster$DailyTempRange = getValues(DailyTempRange)
value.raster$max.temp = getValues(max.temp)
value.raster$AnnualTempRange = getValues(AnnualTempRange)
value.raster$prec.wet = getValues(prec.wet)
value.raster$prec.dry = getValues(prec.dry)
value.raster$GPP = getValues(GPP)
value.raster$kernel100 = getValues(kernel100)
save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/EAfrica_output/rasterstack_noLC_noFriction.RData")

EOF
