#!/bin/bash                                        
#SBATCH -p day                                     
#SBATCH --mem=100g                                
#SBATCH -t 24:00:00                                
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc17_rasterstack_allKernels.sh.%J.out                   
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc17_rasterstack_allKernels.sh.%J.err                  
#SBATCH --mail-type=ALL                            
#SBATCH --mail-user=evlyn.pless@yale.edu           
#SBATCH --job-name=sc17_rasterstack.sh                                                                                        

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc17_rasterstack_allKernels.sh                                                     
                 

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

aridI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2_Int16.tif")
arid = aridI*1
proj4string(arid) <- crs.geo

accessI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif")
access <- accessI*1
proj4string(access) <- crs.geo

precI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip2_Int16.tif")
prec <- precI*1
proj4string(prec) <- crs.geo

mean.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/NAm_clip/bio1_NAmClip2_Int16.tif")
mean.temp <- mean.tempI*1
proj4string(mean.temp) <- crs.geo

human.densityI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif")
human.density <- human.densityI*1
proj4string(human.density) <- crs.geo

frictionI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip2.tif")
friction <- frictionI*1
proj4string(friction) <- crs.geo

min.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip2_Int16.tif")
min.temp <- min.tempI*1
proj4string(min.temp) <- crs.geo

NeedleleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip2.tif")
Needleleaf = NeedleleafI*1
proj4string(Needleleaf) <- crs.geo

EvBroadleafI  = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip2.tif")
EvBroadleaf = EvBroadleafI*1
proj4string(EvBroadleaf) <- crs.geo

DecBroadleafI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip2.tif")
DecBroadleaf = DecBroadleafI*1
proj4string(DecBroadleaf) <- crs.geo

MiscTreesI =  raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip2.tif")
MiscTrees = MiscTreesI*1
proj4string(MiscTrees) <- crs.geo

ShrubsI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip2.tif")
Shrubs = ShrubsI*1
proj4string(Shrubs) <- crs.geo

HerbI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip2.tif")
Herb = HerbI*1
proj4string(Herb) <- crs.geo

CropI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip2.tif")
Crop <- CropI*1
proj4string(Crop) <- crs.geo

FloodI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip2.tif")
Flood = FloodI*1
proj4string(Flood) <- crs.geo

UrbanI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip2.tif")
Urban <- UrbanI*1
proj4string(Urban) <- crs.geo

SnowI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip2.tif")
Snow = SnowI*1
proj4string(Snow) <- crs.geo

BarrenI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip2.tif")
Barren = BarrenI*1
proj4string(Barren) <- crs.geo

WaterI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip2.tif")
Water = WaterI*1
proj4string(Water) <- crs.geo

SlopeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip2.tif")
Slope = SlopeI*1
proj4string(Slope) <- crs.geo

AltitudeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif")
Altitude = AltitudeI*1
proj4string(Altitude) <- crs.geo

PETI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET/NAm_clip/pet_mean_NAmClip2.tif")
PET = PETI*1
proj4string(PET) <- crs.geo

DailyTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio2/NAm_clip/bio2_mean_NAmClip2_Int16.tif")
DailyTempRange = DailyTempRangeI*1
proj4string(DailyTempRange) <- crs.geo

max.tempI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5/NAm_clip/bio5_mean_NAmClip2_Int16.tif")
max.temp = max.tempI*1
proj4string(max.temp) <- crs.geo

AnnualTempRangeI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio7/NAm_clip/bio7_mean_NAmClip2_Int16.tif")
AnnualTempRange = AnnualTempRangeI*1
proj4string(AnnualTempRange) <- crs.geo

prec.wetI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio13/NAm_clip/bio13_mean_NAmClip2_Int16.tif")
prec.wet = prec.wetI*1
proj4string(prec.wet) <- crs.geo

prec.dryI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio14/NAm_clip/bio14_mean_NAmClip2.tif")
prec.dry = prec.dryI*1
proj4string(prec.dry) <- crs.geo

GPPI = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly_mean/NAm_clip/GPP_mean_NAmClip2_Int16.tif")
GPP = GPPI*1
proj4string(GPP) <- crs.geo

kernel50I = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/kernel/KernelRas_50m_fnl.tif")
kernel50 = kernel50I*1
proj4string(kernel50) <- crs.geo

kernel100I = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/kernel/KernelRas_100m_fnl.tif")
kernel100 = kernel100I*1
proj4string(kernel100) <- crs.geo

kernel150I = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/kernel/KernelRas_150m_fnl.tif")
kernel150 = kernel150I*1
proj4string(kernel150) <- crs.geo

kernel200I = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/kernel/KernelRas_200m_fnl.tif")
kernel200 = kernel200I*1
proj4string(kernel200) <- crs.geo

#Create raster stack named env
env=stack(arid, access, prec, mean.temp, human.density, friction, min.temp, Needleleaf, EvBroadleaf, DecBroadleaf, MiscTrees, Shrubs, Herb, Crop, Flood, Urban, Snow, Barren, Water, Slope, Altitude, PET, DailyTempRange, max.temp, AnnualTempRange, prec.wet, prec.dry, GPP, kernel50, kernel100, kernel150, kernel200)

#Rename each raster in env so they can be referenced later
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
names(env) [29] <- "kernel50"
names(env) [30] <- "kernel100" 
names(env) [31] <- "kernel150" 
names(env) [32] <- "kernel200" 

print("raster stack done")

value.raster = as.data.frame(getValues(arid))
names(value.raster) = "arid"
value.raster$access = getValues(access)
value.raster$prec = getValues(prec)
value.raster$mean.temp = getValues(mean.temp)
value.raster$human.density = getValues(human.density)
value.raster$friction = getValues(friction)
value.raster$min.temp = getValues(min.temp)
value.raster$Needleleaf = getValues(Needleleaf)
value.raster$EvBroadleaf = getValues(EvBroadleaf)
value.raster$DecBroadleaf = getValues(DecBroadleaf)
value.raster$MiscTrees = getValues(MiscTrees)
value.raster$Shrubs = getValues(Shrubs)
value.raster$Herb = getValues(Herb)
value.raster$Crop = getValues(Crop)
value.raster$Flood = getValues(Flood)
value.raster$Urban = getValues(Urban)
value.raster$Snow = getValues(Snow)
value.raster$Barren = getValues(Barren)
value.raster$Water = getValues(Water)
value.raster$Slope = getValues(Slope)
value.raster$Altitude = getValues(Altitude)
value.raster$PET = getValues(PET)
value.raster$DailyTempRange = getValues(DailyTempRange)
value.raster$max.temp = getValues(max.temp)
value.raster$AnnualTempRange = getValues(AnnualTempRange)
value.raster$prec.wet = getValues(prec.wet)
value.raster$prec.dry = getValues(prec.dry)
value.raster$GPP = getValues(GPP)
value.raster$kernel50 = getValues(kernel50)
value.raster$kernel100 = getValues(kernel100) 
value.raster$kernel150 = getValues(kernel150) 
value.raster$kernel200 = getValues(kernel200) 
save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc17_rasterstack_allKernels.RData")

EOF
