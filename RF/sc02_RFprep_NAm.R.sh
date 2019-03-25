#!/bin/bash                                                                                                                    
#SBATCH -p day                                                                                                                 
#SBATCH -n 1 -c 20  -N 1                                                                                                       
#SBATCH -t 10:00:00                                                                                                            
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc02_RFprep_FL.sh.%J.out                                                    
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stdout/sc02_RFprep_FL.sh.%J.err                                                    
#SBATCH --mail-type=ALL                                                                                                        
#SBATCH --mail-user=email                                                                                                      
#SBATCH --job-name=sc02_RFprep_FL.sh                                                                                           

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc02_RFprep_FL.R.sh                                                      
                    

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla --no-readline   -q  <<'EOF'  


library(raster)
#library (dismo)
#library (XML)
#library (rgeos)
#library (maptools)
#library(sp)
#library(rasterVis)
# library(rgdal)
#library(rJava)

OUTDIR="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/"

dat_orig <- read.csv(paste0(OUTDIR,"NAm_coors.csv"), header=TRUE)

head(dat_orig)

dat_orig$X

dat <- data.frame(dat_orig, lat=dat_orig$Y, lon=dat_orig$X)

coordinates(dat) <- c('X', 'Y')

projection(dat)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#Florida data
abshum = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/NAm_clip/ABS50_res_NAmClip.tif")
#access = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip.tif")
arid = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/NAm_clip/AI_annual_NAmClip.tif ")
mean.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/NAm_clip/bio1_NAmClip.tif")                   
min.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip.tif")
prec = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip.tif")
GSHL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip.tif")
Needleleaf = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip.tif")
Evergreen  = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip.tif")
Deciduous = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip.tif")
OtherTrees = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip.tif")
Shrubs = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip.tif")
Herb = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip.tif")
Crop = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip.tif")
Flood = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip.tif")
Urban = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip.tif")
Snow = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip.tif")
Barren = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip.tif")
Water = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip.tif")
altitude = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip.tif")
slope = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip.tif")


env=stack(abshum, arid, mean.temp, min.temp, prec, GSHL, Needleleaf, Evergreen, Deciduous, OtherTrees, Shrubs, Herb, Crop, Flood, Urban, Snow, Barren, Water, altitude, slope )

points = extract(env, dat, sp=T)

points_buffer = extract(env, dat, buffer=3, fun=mean, sp=T)

head(points)

write.table(points, paste0(OUTDIR,"RF_out_NAm.txt"))

write.table(points_buffer, paste0(OUTDIR, "RF_out_buffer_NAm.txt"))

EOF
