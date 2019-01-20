#!/bin/bash                                                                                                                                                     
#SBATCH -p day                                                                                                                                                  
#SBATCH -n 1 -c 20  -N 1                                                                                                                                        
#SBATCH -t 10:00:00                                                                                                                                             
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_RFprep_FL.sh.%J.out                                                                         
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stdout/sc01_RFprep_FL.sh.%J.err                                                                         
#SBATCH --mail-type=ALL                                                                                                                                         
#SBATCH --mail-user=email                                                                                                                                       
#SBATCH --job-name=sc01_RFprep_FL.sh                                                                                                                     

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc01_RFprep_FL.sh                                                                           


module load Apps/R/3.3.2-generic


R --vanilla --no-readline   -q  <<'EOF'  

library(raster)
#library (dismo)
#library (XML)
#library (rgeos)
#library (maptools)
#library(sp)
#library(rasterVis)
#library(rgdal)
#library(rJava)

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF

dat_orig <- read.table($OUTDIR/FL_coors.txt, header=TRUE)

head(dat_orig)

dat_orig$X

data.frame(src="gbif",obs=1,lat=gbif_points$lat,lon=gbif_points$lon)

#dat <- data.frame(src="dat_orig", obs=1, lat=dat_orig$Y, lon=dat_orig$X)
dat <- data.frame(dat_orig, lat=dat_orig$Y, lon=dat_orig$X)

coordinates(dat) <- c('X', 'Y')
projection(dat)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#World data
#access = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/accessibility_to_cities_2015_v1.0.tif")  ; NAvalue(access) <- -999
#arid = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/AI_annual.tif")  ; NAvalue(arid) <- -999
#GSHL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif")  ; NAvalue(GSHL) <- 0
#mean.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/CHELSA_bio10_1.tif")  ; NAvalue(mean.temp) <- -999
#min.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/CHELSA_bio10_6.tif")  ; NAvalue(min.temp) <- -999
#mean.prec = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/CHELSA_bio10_12.tif")  ; NAvalue(mean.prec) <- -999
#urban = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/consensus_full_class_9.tif")  ; NAvalue(urban) <- 0

#Florida data
access_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips/accessibility_to_cities_2015_v1.0_FloridaClip.tif")  ; NAvalue(access_FL) <- -999
arid_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/Florida_clips/AI_annual_FloridaClip.tif ")  ; NAvalue(arid_FL) <- -9999
GSHL_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/Florida_clips/AI_annual_FloridaClip.tif ")  ; NAvalue(GSHL_FL) <- 0
mean.temp_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/Florida_clips/CHELSA_bio10_1_FloridaClip.tif")  ; NAvalue(mean.temp_FL) <- -999
min.temp_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/Florida_clips/CHELSA_bio10_6_FloridaClip.tif")  ; NAvalue(min.temp_FL) <- -999
urban_FL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/Florida_clips/consensus_full_class_9_FloridaClip.tif")  ; NAvalue(urban_FL) <- 0


# set the extend 
#do i need to do this if they're all the same extent already?
#ext.glob <- extent (-180,180,-90,84)
#setExtent (prec.glob , ext.glob)

# crop and reset the extend
#prob not necessary
#ext.crop <- extent (-60,-30,-30,0)
#prec = crop(prec.glob , ext.crop) ; setExtent (prec , ext.crop)

env=stack(access_FL, arid_FL, GSHL_FL, mean.temp_FL, min.temp_FL, urban_FL)

points = extract(env, dat, sp=T)

write.table(points, file=$OUTDIR/RF_out.txt)

EOF
