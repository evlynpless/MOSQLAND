#!/bin/bash                                                                         
                                                                                  
#SBATCH -p day                                                                       
                                                                                  
#SBATCH -n 1 -c 20  -N 1                                                            
                                                                               
#SBATCH -t 10:00:00                                                                  
                                                                              
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc03_RF_FST.R.sh.%J.out          
                                                                                  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stdout/sc02_RF_FST.R.sh.%J.err        
                                                                                     
#SBATCH --mail-type=ALL                                                              
#SBATCH --mail-user=email                                                           
#SBATCH --job-name=sc03_RF_FST.sh                                                  

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc03_RF_FST.R.sh  

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla --no-readline   -q  <<'EOF'

library("sp")
library("spatstat")
library("maptools")
library("raster")

OUTDIR="/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/" 

###############################################
#Plot points as SpatialPoints
###############################################

G.table <- read.csv(paste0(OUTDIR,"FST_list_Florida.csv"), header=TRUE) 

G.coordinates1 <- G.table[,c(5,3)] 
G.points1 <- SpatialPoints(G.coordinates1)  # ... converted into a spatial object

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system
proj4string(G.points1) <- crs.geo  # define projection system of our data

G.coordinates2 <- G.table[,c(6,4)]
G.points2 <- SpatialPoints(G.coordinates2)
proj4string(G.points2) <- crs.geo


###############################################
#Plot lines as SpatialLines:
###############################################
#create dataframes of begin and end coordinates from a file:

begin.table <- G.table[,c(5,3)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- G.table[,c(6,4)]
end.coord <- end.table
coordinates(end.coord) <- c("long2", "lat2")

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

plot(p,add=T)

spatial.p <- as(p, "SpatialLines")
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system
proj4string(spatial.p) <- crs.geo  # define projection system of our data


abshum = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/NAm_clip/ABS50_res_NAmClip.tif")                        
proj4string(abshum) <- crs.geo      
access = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip.tif")
proj4string(access) <- crs.geo    
arid = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/NAm_clip/AI_annual_NAmClip.tif ")                         
proj4string(arid) <- crs.geo     
mean.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/NAm_clip/bio1_NAmClip.tif")                      
proj4string(mean.temp) <- crs.geo     
min.temp = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip.tif")                  
proj4string(min.temp) <- crs.geo     
prec = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip.tif")                    
proj4string(prec) <- crs.geo     
GSHL = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip.tif")              
proj4string(GSHL) <- crs.geo                                                                                                                 
Needleleaf = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip.tif")           
proj4string(Needleleaf) <- crs.geo 
Evergreen  = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip.tif")           
proj4string(Evergreen) <- crs.geo 
Deciduous = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip.tif")            
proj4string(Deciduous) <- crs.geo 
OtherTrees = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip.tif")           
proj4string(OtherTrees) <- crs.geo 
Shrubs = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip.tif")               
proj4string(Shrubs) <- crs.geo 
Herb = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip.tif")                 
proj4string(Herb) <- crs.geo 
Crop = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip.tif")                 
proj4string(Crop) <- crs.geo 
Flood = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip.tif")                
proj4string(Flood) <- crs.geo 
Urban = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip.tif")                
proj4string(Urban) <- crs.geo 
Snow = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip.tif")                
proj4string(Snow) <- crs.geo 
Barren = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip.tif")              
proj4string(Barren) <- crs.geo 
Water = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip.tif")               
proj4string(Water) <- crs.geo 
altitude = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip.tif")    
proj4string(altitude) <- crs.geo 
slope = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip.tif")  
proj4string(slope) <- crs.geo 
friction = raster("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip.tif")
proj4string(friction) <- crs.geo  

env=stack(abshum, arid, mean.temp, min.temp, prec, GSHL, Needleleaf, Evergreen, Deciduous, OtherTrees, Shrubs, Herb, Crop, Flood, Urban, Snow, Barren, Water, altitude, slope, friction ) 


####Extract final line sums

line_sum = raster::extract(env, spatial.p, fun=sum, na.rm=TRUE)

write.csv(line_sum, "line_means.csv")

EOF
