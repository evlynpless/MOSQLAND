#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc27_edit_bio13.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc27_edit_bio13.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc27_edit_bio13.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio13

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio13/Florida_clips

for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  
gdal_edit.py -a_ullr -180   84  180  -90 bio13_${n}.tif 

done 


gdalbuildvrt -separate -overwrite  bio13.vrt    bio13_2000.tif   bio13_2001.tif   bio13_2002.tif   bio13_2003.tif   bio13_2004.tif   bio13_2005.tif   bio13_2006.tif   bio13_2007.tif   bio13_2008.tif   bio13_2009.tif    bio13_2010.tif    bio13_2011.tif    bio13_2012.tif    bio13_2013.tif

pkstatprofile  -f mean  -i  bio13.vrt  -o  bio13_mean.tif

gdal_edit.py -a_nodata -9999 bio13_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio13_mean.tif  $OUTDIR/bio13_mean_FloridaClip.tif

gdal_translate  -projwin -113.1 35.8 -100.5 28.7  bio13_mean.tif  SW_clip/bio13_mean_SWclip.tif

#gdal_translate -of AAIGrid bio13_mean_FloridaClip.tif bio13_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio13: precipitation of the wettest month
#It also averages over the years 2000-2013
