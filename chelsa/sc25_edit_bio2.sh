#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc25_edit_bio2.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc25_edit_bio2.sh.%J.err
                                                                                                              
##   sbatch  ~/scripts/MOSQLAND/sc25_edit_bio2.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio2                          

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio2/Florida_clips

#for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  

#gdal_edit.py -a_ullr -180   84  180  -90 bio1_${n}.tif 

#done 


gdalbuildvrt -separate -overwrite  bio2.vrt    bio2_2000.tif   bio2_2001.tif   bio2_2002.tif   bio2_2003.tif   bio2_2004.tif   bio2_2005.tif   bio2_2006.tif   bio2_2007.tif   bio2_2008.tif   bio2_2009.tif    bio2_2010.tif    bio2_2011.tif    bio2_2012.tif    bio2_2013.tif

pkstatprofile  -f mean  -i  bio2.vrt  -o  bio2_mean.tif

gdal_edit.py -a_nodata -9999 bio2_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio2_mean.tif  $OUTDIR/bio2_mean_FloridaClip.tif

gdal_translate  -projwin -113.1 35.8 -100.5 28.7  bio2_mean.tif   SW_clip/ bio2_mean_SWclip.tif

#gdal_translate -of AAIGrid bio2_mean_FloridaClip.tif bio2_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio2: Mean diurnal (temp) range
#It also averages over the years 2000-2013
