#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc18_edit_bio1.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc18_edit_bio1.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc18_edit_bio1.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1                          \

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio1/Florida_clips

#for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  

#gdal_edit.py -a_ullr -180   84  180  -90 bio1_${n}.tif 

#done 


gdalbuildvrt -separate -overwrite  bio1.vrt    bio1_2000.tif   bio1_2001.tif   bio1_2002.tif   bio1_2003.tif   bio1_2004.tif   bio1_2005.tif   bio1_2006.tif   bio1_2007.tif   bio1_2008.tif   bio1_2009.tif    bio1_2010.tif    bio1_2011.tif    bio1_2012.tif    bio1_2013.tif

pkstatprofile  -f mean  -i  bio1.vrt  -o  bio1_mean.tif

gdal_edit.py -a_nodata -9999 bio1_mean_FloridaClip.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio1_mean.tif  $OUTDIR/bio1_mean_FloridaClip.tif

gdal_translate  -projwin -113.1 35.8 -100.5 28.7 bio1_mean.tif SW_clip/bio1_mean_SWclip.tif

#This script is for editing the chelsa version of bio1: average annual temperature
#It also averages over the years 2000-2013
