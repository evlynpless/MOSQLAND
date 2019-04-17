#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc26_edit_bio7.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc26_edit_bio7.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc26_edit_bio7.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio7                          \

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5/Florida_clips

for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  
gdal_edit.py -a_ullr -180   84  180  -90 bio7_${n}.tif 
done 

    
gdalbuildvrt -separate -overwrite  bio7.vrt    bio7_2000.tif   bio7_2001.tif   bio7_2002.tif   bio7_2003.tif   bio7_2004.tif   bio7_2005.tif   bio7_2006.tif   bio7_2007.tif   bio7_2008.tif   bio7_2009.tif    bio7_2010.tif    bio7_2011.tif    bio7_2012.tif    bio7_2013.tif

pkstatprofile  -f mean  -i  bio7.vrt  -o  bio7_mean.tif

gdal_edit.py -a_nodata -9999 bio7_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio7_mean.tif  $OUTDIR/bio7_mean_FloridaClip.tif

#gdal_translate -of AAIGrid bio5_mean_FloridaClip.tif bio5_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio7: Temperature Annual Range
#It also averages over the years 2000-2013
