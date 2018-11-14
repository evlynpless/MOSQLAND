#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc22_edit_bio5.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc22_edit_bio5.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc22_edit_bio5.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5                          \

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio5/Florida_clips

#for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  

#gdal_edit.py -a_ullr -180   84  180  -90 bio1_${n}.tif 

#done 


gdalbuildvrt -separate -overwrite  bio5.vrt    bio5_2000.tif   bio5_2001.tif   bio5_2002.tif   bio5_2003.tif   bio5_2004.tif   bio5_2005.tif   bio5_2006.tif   bio5_2007.tif   bio5_2008.tif   bio5_2009.tif    bio5_2010.tif    bio5_2011.tif    bio5_2012.tif    bio5_2013.tif

pkstatprofile  -f mean  -i  bio5.vrt  -o  bio5_mean.tif

gdal_edit.py -a_nodata -9999 bio5_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio5_mean.tif  $OUTDIR/bio5_mean_FloridaClip.tif

gdal_translate -of AAIGrid bio5_mean_FloridaClip.tif bio5_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio5: max temperature of hottest month
#It also averages over the years 2000-2013
