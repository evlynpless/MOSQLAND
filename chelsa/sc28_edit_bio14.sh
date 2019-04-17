#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc28_edit_bio14.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc28_edit_bio14.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc28_edit_bio14.sh                                               

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio14

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio14/Florida_clips

for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 ; do                                                                                  

gdal_edit.py -a_ullr -180   84  180  -90 bio14_${n}.tif 

done 


gdalbuildvrt -separate -overwrite  bio14.vrt    bio14_2000.tif   bio14_2001.tif   bio14_2002.tif   bio14_2003.tif   bio14_2004.tif   bio14_2005.tif   bio14_2006.tif   bio14_2007.tif   bio14_2008.tif   bio14_2009.tif    bio14_2010.tif    bio14_2011.tif    bio14_2012.tif    bio14_2013.tif

pkstatprofile  -f mean  -i  bio14.vrt  -o  bio14_mean.tif

gdal_edit.py -a_nodata -9999 bio14_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio14_mean.tif  $OUTDIR/bio14_mean_FloridaClip.tif

#gdal_translate -of AAIGrid bio5_mean_FloridaClip.tif bio5_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio14: precipitation of the driest month
#It also averages over the years 2000-2013
