#!/bin/bash  
#SBATCH -p scavenge
#SBATCH -N 1                                                                                                                     
#SBATCH -c 1                                                                                                                     
#SBATCH -t 4:00:00                                                                                                               
#SBATCH --mail-type=ALL                                                                                           
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                                      
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc20_edit_bio12.sh.%J.out                                                     
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc20_edit_bio12.sh.%J.err                                                                                         
##   sbatch  ~/scripts/MOSQLAND/sc20_edit_bio12.sh
 
cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12                          

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/Florida_clips

#for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013; do                                                                                                                   
#gdal_edit.py -a_ullr -180   84  180  -90 bio12_${n}.tif                                                                                                                                                                                                           
#done                                                                                                                                                                                          

gdalbuildvrt -separate -overwrite  bio12.vrt    bio12_2000.tif   bio12_2001.tif   bio12_2002.tif   bio12_2003.tif   bio12_2004.tif   bio12_2005.tif   bio12_2006.tif   bio12_2007.tif   bio12_2008.tif   bio12_2009.tif    bio12_2010.tif    bio12_2011.tif     bio12_2012.tif     bio12_2013.tif
 
pkstatprofile  -f mean  -i  bio12.vrt  -o  bio12_mean.tif

gdal_edit.py -a_nodata -9999 bio12_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio12_mean.tif  $OUTDIR/bio12_mean_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/bio12_mean_FloridaClip.tif $OUTDIR/bio12_mean_FloridaClip.asc

#This script is for editing the chelsa version of bio12: annual precipitation                                                                                                          
#It also averages over the years 2000-2013 
