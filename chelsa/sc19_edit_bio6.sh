#!/bin/bash                        
#SBATCH -p scavenge                
#SBATCH -N 1                       
#SBATCH -c 1                     
#SBATCH -t 4:00:00                 
#SBATCH --mail-type=ALL            
#SBATCH --mail-user=evlyn.pless@yale.edu                                                              
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc19_edit_bio6.sh.%J.out                          
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc19_edit_bio6.sh.%J.err                          
                                                                    
##   sbatch  ~/scripts/MOSQLAND/sc19_edit_bio6.sh                                                     
                                                                                         
cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6                          

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio6/Florida_clips

#for n in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 ; do                                                                                                                     
#gdal_edit.py -a_ullr -180   84  180  -90 bio6_${n}.tif                                                                                                                                                                                                                                           
#done                                                                                                                

gdalbuildvrt -separate -overwrite  bio6.vrt    bio6_2000.tif   bio6_2001.tif   bio6_2002.tif   bio6_2003.tif   bio6_2004.tif   bio6_2005.tif   bio6_2006.tif   bio6_2007.tif   bio6_2008.tif   bio6_2009.tif    bio6_2010.tif    bio6_2011.tif    bio6_2012.tif     bio6_2013.tif

pkstatprofile  -f mean  -i  bio6.vrt  -o  bio6_mean.tif

gdal_edit.py -a_nodata -9999 bio6_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  bio6_mean.tif  $OUTDIR/bio6_mean_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/bio6_mean_FloridaClip.tif $OUTDIR/bio6_mean_FloridaClip.asc

#making a version with all positive values for circuitscape
gdal_calc.py -A $OUTDIR/bio6_mean_FloridaClip.tif  --outfile=$OUTDIR/bio6_mean_FloridaClip_positive.tif --NoDataValue=-9999  --calc="A+691"

gdal_translate -of AAIGrid $OUTDIR/bio6_mean_FloridaClip_positive.tif $OUTDIR/bio6_mean_FloridaClip_positive.asc


#This script is for editing the chelsa version of bio6: min temperature of coldest month                                                                                   #It also averages over the years 2000-2011  
