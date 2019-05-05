1;95;0c#!/bin/bash                                                                                              
#SBATCH -p scavenge                                                                                      
#SBATCH -N 1                                                                                             
#SBATCH -c 1                                                                                             
#SBATCH -t 4:00:00                                                                                       
#SBATCH --mail-type=ALL                                                                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                 
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_altitude.sh%J.out                             
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_altitude.sh.%J.err                            \
                                                                                                         
##   sbatch  ~/scripts/MOSQLAND/sc02_edit_altitude.sh 

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude/Florida_clips

gdal_edit.py -a_ullr -180   84  180  -90 altitude_1KMmedian_MERIT.tif

gdal_edit.py -a_nodata -999 altitude_1KMmedian_MERIT.tif

pksetmask -i altitude_1KMmedian_MERIT.tif -m altitude_1KMmedian_MERIT.tif -o altitude_1KMmedian_MERIT.tif --msknodata -9999 -nodata -999   
 
gdal_translate  -projwin -85 31.5 -79.8 24.0  altitude_1KMmedian_MERIT.tif  $OUTDIR/altitude_1KMmedian_MERIT_FloridaClip.tif

gdal_translate  -projwin -113.1 35.8 -100.5 28.7 altitude_1KMmedian_MERIT.tif   SW_clip/altitude_1KMmedian_MERIT_SWclip.tif

#gdal_translate -of AAIGrid $OUTDIR/altitude_1KMmedian_MERIT_FloridaClip.tif $OUTDIR/altitude_1KMmedian_MERIT_FloridaClip.asc
