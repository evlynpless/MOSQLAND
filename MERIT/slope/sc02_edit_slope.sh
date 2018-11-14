#!/bin/bash                                                                                              
#SBATCH -p scavenge                                                                                      
#SBATCH -N 1                                                                                             
#SBATCH -c 1                                                                                             
#SBATCH -t 4:00:00                                                                                       
#SBATCH --mail-type=ALL                                                                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                 
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_slope.sh%J.out                             
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_slope.sh.%J.err                           
                                                                                                         
##   sbatch  ~/scripts/MOSQLAND/sc02_edit_slope.sh 

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/slope/Florida_clips

#gdal_edit.py -a_ullr -180   84  180  -90 slope_1KMmedian_MERIT.tif

gdal_edit.py -a_nodata -999 slope_1KMmedian_MERIT.tif

pksetmask -i slope_1KMmedian_MERIT.tif -m slope_1KMmedian_MERIT.tif -o slope_1KMmedian_MERIT.tif --msknodata -9999 -nodata -999   
 
gdal_translate  -projwin -85 31.5 -79.8 24.0  slope_1KMmedian_MERIT.tif  $OUTDIR/slope_1KMmedian_MERIT_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/slope_1KMmedian_MERIT_FloridaClip.tif $OUTDIR/slope_1KMmedian_MERIT_FloridaClip.asc

#Making a raster with all positive values for circuitscape 
gdal_calc.py -A $OUTDIR/slope_1KMmedian_MERIT_FloridaClip.tif  --outfile=$OUTDIR/slope_1KMmedian_MERIT_FloridaClip_positive.tif --NoDataValue=-999  --calc="A+0.01"

gdal_translate -of AAIGrid $OUTDIR/slope_1KMmedian_MERIT_FloridaClip_positive.tif $OUTDIR/slope_1KMmedian_MERIT_FloridaClip_positive.asc

