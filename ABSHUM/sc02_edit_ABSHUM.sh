#!/bin/bash
#SBATCH -p day
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_ABSHUM.sh.%J.out
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_ABSHUM.sh.%J.err
#SBATCH --job-name=sc02_edit_ABSHUM.sh

##   sbatch  ~/scripts/MOSQLAND/ABSHUM/sc02_edit_ABSHUM.sh

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/Florida_clips

#gdal_edit.py -a_ullr -180   90  180  -60 ABS50S.tif

#gdal_edit.py -a_ullr -180   90  180  -60 ABS50.tif

#gdalwarp -r cubic -tr 0.00833333333333333 0.00833333333333333  -te -180 -60 180 90  ABS50_BK.tif ABS50_res_cubic.tif -co BIGTIFF=YES -co COMPRESS=DEFLATE -co ZLEVEL=9 -overwrite

gdalwarp -r cubicspline -tr 0.00833333333333333 -0.00833333333333333  -te -180 -60 180 90  ABS50_BK.tif ABS50_res_cubicSP.tif -co BIGTIFF=YES -co COMPRESS=DEFLATE -co ZLEVEL=9 -overwrite


#gdal_edit.py -a_ullr -180   90  180  -60 ABS50_res.tif

#gdal_translate  -projwin -85 31.5 -79.8 24.0  ABS50S.tif  $OUTDIR/ABS50S_FloridaClip.tif

#gdal_translate  -projwin -85 31.5 -79.8 24.0  ABS50_res_cubic.tif  $OUTDIR/ABS50_cubic_FloridaClip.tif -co COMPRESS=DEFLATE -co ZLEVEL=9               

gdal_translate  -projwin -85 31.5 -79.8 24.0  ABS50_res_cubicSP.tif  $OUTDIR/ABS50_cubicSP_FloridaClip.tif -co COMPRESS=DEFLATE -co ZLEVEL=9


#gdal_translate -of AAIGrid $OUTDIR/ABS50S_FloridaClip.tif $OUTDIR/ABS50S_FloridaClip.asc
