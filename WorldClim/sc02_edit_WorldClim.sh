#!/bin/bash
                                                                                                                               
#SBATCH -p scavenge
                                                                                                                               
#SBATCH -N 1
                                                                                                                               
#SBATCH -c 1
                                                                                                                               
#SBATCH -t 4:00:00
                                                                                                                               
#SBATCH --mail-type=ALL
                                                                                                                               
#SBATCH --mail-user=evlyn.pless@yale.edu
                                                                                                                               
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_WorldClim.sh.%J.out                                                   
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_WorldClim.sh.%J.err                                                   

##   sbatch  ~/scripts/MOSQLAND/sc02_edit_WorldClim.sh                                                                              

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/WorldClim

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/WorldClim/Florida_clips

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 ; do

#gdal_edit.py -a_ullr -180   84  180  -90 wc2.0_bio_30s_${n}.tif                                                    

gdal_edit.py -a_nodata -999 wc2.0_bio_30s_${n}.tif                                                                    

pksetmask -i wc2.0_bio_30s_${n}.tif -m wc2.0_bio_30s_${n}.tif -o wc2.0_bio_30s_${n}.tif --msknodata -1.69999999999999994e+308 -nodata -999                                                                                                       
gdal_translate  -projwin -85 31.5 -79.8 24.0  wc2.0_bio_30s_${n}.tif  $OUTDIR/wc2.0_bio_30s_${n}.tif

gdal_translate -of AAIGrid  $OUTDIR/wc2.0_bio_30s_${n}.tif $OUTDIR/wc2.0_bio_30s_${n}.asc                                                                                                                         
done
