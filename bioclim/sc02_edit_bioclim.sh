#!/bin/bash                                                                                                 
#SBATCH -p scavenge                                                                                         
#SBATCH -N 1                                                                                                
#SBATCH -c 1                                                                                                
#SBATCH -t 4:00:00                                                                                          
#SBATCH --mail-type=ALL                                                                                     
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                    
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_bioclim.sh.%J.out                                
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_bioclim.sh.%J.err                               \
                                                                                                            
##   sbatch  ~/scripts/MOSQLAND/sc02_edit_bioclim.sh                                                           

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim/Florida_clips

for n in $(seq 1 19 ); do

gdal_edit.py -a_ullr -180   84  180  -90 CHELSA_bio10_${n}.tif

gdal_edit.py -a_nodata -999 CHELSA_bio10_${n}.tif

pksetmask -i CHELSA_bio10_${n}.tif -m CHELSA_bio10_${n}.tif -o CHELSA_bio10_${n}.tif --msknodata -32768 -nodata -999                                   \

gdal_translate  -projwin -85 31.5 -79.8 24.0  CHELSA_bio10_${n}.tif  $OUTDIR/CHELSA_bio10_${n}_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/CHELSA_bio10_${n}_FloridaClip.tif $OUTDIR/CHELSA_bio10_${n}_FloridaClip.asc                                                                                                      

done
