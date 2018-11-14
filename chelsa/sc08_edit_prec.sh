#!/bin/bash                                                                                                                      
#SBATCH -p scavenge                                                                                                              
#SBATCH -N 1                                                                                                                     
#SBATCH -c 1                                                                                                                     
#SBATCH -t 4:00:00                                                                                                               
#SBATCH --mail-type=ALL                                                                                                          
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                                         
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc08_edit_prec.sh.%J.out                                                
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc08_edit_prec.sh.%J.err                                                

##   sbatch  ~/scripts/MOSQLAND/sc08_edit_prec.sh

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/prec

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/prec/Florida_clips

for n in $(seq 1 12 ) ; do

#gdal_edit.py -a_ullr -180   84  180  -90 CHELSA_prec_${n}_V1.2_land.tif                                      

#gdal_edit.py -a_nodata -999 CHELSA_prec_${n}_V1.2_land.tif                                                   

#pksetmask -i CHELSA_prec_${n}_V1.2_land.tif -m CHELSA_prec_${n}_V1.2_land.tif -o CHELSA_prec_${n}_V1.2_land.tif --msknodata -32767 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0  CHELSA_prec_${n}_V1.2_land.tif  $OUTDIR/CHELSA_prec_${n}_V1.2_land_FloridaClip.tif

#gdal_translate -of AAIGrid  $OUTDIR/CHELSA_prec_${n}_V1.2_land_FloridaClip.tif $OUTDIR/CHELSA_prec_${n}_V1.2_land_FloridaClip.asc

done
