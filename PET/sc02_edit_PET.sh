#!/bin/bash                                                                                              
#SBATCH -p scavenge                                                                                      
#SBATCH -N 1                                                                                             
#SBATCH -c 1                                                                                             
#SBATCH -t 4:00:00                                                                                       
#SBATCH --mail-type=ALL                                                                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                                                 
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc02_edit_PET.sh%J.out                             
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc02_edit_PET.sh.%J.err                            \
                                                                                                        
##   sbatch  ~/scripts/MOSQLAND/sc02_edit_PET.sh                                                        

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET/Florida_clips

gdalbuildvrt -separate -overwrite  pet.vrt    pet_he_1.tif   pet_he_2.tif   pet_he_3.tif   pet_he_4.tif   pet_he_5.tif   pet_he_6.tif   pet_he_7.tif   pet_he_8.tif   pet_he_9.tif    pet_he_10.tif    pet_he_11.tif   pet_he_12.tif

pkstatprofile  -f mean  -i  pet.vrt  -o  pet_mean.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  pet_mean.tif  $OUTDIR/pet_mean_FloridaClip.tif

pksetmask -i $OUTDIR/pet_mean_FloridaClip.tif -m $OUTDIR/pet_mean_FloridaClip.tif -o $OUTDIR/pet_mean_FloridaClip.tif --msknodata -16257 -nodata -999

gdal_edit.py -a_nodata -999 $OUTDIR/pet_mean_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/pet_mean_FloridaClip.tif $OUTDIR/pet_mean_FloridaClip.asc

gdal_translate  -projwin -113.1 35.8 -100.5 28.7 pet_mean.tif   SW_clip/pet_mean_SWclip.tif
