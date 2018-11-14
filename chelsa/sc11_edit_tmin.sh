#!/bin/bash
                                                                                                        
#SBATCH -p scavenge 
                                                                                                
#SBATCH -N 1    
                                                                                                    
#SBATCH -c 1  
                                                                                                      
#SBATCH -t 4:00:00    
                                                                                              
#SBATCH --mail-type=ALL 
                                                                                            
#SBATCH --mail-user=evlyn.pless@yale.edu 
                                                                           
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc11_edit_tmin.sh.%J.out 
                                       
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc11_edit_tmin.sh.%J.err                                       
                                                                                                                    
##   sbatch  ~/scripts/MOSQLAND/sc11_edit_tmin.sh                                                                   

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmin
                                                                  
OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmin/Florida_clips

for n in $(seq 1 12 ) ; do

gdal_edit.py -a_ullr -180   84  180  -90 CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif                             

gdal_edit.py -a_nodata -999 CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif

pksetmask -i CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif -m CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif -o CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif --msknodata -32768 -nodata -999                                                           

gdal_translate  -projwin -85 31.5 -79.8 24.0  CHELSA_tmin10_${n}_1979-2013_V1.2_land.tif  $OUTDIR/CHELSA_tmin10_${n}_FloridaClip.tif

#gdal_translate -of AAIGrid $OUTDIR/CHELSA_tmin10_${n}_FloridaClip.tif $OUTDIR/CHELSA_tmin10_${n}_FloridaClip.asc

done
