#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc24_reclass_bio12.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc24_reclass_bio12.sh.%J.err                   
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc24_reclass_bio12.sh 

export DIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/bio12/Florida_clips



echo 1   10   20   30   40   50   60   70   80   90  100   linpos    >  $DIR/weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     linneg    >> $DIR/weight.txt
echo 100 80   60   40   20   1    20   40   60   80  100   linpar    >> $DIR/weight.txt
echo 1   25   50   75   100  150  200  250  325  400 500   exppos    >> $DIR/weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     expneg    >> $DIR/weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     exppar    >> $DIR/weight.txt



export min=$(pkstat -nodata -9999 -mm -i  $DIR/bio12_mean_FloridaClip.tif   | awk '{ print $2    }' )
export binrange=$(pkstat -nodata -9999 -mm -i  $DIR/bio12_mean_FloridaClip.tif   | awk '{ print ($4-$2) / 11   }' )


cat $DIR/weight.txt | xargs -n 12 -P 2  bash -c $'
echo  lastcol ${12} 

pkgetmask    -min $min                                  -max $(echo   $min + $binrange       | bc ) -ot Int16  -data $1  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_a_${12}.tif
pkgetmask    -min $(echo  $min + $binrange       | bc ) -max $(echo   $min + $binrange "*" 2 | bc ) -ot Int16  -data $2  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_b_${12}.tif
pkgetmask    -min $(echo  $min + $binrange "*" 2 | bc ) -max $(echo   $min + $binrange "*" 3 | bc ) -ot Int16  -data $3  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_c_${12}.tif
pkgetmask    -min $(echo  $min + $binrange "*" 3 | bc ) -max $(echo   $min + $binrange "*" 4 | bc ) -ot Int16  -data $4  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_d_${12}.tif
pkgetmask    -min $(echo  $min + $binrange "*" 4 | bc ) -max $(echo   $min + $binrange "*" 5 | bc ) -ot Int16  -data $5  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_e_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 5 | bc ) -max $(echo   $min + $binrange "*" 6 | bc ) -ot Int16  -data $6  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_f_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 6 | bc ) -max $(echo   $min + $binrange "*" 7 | bc ) -ot Int16  -data $7  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_g_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 7 | bc ) -max $(echo   $min + $binrange "*" 8 | bc ) -ot Int16  -data $8  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_h_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 8 | bc ) -max $(echo   $min + $binrange "*" 9 | bc ) -ot Int16  -data $9  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_i_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 9 | bc ) -max $(echo   $min + $binrange "*" 10 | bc ) -ot Int16  -data ${10}  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_j_${12}.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 10 | bc ) -max $(echo   $min + $binrange "*" 11 | bc ) -ot Int16  -data ${11}  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i $DIR/bio12_mean_FloridaClip.tif -o $DIR/input_k_${12}.tif 

gdalbuildvrt -separate -overwrite  $DIR/input${12}.vrt    $DIR/input_?_${12}.tif 

pkstatprofile  -f max -ot Int32  -i  $DIR/input${12}.vrt  -o  $DIR/combined_${12}_tmp.tif

pksetmask -ot Int32 -m  $DIR/bio12_mean_FloridaClip.tif -msknodata -9999 -nodata -9999    -co COMPRESS=DEFLATE -co ZLEVEL=9   -i  $DIR/combined_${12}_tmp.tif -o  $DIR/combined_${12}.tif

gdal_translate -of AAIGrid $DIR/combined_${12}.tif $DIR/combined_${12}.asc

rm -f $DIR/combined_${12}_tmp.tif

rm -f  $DIR/input_?_${12}.tif 

rm -f $DIR/input${12}.vrt
' _ 



