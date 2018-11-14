#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc04_test_bins_A.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc04_test_bins_A.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc04_test_bins_A.sh 

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips


echo 100 10 8  ex      > weight.txt
echo 4   2  1  pa    >> weight.txt



export min=$(pkstat -nodata -999 -mm -i   accessibility_to_cities_2015_v1.0_FloridaClip.tif    | awk '{ print $(2)    }' )
export binrange=$(pkstat -nodata -999 -mm -i  accessibility_to_cities_2015_v1.0_FloridaClip.tif    | awk '{ print ($4-$2) / 3   }' )


cat weight.txt | xargs -n 4 -P 2 bash -c $' 

pkgetmask    -min $min                                  -max $(echo   $min + $binrange       | bc )  -data $(1)  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$(1)_a.tif 
pkgetmask    -min $(echo  $min + $binrange       | bc ) -max $(echo   $min + $binrange "*" 2 | bc )  -data $(2)  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$(2)_b.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 2 | bc ) -max $(echo   $min + $binrange "*" 3 | bc )  -data $(3)  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$(3)_c.tif 

gdalbuildvrt -separate  input$(4).vrt    input$(1)_a.tif    input$(2)_b.tif   input$(3)_c.tif 

pkstatprofile  -f max  -i  input$(4).vrt  -o  input$(4)_max.tif

' _ 


#the goal here is to test if $# and $(#) behave the same
