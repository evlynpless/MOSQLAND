#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc04_test_bins_B.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc04_test_bins_B.sh.%J.err                   \
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc04_test_bins_B.sh 

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips


echo 100 10 8  4  ex      > weight.txt
echo 2   4  2  1  pa    >> weight.txt
echo 2   3  4  5  lin >> weight.txt



export min=$(pkstat -nodata -999 -mm -i   accessibility_to_cities_2015_v1.0_FloridaClip.tif    | awk '{ print $2    }' )
export binrange=$(pkstat -nodata -999 -mm -i  accessibility_to_cities_2015_v1.0_FloridaClip.tif    | awk '{ print ($4-$2) / 3   }' )


cat weight.txt | xargs -n 5 -P 2 bash -c $' 

pkgetmask    -min $min                                  -max $(echo   $min + $binrange       | bc )  -data $1  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$1_a.tif 
pkgetmask    -min $(echo  $min + $binrange       | bc ) -max $(echo   $min + $binrange "*" 2 | bc )  -data $2  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$2_b.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 2 | bc ) -max $(echo   $min + $binrange "*" 3 | bc )  -data $3  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$3_c.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 3 | bc ) -max $(echo   $min + $binrange "*" 4 | bc )  -data $4  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif  -o input$4_d.tif 

gdalbuildvrt -separate  input$5.vrt    input$1_a.tif    input$2_b.tif   input$3_c.tif   input$4_d.tif

pkstatprofile  -f max  -i  input$5.vrt  -o  input$5_max.tif

' _ 

#this is to test what happens when two input files have the same name
