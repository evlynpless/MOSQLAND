#!/bin/bash                                                                                     
#SBATCH -p scavenge                                                                             
#SBATCH -N 1                                                                                    
#SBATCH -c 1                                                                                    
#SBATCH -t 4:00:00                                                                              
#SBATCH --mail-type=ALL                                                                         
#SBATCH --mail-user=evlyn.pless@yale.edu                                                        
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc05_reclass_access.sh.%J.out                    
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc05_reclass_access.sh.%J.err                   
                                                                                                
##   sbatch  ~/scripts/MOSQLAND/sc05_reclass_access.sh 

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips


echo 1   10   20   30   40   50   60   70   80   90  100   linpos    > weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     linneg    >> weight.txt
echo 100 80   60   40   20   1    20   40   60   80  100   linpar    >> weight.txt
echo 1   25   50   75   100  150  200  250  325  400 500   exppos    >> weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     expneg    >> weight.txt
echo 100 90   80   70   60   50   40   30   20   10  1     exppar    >> weight.txt



export min=$(pkstat -nodata -999 -mm -i  accessibility_to_cities_2015_v1.0_FloridaClip.tif   | awk '{ print $2    }' )
export binrange=$(pkstat -nodata -999 -mm -i  accessibility_to_cities_2015_v1.0_FloridaClip.tif   | awk '{ print ($4-$2) / 11   }' )


cat weight.txt | xargs -n 12 -P 2 bash -c $' 

pkgetmask    -min $min                                  -max $(echo   $min + $binrange       | bc )  -data $1  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$1_a.tif 
pkgetmask    -min $(echo  $min + $binrange       | bc ) -max $(echo   $min + $binrange "*" 2 | bc )  -data $2  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$2_b.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 2 | bc ) -max $(echo   $min + $binrange "*" 3 | bc )  -data $3  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$3_c.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 3 | bc ) -max $(echo   $min + $binrange "*" 4 | bc )  -data $4  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$4_d.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 4 | bc ) -max $(echo   $min + $binrange "*" 5 | bc )  -data $5  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$5_e.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 5 | bc ) -max $(echo   $min + $binrange "*" 6 | bc )  -data $6  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$6_f.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 6 | bc ) -max $(echo   $min + $binrange "*" 7 | bc )  -data $7  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$7_g.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 7 | bc ) -max $(echo   $min + $binrange "*" 8 | bc )  -data $8  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$8_h.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 8 | bc ) -max $(echo   $min + $binrange "*" 9 | bc )  -data $9  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$9_i.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 9 | bc ) -max $(echo   $min + $binrange "*" 10 | bc )  -data $10  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$10_j.tif 
pkgetmask    -min $(echo  $min + $binrange "*" 10 | bc ) -max $(echo   $min + $binrange "*" 11 | bc )  -data $11  -nodata 0  -co COMPRESS=DEFLATE -co ZLEVEL=9 -i accessibility_to_cities_2015_v1.0_FloridaClip.tif -o input$11_k.tif 


gdalbuildvrt -separate  -input$12.vrt    input$1_a.tif    input$2_b.tif   input$3_c.tif   input$4_d.tif   input$5_e.tif   input$6_f.tif   input$7_g.tif   input$8_h.tif   input$9_i.tif   input$10_j.tif   input$11_k.tif

pkstatprofile  -f max  -i  input$12.vrt  -o  combined$12_max.tif

' _ 

