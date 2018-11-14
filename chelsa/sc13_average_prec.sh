#Create annual average raster from the monthly average data (prec_
#Crop Florida clip

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/prec
 
OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/prec/Florida_clips
             
#gdal_calc.py -A CHELSA_prec_1_V1.2_land.tif -B CHELSA_prec_2_V1.2_land.tif -C CHELSA_prec_3_V1.2_land.tif -D CHELSA_prec_4_V1.2_land.tif -E CHELSA_prec_5_V1.2_land.tif -F CHELSA_prec_6_V1.2_land.tif -G CHELSA_prec_7_V1.2_land.tif -H CHELSA_prec_8_V1.2_land.tif -I CHELSA_prec_9_V1.2_land.tif -J CHELSA_prec_10_V1.2_land.tif -K CHELSA_prec_11_V1.2_land.tif -L CHELSA_prec_12_V1.2_land.tif --outfile=CHELSA_prec_annual_average.tif --calc="A+B+C+D+E+F+G+H+I+J+K+L/12"

gdal_edit.py -a_nodata -999 CHELSA_prec_annual_average.tif 

pksetmask -i CHELSA_prec_annual_average.tif -m CHELSA_prec_annual_average.tif -o CHELSA_prec_annual_average.tif --msknodata -32767 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0 CHELSA_prec_annual_average.tif  $OUTDIR/CHELSA_prec_annual_average_FloridaClip.tif

gdal_translate -of GTiff CHELSA_prec_annual_average_FloridaClip.tif CHELSA_prec_annual_average_FloridaClip.asc

