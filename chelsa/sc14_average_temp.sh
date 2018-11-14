#Create annual average raster from the monthly average data  (temp)                                
#Crop Florida clip                                                                             

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/temp

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/temp/Florida_clips

#gdal_calc.py -A CHELSA_temp10_1_1979-2013_V1.2_land.tif -B CHELSA_temp10_2_1979-2013_V1.2_land.tif -C CHELSA_temp10_3_1979-2013_V1.2_land.tif -D CHELSA_temp10_4_1979-2013_V1.2_land.tif -E CHELSA_temp10_5_1979-2013_V1.2_land.tif -F CHELSA_temp10_6_1979-2013_V1.2_land.tif -G CHELSA_temp10_7_1979-2013_V1.2_land.tif -H CHELSA_temp10_8_1979-2013_V1.2_land.tif -I CHELSA_temp10_9_1979-2013_V1.2_land.tif -J CHELSA_temp10_10_1979-2013_V1.2_land.tif -K CHELSA_temp10_11_1979-2013_V1.2_land.tif -L CHELSA_temp10_12_1979-2013_V1.2_land.tif --outfile=CHELSA_temp_annual_average.tif --calc="A+B+C+D+E+F+G+H+I+J+K+L/12"      

gdal_edit.py -a_nodata -999 CHELSA_temp_annual_average.tif

pksetmask -i CHELSA_temp_annual_average.tif -m CHELSA_temp_annual_average.tif -o CHELSA_temp_annual_average.tif --msknodata -32767 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0 CHELSA_temp_annual_average.tif  $OUTDIR/CHELSA_temp_annual_average_FloridaClip.tif

gdal_translate -of GTiff $OUTDIR/CHELSA_temp_annual_average_FloridaClip.tif $OUTDIR/CHELSA_temp_annual_average_FloridaClip.asc
