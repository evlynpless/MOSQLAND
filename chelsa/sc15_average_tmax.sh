cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmax

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmax/Florida_clips

#gdal_calc.py -A CHELSA_tmax10_1_1979-2013_V1.2_land.tif -B CHELSA_tmax10_2_1979-2013_V1.2_land.tif -C CHELSA_tmax10_3_1979-2013_V1.2_land.tif -D CHELSA_tmax10_4_1979-2013_V1.2_land.tif -E CHELSA_tmax10_5_1979-2013_V1.2_land.tif -F CHELSA_tmax10_6_1979-2013_V1.2_land.tif -G CHELSA_tmax10_7_1979-2013_V1.2_land.tif -H CHELSA_tmax10_8_1979-2013_V1.2_land.tif -I CHELSA_tmax10_9_1979-2013_V1.2_land.tif -J CHELSA_tmax10_10_1979-2013_V1.2_land.tif -K CHELSA_tmax10_11_1979-2013_V1.2_land.tif -L CHELSA_tmax10_12_1979-2013_V1.2_land.tif --outfile=CHELSA_tmax_annual_average.tif --calc="A+B+C+D+E+F+G+H+I+J+K+L/12"                                        

gdal_edit.py -a_nodata -999 CHELSA_tmax_annual_average.tif

pksetmask -i CHELSA_tmax_annual_average.tif -m CHELSA_tmax_annual_average.tif -o CHELSA_tmax_annual_average.tif --msknodata -32768 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0 CHELSA_tmax_annual_average.tif  $OUTDIR/CHELSA_tmax_annual_average_FloridaClip.tif

gdal_translate -of GTiff $OUTIR/CHELSA_tmax_annual_average_FloridaClip.tif $OUTDIR/CHELSA_tmax_annual_average_FloridaClip.asc
