#CHELSA cold
#This is to adjust corners, set no data to -999, and clip an asc of Florida

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/cold                 

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/cold/Florida_cli\
ps

gdal_edit.py -a_ullr -180   84  180  -90 CHELSA_bio_6.tif                             

pksetmask -i CHELSA_bio_6.tif -m CHELSA_bio_6.tif -o CHELSA_bio_6.tif --msknodata -999\
99 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0 CHELSA_bio_6.tif  $OUTDIR/CHELSA_bio_6_F\
loridaClip.tif                                                                         

#gdal_translate -of AAIGrid CHELSA_bio_6_FloridaClip.tif CHELSA_bio_6_FloridaClip.asc
